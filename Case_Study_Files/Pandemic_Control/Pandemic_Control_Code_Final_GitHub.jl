using DataFrames, Interpolations, JSON, Statistics, Plots, DifferentialEquations
using InfiniteOpt, JuMP, Juniper, Ipopt, MathOptInterface, SCIP, Gurobi, GAMS

function optimizer_type(optimizer)
    if optimizer == "Juniper-Ipopt"
        nl_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
        minlp_solver = optimizer_with_attributes(Juniper.Optimizer, "nl_solver" => nl_solver)
        model = InfiniteModel(minlp_solver)
    elseif optimizer == "Ipopt"
        model = InfiniteModel(Ipopt.Optimizer)
    elseif optimizer == "SCIP"
        model = InfiniteModel(SCIP.Optimizer)
    elseif optimizer == "Gurobi" # Best for MIPs or LPs with linear constraints, linear or quadratic objectives, and discrete integer variables. 
        model = InfiniteModel(Gurobi.Optimizer)
        # set_optimizer_attribute(model, "MIPGap", 0.01)
    else
        model = InfiniteModel(GAMS.Optimizer)
        set_optimizer_attribute(model, "solver", optimizer)
        set_optimizer_attribute(model, "reslim", 3600*9)  # Set time limit to 9 hours
    return model
    end
end
i_max = 0.02 # This is the limit of maximum fraction of infected population --> make global variable 
function solve_ode_system()
    η = 0.303
    ρ = 0.727
    N = 1e5
    ξ = 0.3
    function seir!(du, u, p, t)
        s, e, i, r = u
        u0 = 0.5  # This should be replaced with the actual control function if available
        du[1] = -(1 - u0) * ρ * s * i
        du[2] = (1 - u0) * ρ * s * i - ξ * e
        du[3] = ξ * e - η * i
        du[4] = η * i
    end
    y0 = [1 - 1/N, 1/N, 0, 0]
    tspan = (0, 200)
    t = range(tspan[1], tspan[2], length=111)
    prob = ODEProblem(seir!, y0, tspan)
    sol = solve(prob, Tsit5(), saveat=t)

    return sol
end
sol = solve_ode_system()
tspan = (0, 200)
tss = range(tspan[1], tspan[2], length=111)
s_trajectory = [sol.u[i][1] for i in 1:length(sol.t)]
e_trajectory = [sol.u[i][2] for i in 1:length(sol.t)]
i_trajectory = [sol.u[i][3] for i in 1:length(sol.t)]
r_trajectory = [sol.u[i][4] for i in 1:length(sol.t)]
u_trajectory = [0.5 for i in 1:length(sol.t)]
s_inter = linear_interpolation(tss, s_trajectory)
e_inter = linear_interpolation(tss, e_trajectory)
i_inter = linear_interpolation(tss, i_trajectory)
r_inter = linear_interpolation(tss, r_trajectory)
u_inter = linear_interpolation(tss, u_trajectory)
s_s = t -> s_inter(t)
e_s = t -> e_inter(t)
i_s = t -> i_inter(t)
r_s = t -> r_inter(t)
u_s = t -> u_inter(t)
Func = [s_s, e_s, i_s, r_s, u_s]

function base_model(optimizer, func)
    model = optimizer_type(optimizer)

    # Set the SEIR parameters
    η = 0.303
    ρ = 0.727
    N = 1e5
    ξ = 0.3
    
    t_start = 0
    t_end = 200
    extra_ts = [0.001, 0.002, 0.004, 0.008, 0.02, 0.04, 0.08, 0.2, 0.4, 0.8]
    e0 = 1 / N
    i0 = 0
    r0 = 0
    s0 = 1 - e0;

    # Add the infinite parameter
    @infinite_parameter(model, t ∈ [t_start, t_end], num_supports = 101) # t ∈ [0, 200] with 101 supports
    add_supports(t, extra_ts)

    # Add the variables
    @variable(model, 0 ≤ s ≤ 1, Infinite(t), start = func[1]) # 0 ≤ s(t) ≤ 1
    @variable(model, 0 ≤ e ≤ 1, Infinite(t), start = func[2]) # 0 ≤ e(t) ≤ 1
    @variable(model, 0 ≤ i ≤ 1, Infinite(t), start = func[3]) # 0 ≤ i(t) ≤ 1
    @variable(model, 0 ≤ r ≤ 1, Infinite(t), start = func[4]) # 0 ≤ r(t) ≤ 1
    @variable(model, 0 ≤ u ≤ 0.8, Infinite(t), start = func[5]) # 0 ≤ u(t) ≤ 0.8 with guess of 0.5

    # Add the objective
    @objective(model, Min, ∫(u, t, t_start, t_end)) # minimize ∫_{t ∈ [0, 200]} u(t) dt

    # Add the initial conditions
    fix(s(0), s0, force = true) # add s(0)
    fix(e(0), e0, force = true) # add e(0)
    fix(i(0), i0, force = true) # add i(0)
    fix(r(0), r0, force = true) # add r(0)

    # Define the SEIR equations
    @constraint(model, ∂(s,t) == -(1-u)*ρ*s*i) # d/dt[s(t)] = -(1 - u(t))ρs(t)i(t)
    @constraint(model, ∂(e,t) == (1-u)*ρ*s*i-ξ*e) # d/dt[e(t)] = (1 - u(t))ρs(t)i(t) - ξe(t)
    @constraint(model, ∂(i,t) == ξ*e-η*i) # d/dt[i(t)] = ξe(t) - ηi(t)
    @constraint(model, ∂(r,t) == η*i) # d/dt[r(t)] = ηi(t)

    ########################################################################

    # Define the left-hand-side of the constraint we want to enforce of the form h(t) <= i-imax
    @expression(model, h, i - i_max)
    return model
end
function execute_model(model)
    # set_silent(model)
    optimize!(model)
    println("Optimal Objective:  ", objective_value(model))
    println("Time for solver to compute:    ", round(JuMP.solve_time(model), digits = 5), " s")
    println(termination_status(model))
    return 
end
function hard_constraint()
    func = [0,0,0,0,0.5]
    model = base_model("Ipopt", func)
    @constraint(model, model[:h] ≤ 0)
    execute_model(model)
    return model
end
function interpolate_data(method, α = 1)
    if method == "hard"
        model = hard_constraint()
        ts = value(model[:t])
        s_vals = value(model[:s])
        e_vals = value(model[:e])
        i_vals = value(model[:i])
        r_vals = value(model[:r])
        u_vals = value(model[:u])
        s_interp = linear_interpolation(ts, s_vals)
        e_interp = linear_interpolation(ts, e_vals)
        i_interp = linear_interpolation(ts, i_vals)
        r_interp = linear_interpolation(ts, r_vals)
        u_interp = linear_interpolation(ts, u_vals)
        s_start = t -> s_interp(t)
        e_start = t -> e_interp(t)
        i_start = t -> i_interp(t)
        r_start = t -> r_interp(t)
        u_start = t -> u_interp(t)
        func = [s_start, e_start, i_start, r_start, u_start]
    elseif method == "CVaR"
        model = CVaR(α)
        #interpolate:
        ts = value(model[:t])
        ss = value(model[:s])
        es = value(model[:e])
        is = value(model[:i])
        rs = value(model[:r])
        us = value(model[:u])
        ϕs = value(model[:ϕ])
        s_interp = linear_interpolation(ts,ss)
        s_start = ts -> s_interp(ts)
        e_interp = linear_interpolation(ts,es)
        e_start = ts -> e_interp(ts)
        i_interp = linear_interpolation(ts,is)
        i_start = ts -> i_interp(ts)
        r_interp = linear_interpolation(ts,rs)
        r_start = ts -> r_interp(ts)
        u_interp = linear_interpolation(ts,us)
        u_start = ts -> u_interp(ts)
        ϕ_interp = linear_interpolation(ts,ϕs)
        ϕ_start = ts -> ϕ_interp(ts)
        func = [s_start, e_start, i_start, r_start, u_start, ϕ_start]
    end
    return model, func
end
function big_M(α) #This is an exact formulation
    # requires a mixed-integer solver (i.e., Juniper or Gurobi depending on if it's linear or not)
    func = [0,0,0,0,0.5]
    # func = interpolate_data("hard")[2] # initialize the parameters with the solution of the hard constraint
    model = base_model("Juniper-Ipopt", func)
    M = 1 #pick a value that is as big as `h(t)` can be 
    @variable(model, q, Bin, Infinite(model[:t]))
    @constraint(model, model[:h] <= (1 - q) * M) 
    @constraint(model, 𝔼(q, model[:t]) ≥ α)
    execute_model(model)
    return model
end
function CVaR(α)  #This is a continuous approximation
    # The SigVaR paper describes this: https://arxiv.org/abs/2004.02402 (page 17)
    func = [0,0,0,0,0.5]
    # func = interpolate_data("hard")[2]
    model = base_model("Ipopt", func)
    @variable(model, λ ≥ 0) # the optimal value is needed for initializing the SigVaR approach
    @variable(model, ϕ ≥ 0, Infinite(model[:t]))
    @constraint(model, ϕ ≥ model[:h] - λ) 
    @constraint(model, 𝔼(ϕ, model[:t]) ≤ -λ * (1 - α))
    execute_model(model)
    return model
end
function Create_SigVaR(μ, τ, α, func) # Enforce the continuous SigVaR approximation
    model = base_model("Ipopt", func)
    @variable(model, ϕ ≥ 0, Infinite(model[:t]), start = func[6])
    @constraint(model, ϕ ≥ (2*(1 + μ)) / (μ + exp(-τ * model[:h])) - 1)
    @constraint(model, 𝔼(ϕ, model[:t]) ≤ 1 - α)
    execute_model(model)
    return model
end
#create storage arrays for all the optimal values through each iteration
ϕ_vals = Vector{Any}()
s_vals = Vector{Any}()
e_vals = Vector{Any}()
i_vals = Vector{Any}()
r_vals = Vector{Any}()
u_vals = Vector{Any}()
h_vals = Vector{Any}()
μ_vals = Vector{Any}()
τ_vals = Vector{Any}()
iters = Vector{Any}()
obj_vals = Vector{Any}()
status_vals = Vector{Any}()
λ_vals = Vector{Any}()
runtime_vals = Vector{Any}() 

function SigVaR_algorithm(α, num_iterations)
    # Described in https://arxiv.org/abs/2004.02402 (page 16)
    # interpolate:
    cvar = interpolate_data("CVaR", α)
    model = cvar[1] # initialize the parameters with the solution of the CVaR
    func = cvar[2]
    lastobj = objective_value(model)
    println("CVaR Objective = ", lastobj)
    λ_val = value(model[:λ])
    Γ = -1/(λ_val - i_max)
    push!(λ_vals, λ_val)

    μ = 1.55026018
    τ = ((μ + 1)/2)*Γ

    last_feasible_model = model

    for ł in 1:num_iterations
        SigVaR_model = Create_SigVaR(μ, τ, α, func)
        status = termination_status(SigVaR_model)
        if status != MathOptInterface.LOCALLY_SOLVED
            break
        end        
        last_feasible_model = SigVaR_model
        newobj = objective_value(SigVaR_model)
        offset = abs(newobj-lastobj)
        ts = value(SigVaR_model[:t])

        lastobj = newobj
        #store optimal values in empty arrays:
        push!(status_vals, status)
        push!(obj_vals, objective_value(SigVaR_model))
        push!(iters, ł)
        push!(ϕ_vals, value(SigVaR_model[:ϕ]))
        push!(s_vals, value(SigVaR_model[:s]))
        push!(e_vals, value(SigVaR_model[:e]))
        push!(i_vals, value(SigVaR_model[:i]))
        push!(r_vals, value(SigVaR_model[:r]))
        push!(u_vals, value(SigVaR_model[:u]))
        push!(h_vals, value(SigVaR_model[:h]))
        push!(μ_vals, μ)
        push!(τ_vals, τ)
        push!(runtime_vals, JuMP.solve_time(SigVaR_model))
        if offset <= 0.001
            break
        end
        s_interp = linear_interpolation(ts,s_vals[ł])
        s_start = ts -> s_interp(ts)
        e_interp = linear_interpolation(ts,e_vals[ł])
        e_start = ts -> e_interp(ts)
        i_interp = linear_interpolation(ts,i_vals[ł])
        i_start = ts -> i_interp(ts)
        r_interp = linear_interpolation(ts,r_vals[ł])
        r_start = ts -> r_interp(ts)
        u_interp = linear_interpolation(ts,u_vals[ł])
        u_start = ts -> u_interp(ts)
        ϕ_interp = linear_interpolation(ts,ϕ_vals[ł])
        ϕ_start = ts -> ϕ_interp(ts)

        func = [s_start, e_start, i_start, r_start, u_start, ϕ_start]
        
        println("μ = ", μ)
        println("τ = ", τ)
        w = 2
        μ *= w
        τ = ((μ + 1)/2)*Γ
        # plot_iterations(ts, value(SigVaR_model[:i]), value(SigVaR_model[:u]))
        println("iteration # ", ł)
        println("offset from last iteration = ", offset)
    end
    return last_feasible_model
end
function complimentarity(α, ϵ, func) 
    model = base_model("ipopt", func)
    @variable(model, 0 ≤ y1 ≤ 1, Infinite(model[:t]), start = 0.5)
    @variable(model, 0 ≤ y2 ≤ 1, Infinite(model[:t]), start = 0.5)
    @constraint(model, (1/2)*(y1 + y2 - √((y1-y2)^2 + ϵ^2)) ≤ ϵ) #max upper bound
    @constraint(model, (1/2)*(y1 + y2 - √((y1-y2)^2 + ϵ^2)) ≥ -ϵ) #max lower bound
    @constraint(model, y1 + y2 ≤ 1 + ϵ) #linear upper bound
    @constraint(model, y1 + y2 ≥ 1 - ϵ) #linear lower bound
    tol = 1e-9
    @constraint(model, model[:h] ≤ (1 - y1) * (1 - i_max + 5*tol)) 
    @constraint(model, 𝔼(y1, model[:t]) ≥ α)
    execute_model(model)
    return model
end
function plot_iterations(ts, i, u)
    p0 = plot(
        legend=:best,
        xlabel = "Time (Days)",
        ylabel = "i(t)",
        gridalpha = 0.2,
        legendfontsize = 10,
    )
    plot!(ts, i, label = "Complimentarity", linecolor = :red, linewidth = 4)

    # Plot u(t):
    p1 = plot(
        legend=:best,
        xlabel = "Time (Days)",
        ylabel = "u(t)",
        gridalpha = 0.2,
        legendfontsize = 10,
    )
    plot!(ts, u, label = "Complimentarity", linecolor = :orange, linewidth = 4)
    ylims!(-0.04, 1)

    plt=plot(p0, p1, layout = (2, 1), size=(800, 600), dpi = 400)
    display(plt)
end
comp_obj_vals = Vector{Any}()
comp_status_vals = Vector{Any}()
comp_runtime_vals = Vector{Any}() 

function compliment_iterate(α)
    ϵs = [1.00000000e+00, 9.55000000e-01, 9.10000000e-01, 8.65000000e-01,
    8.20000000e-01, 7.75000000e-01, 7.30000000e-01, 6.85000000e-01,
    6.40000000e-01, 5.95000000e-01, 5.50000000e-01, 5.05000000e-01,
    4.60000000e-01, 4.15000000e-01, 3.70000000e-01, 3.25000000e-01,
    2.80000000e-01, 2.35000000e-01, 1.90000000e-01, 1.45000000e-01,
    1.00000000e-01, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 4.46683592e-02, 1.99526231e-02, 8.91250938e-03,
    3.98107171e-03, 1.77827941e-03, 7.94328235e-04, 3.54813389e-04,
    1.58489319e-04, 7.07945784e-05, 3.16227766e-05, 1.41253754e-05,
    6.30957344e-06, 2.81838293e-06, 1.25892541e-06, 5.62341325e-07,
    2.51188643e-07, 1.12201845e-07, 5.01187234e-08, 2.23872114e-08,
    1.00000000e-08]
    # Use hard constraint as initialization
    hard = interpolate_data("hard")
    initial_solution = hard[1]
    func = hard[2]
    
    last_feasible_model = initial_solution
    for ϵ in ϵs
        model = complimentarity(α, ϵ, func)
        status = termination_status(model)
        if status != MathOptInterface.LOCALLY_SOLVED
            println("Infeasible solution for ϵ = ", ϵ)
            break
        end
        last_feasible_model = model
        ts = value(model[:t])
        initial_values_s = value(model[:s])
        initial_values_e = value(model[:e])
        initial_values_i = value(model[:i])
        initial_values_r = value(model[:r])
        initial_values_u = value(model[:u])

        s_interp = linear_interpolation(ts,initial_values_s)
        s_start = ts -> s_interp(ts)
        e_interp = linear_interpolation(ts,initial_values_e)
        e_start = ts -> e_interp(ts)
        i_interp = linear_interpolation(ts,initial_values_i)
        i_start = ts -> i_interp(ts)
        r_interp = linear_interpolation(ts,initial_values_r)
        r_start = ts -> r_interp(ts)
        u_interp = linear_interpolation(ts,initial_values_u)
        u_start = ts -> u_interp(ts)
        func = [s_start, e_start, i_start, r_start, u_start]

        push!(comp_status_vals, status)
        push!(comp_obj_vals, objective_value(model))
        push!(comp_runtime_vals, solve_time(model))
        
        # plot_iterations(ts, initial_values_i, initial_values_u)
        println("epsilon = ", ϵ)
        println("------------------- Complimentarity Method -------------------")
    end
    return last_feasible_model
end 
function runmodel(method, α = 1) 
    if method == "Hard"
        model = hard_constraint()
    elseif method == "bigM"
        model = big_M(α)
    elseif method == "CVaR"
        model = CVaR(α)
    elseif method == "SigVaR"
        num_iterations = 15
        model = SigVaR_algorithm(α, num_iterations)
    end
    return model
end
function clear_global_vectors()
    global ϕ_vals = Vector{Any}()
    global s_vals = Vector{Any}()
    global e_vals = Vector{Any}()
    global i_vals = Vector{Any}()
    global r_vals = Vector{Any}()
    global u_vals = Vector{Any}()
    global h_vals = Vector{Any}()
    global μ_vals = Vector{Any}()
    global τ_vals = Vector{Any}()
    global iters = Vector{Any}()
    global obj_vals = Vector{Any}()
    global status_vals = Vector{Any}()
    global λ_vals = Vector{Any}()
    global runtime_vals = Vector{Any}()
    global comp_obj_vals = Vector{Any}()
    global comp_status_vals = Vector{Any}()
    global comp_runtime_vals = Vector{Any}()
end
function master(α, method)
    if method == "Hard"
        #---------------- Hard Constraint Method -------------------#
        model = runmodel("Hard")
        ts = value(model[:t])
        s_opt = value(model[:s])
        e_opt = value(model[:e])
        i_opt = value(model[:i])
        r_opt = value(model[:r])
        u_opt = value(model[:u])
        h_opt = value(model[:h])
        solve_time = JuMP.solve_time(model)

        α_min = mean(h_opt .< 1e-5)
        println("The minimum α value = ", α_min)

        data = Dict(
            # Save Hard Constraint Results:
            "s_opt" => s_opt, "e_opt" => e_opt,
            "i_opt" => i_opt*100, "r_opt" => r_opt, 
            "h_opt" => h_opt, "u_opt" => u_opt, "α_min" => α_min, 
            "ts" => ts, "objective value" => objective_value(model),
            "solver time" => round(solve_time, digits = 5),
            "termination status" => termination_status(model),
            "solver name" => solver_name(model)
        )
        json_data = JSON.json(data)
        file_path = "data_folder/hard.json"
        open(file_path, "w") do file
            println(file, json_data)
        end

    elseif method == "CVaR"
        #--------------------- CVaR Method --------------------------#
        model = runmodel("CVaR", α)
        ts = value(model[:t])
        s_opt = value(model[:s])
        e_opt = value(model[:e])
        i_opt = value(model[:i])
        r_opt = value(model[:r])
        u_opt = value(model[:u])
        h_opt = value(model[:h])
        ϕ_opt = value(model[:ϕ])
        solve_time = JuMP.solve_time(model)

        data = Dict(
            "s_opt" => s_opt, "e_opt" => e_opt,
            "i_opt" => i_opt*100, "r_opt" => r_opt, 
            "h_opt" => h_opt, "u_opt" => u_opt, 
            "λ values" => λ_vals, "ts" => ts,
            "solver time" => round(solve_time, digits = 5),
            "objective value" => objective_value(model),
            "termination status" => termination_status(model),
            "solver name" => solver_name(model)
        )
        json_data = JSON.json(data)
        file_path = "data_folder/CVaR_alpha_$α.json"
        open(file_path, "w") do file
            println(file, json_data)
        end

    elseif method == "SigVaR"
        #--------------------- SigVaR Method --------------------------#
        model = runmodel("SigVaR", α)
        ts = value(model[:t])
        s_opt = value(model[:s])
        e_opt = value(model[:e])
        i_opt = value(model[:i])
        r_opt = value(model[:r])
        u_opt = value(model[:u])
        ϕ_opt = value(model[:ϕ])

        data = Dict(
            # Save SigVaR Results:
            "τ values" => τ_vals, "μ values" => μ_vals, 
            "ϕ values" => ϕ_vals, "h values" => h_vals,
            "i values" => i_vals*100, "u values" => u_vals,
            "iterations" => iters, "termination status" => status_vals, 
            "i_opt" => i_opt*100, 
            "ϕ_opt" => ϕ_opt, "u_opt" => u_opt, 
            "ts" => ts, "solver name" => solver_name(model),
            "solver times" => round.(runtime_vals, digits = 5),
            "objective value" => obj_vals
        )
        json_data = JSON.json(data)
        file_path = "data_folder/SigVaR_alpha_$α.json"
        open(file_path, "w") do file
            println(file, json_data)
        end

    elseif method == "bigM"
        #--------------------- BigM Method --------------------------#
        model = runmodel("bigM", α)
        ts = value(model[:t])
        s_opt = value(model[:s])
        e_opt = value(model[:e])
        i_opt = value(model[:i])
        r_opt = value(model[:r])
        u_opt = value(model[:u])
        h_opt = value(model[:h])
        q_opt = value(model[:q])
        solve_time = JuMP.solve_time(model)

        data = Dict(
            # Save BigM Results:
            "s_opt" => s_opt, "e_opt" => e_opt,
            "i_opt" => i_opt*100, "r_opt" => r_opt, 
            "h_opt" => h_opt, "u_opt" => u_opt, 
            "q_opt" => q_opt, "ts" => ts, 
            "solver time" => round(solve_time, digits = 5),
            "solver name" => solver_name(model), "termination status" => termination_status(model),
            "objective value" => objective_value(model)
        )
        json_data = JSON.json(data)
        file_path = "data_folder/bigM_alpha_$α.json"
        open(file_path, "w") do file
            println(file, json_data)
        end

    elseif method == "complimentarity"
        # -------------- Complimentarity Method -------------- #
        model = compliment_iterate(α)
        ts = value(model[:t])
        s_opt = value(model[:s])
        e_opt = value(model[:e])
        i_opt = value(model[:i])
        r_opt = value(model[:r])
        u_opt = value(model[:u])
        h_opt = value(model[:h])
        y1_opt_compliment = value(model[:y1])
        y2_opt_compliment = value(model[:y2])
        solve_time = JuMP.solve_time(model)

        data = Dict(
            # Save Complimentarity Results:
            "s_opt" => s_opt, "e_opt" => e_opt,
            "i_opt" => i_opt*100, "r_opt" => r_opt, 
            "h_opt" => h_opt, "u_opt" => u_opt, 
            "ts" => ts, 
            "y1_opt_compliment" => y1_opt_compliment, 
            "y2_opt_compliment" => y2_opt_compliment, 
            "termination status" => comp_status_vals, 
            "solver times" => round.(comp_runtime_vals, digits = 5),
            "solver name" => solver_name(model),
            "objective values" => comp_obj_vals,
            "objective value" => comp_obj_vals[end],
            )

        json_data = JSON.json(data)
        file_path = "data_folder/compliment_alpha_$α.json"
        open(file_path, "w") do file
            println(file, json_data)
        end
    end
    clear_global_vectors()
end
αs = [0.85, 0.9, 0.95, 0.96, 0.97, 0.99] # Chosen α values

for α in αs
    master(α, "Hard")
    master(α, "CVaR")
    master(α, "SigVaR")
    master(α, "complimentarity")
end

for α in αs # Run the bigM method seperately as it can take magnitudes longer to converge
    master(α, "bigM")
end






