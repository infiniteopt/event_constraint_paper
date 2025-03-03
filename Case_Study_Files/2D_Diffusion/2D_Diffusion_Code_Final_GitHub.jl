using DataFrames, Interpolations, JSON, Statistics, Plots
using InfiniteOpt, JuMP, Juniper, Ipopt, MathOptInterface, Gurobi, GAMS

function optimizer_type(optimizer)
    if optimizer == "Juniper-Ipopt" # Best for Mixed-Integer Nonlinear Programming (MINLPs)
        nl_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
        minlp_solver = optimizer_with_attributes(Juniper.Optimizer, "nl_solver" => nl_solver)
        model = InfiniteModel(minlp_solver)
    elseif optimizer == "Ipopt" # Best for Nonlinear Programming (NLP) with continuous variables
        model = InfiniteModel(Ipopt.Optimizer)

    elseif optimizer == "Gurobi" # Best for MIPs or LPs with linear constraints, linear or quadratic objectives, and discrete integer variables. 
        model = InfiniteModel(Gurobi.Optimizer)
        # set_optimizer_attribute(model, "MIPGap", 0.01)
    else
        model = InfiniteModel(GAMS.Optimizer)
        set_optimizer_attribute(model, "solver", optimizer)
        set_optimizer_attribute(model, "reslim", 3600*15)  # Set time limit to 15 hours
    return model
    end
end
Tmax = 1.1
Tsp = 1.0
num_1d_grid_pts = 62
D = 0.05
xs = collect(LinRange(-1, 1, num_1d_grid_pts))
num_1d_heaters = 6 # choose such that (num_1d_grid_pts - num_1d_heaters) / (num_1d_heaters + 1) is an integer
# Determine heater placements
spacing = Int((num_1d_grid_pts - num_1d_heaters) / (num_1d_heaters + 1)) + 1
heater_positions_1d = xs[spacing:spacing:num_1d_grid_pts] 
function base_model(optimizer, func, linear = true)
    # Initialize the model (Gurobi if linear, Ipopt if not) and infinite parameter
    model = optimizer_type(optimizer)
    @infinite_parameter(
        model,
        x[1:2] ∈ [-1, 1],
        independent = true,
        supports = xs,
        derivative_method = FiniteDifference(Central())
    )

    # Create the variables
    @variable(model, 0 <= T <= 5.5, Infinite(x), start = func[1])
    @variable(model, u == 0, Infinite(x), start = func[2]) # heater 
    # Make it so heater only heats at certain points
    for x1 in heater_positions_1d
        for x2 in heater_positions_1d
            u_pt = u([x1, x2])
            unfix(u_pt)
            set_lower_bound(u_pt, 0.0)
            set_upper_bound(u_pt, linear ? 50 : 50^2)
            set_start_value(u_pt, linear ? 25 : 720)
        end
    end

    # Set the objective function
    @objective(model, Min, ∫(∫((T - Tsp)^2, x[1]), x[2]))

    heating = linear ? u : sqrt(u + 0.00001)
    
    @constraint(model, ODE, D * (@∂(T, x[1]^2) + @∂(T, x[2]^2)) + heating - 0.1 == 0)
    @constraint(model, T([-1, x[2]]) == 0)
    @constraint(model, T([1, x[2]]) == 0)
    @constraint(model, T([x[1], -1]) == 0)
    @constraint(model, T([x[1], 1]) == 0)

    # Define the LHS of the spatial constraint we seek to relax
    @expression(model, h, T - Tmax) #dont want to go past 1.1 Temp max
    return model
end
function execute_model(model)
    # set_silent(model)
    optimize!(model)
    println("Optimal Objective:  ", objective_value(model))
    println("Time for solver to compute:    ", round(JuMP.solve_time(model), digits = 5), " s")
    println("termination status: ", termination_status(model))
    return 
end
function hard_constraint(linear)
    func = [0,0,0]
    model = base_model("Ipopt", func, linear)
    @constraint(model, hard, model[:h] ≤ 0)
    execute_model(model)
    println("------------------ Hard Constraint Method -------------------")
    return model
end
function interpolate_data(linear, method, α = 1)
    if method == "hard"
        model = hard_constraint(linear)
        xs1 = value.(model[:x][1])
        xs2 = value.(model[:x][2])
        Ts = value.(model[:T])
        us = value.(model[:u])
    
        T_interp = interpolate((xs1,xs2), Ts, Gridded(Linear()))
        T_start = x -> T_interp(x[1], x[2])
        u_interp = interpolate((xs1,xs2), us, Gridded(Linear()))
        u_start = x -> u_interp(x[1], x[2])
    
        func = [T_start, u_start]
    elseif method == "CVaR"
        model = CVaR(α, linear)
        #interpolate:
        xs1 = value.(model[:x][1])
        xs2 = value.(model[:x][2])
        Ts = value.(model[:T])
        us = value.(model[:u])
        ϕs = value(model[:ϕ])
    
        T_interp = interpolate((xs1,xs2), Ts, Gridded(Linear()))
        T_start = x -> T_interp(x[1], x[2])
        u_interp = interpolate((xs1,xs2), us, Gridded(Linear()))
        u_start = x -> u_interp(x[1], x[2])
        ϕ_interp = interpolate((xs1,xs2), ϕs, Gridded(Linear()))
        ϕ_start = x -> ϕ_interp(x[1], x[2])
    
        func = [T_start, u_start, ϕ_start]
    end
    return model, func
end
function big_M(α, linear, solver = "Gurobi") #This is an exact formulation
    # requires a mixed-integer solver (i.e., Juniper or Gurobi depending on if it's linear or not)
    func = [0,0,0]
    # func = interpolate_data(linear, "hard")[2]
    model = base_model(linear ? "Gurobi" : solver, func, linear) # we want to solve the linear problem with Gurobi if linear is true, otherwise with dicopt

    # One-sided big-M
    @variable(model, q, Bin, Infinite(model[:x]))
    @constraint(model, model[:h] <= q * 5.5)
    prob = 𝔼(𝔼(1 - q, model[:x][1]), model[:x][2])
    @constraint(model, prob ≥ α)

    execute_model(model)
    println("------------------- Big-M Method -------------------")

    return model
end
function CVaR(α, linear)  #This is a continuous approximation
    # OPTION 3 (Enforce the continuous CVaR approximation)
    # The SigVaR paper describes this: https://arxiv.org/abs/2004.02402 (page 17)
    # CVaR approximation is as follows:
    func = [0,0,0]
    # func = interpolate_data(linear, "hard")[2]

    model = base_model(linear ? "Gurobi" : "Ipopt", func, linear) # we want to solve the linear problem with Gurobi if linear is true, otherwise with dicopt
    @variable(model, λ ≥ 0) # the optimal value is needed for initializing the SigVaR approach
    @variable(model, ϕ ≥ 0, Infinite(model[:x]))
    @constraint(model, ϕ ≥ model[:h] - λ) 
    @constraint(model, 𝔼(𝔼(ϕ, model[:x][1]),model[:x][2]) ≤ -λ * (1 - α)) #Eq 2.11 in your paper #Eq 36d in his paper
    execute_model(model)
    λ_val = value(model[:λ])
    Γ = -1/(λ_val - Tmax)
    println("------------------- CVaR Method -------------------")
    println("λ = ", λ_val)
    println("Γ = ", Γ)
    return model
end
function Create_SigVaR(μ, τ, α, func, linear) # Enforce the continuous SigVaR approximation
    model = base_model("Ipopt", func, linear)
    @variable(model, ϕ ≥ 0, Infinite(model[:x]), start = func[3])
    @constraint(model, ϕ ≥ (2*(1 + μ)) / (μ + exp(-τ * model[:h])) - 1)
    @constraint(model, 𝔼(𝔼(ϕ, model[:x][1]), model[:x][2]) ≤ 1 - α)
    execute_model(model)
    return model
end
#create storage arrays for all the optimal values through each iteration
ϕ_vals = Vector{Any}()
T_vals = Vector{Any}()
u_vals = Vector{Any}()
h_vals = Vector{Any}()
μ_vals = Vector{Any}()
τ_vals = Vector{Any}()
iters = Vector{Any}()
obj_vals = Vector{Any}()
status_vals = Vector{Any}()
λ_vals = Vector{Any}()
runtime_vals = Vector{Any}() 
function SigVaR_algorithm(α, num_iterations, linear)
    # Described in https://arxiv.org/abs/2004.02402 (page 16)
    # Formulation is equation 35
    # Note that our version of α is 1-α in their paper
    # Here `h` is the same thing as the `z` used their paper
    cvar = interpolate_data(linear, "CVaR", α)
    model = cvar[1] # initialize the parameters with the solution of the CVaR
    func = cvar[2]    
    lastobj = objective_value(model)
    println("CVaR Objective = ", lastobj)
    λ_val = value(model[:λ])
    Γ = -1/(λ_val - Tmax)
    push!(λ_vals, λ_val) 

    μ = 15.5026018
    τ = ((μ + 1)/2)*Γ
    println("μ = ", μ)
    println("τ = ", τ)  
    last_feasible_model = model

    for ł in 1:num_iterations
        SigVaR_model = Create_SigVaR(μ, τ, α, func, linear)
        status = termination_status(SigVaR_model)
        if status != MathOptInterface.LOCALLY_SOLVED
            break
        end
        last_feasible_model = SigVaR_model
        newobj = objective_value(SigVaR_model)
        offset = abs(newobj-lastobj)
        lastobj = newobj
        #store optimal values in empty arrays:
        push!(status_vals, status)
        push!(obj_vals, objective_value(SigVaR_model))
        push!(iters, ł)
        push!(ϕ_vals, value.(SigVaR_model[:ϕ]))
        push!(T_vals, value.(SigVaR_model[:T]))
        push!(u_vals, value.(SigVaR_model[:u]))
        push!(h_vals, value.(SigVaR_model[:h]))
        push!(μ_vals, μ)
        push!(τ_vals, τ)
        push!(runtime_vals, JuMP.solve_time(SigVaR_model))
        # push!(num_iterations_vals, result_count(SigVaR_model, :iteration))
        
        if offset <= 0.003
            break
        end
        xs1 = value.(SigVaR_model[:x][1])
        xs2 = value.(SigVaR_model[:x][2])

        T_interp = interpolate((xs1,xs2), T_vals[ł], Gridded(Linear()))
        T_start = x -> T_interp(x[1],x[2])
        u_interp = interpolate((xs1,xs2), u_vals[ł], Gridded(Linear()))
        u_start = x -> u_interp(x[1],x[2])
        ϕ_interp = interpolate((xs1,xs2), ϕ_vals[ł], Gridded(Linear()))
        ϕ_start = x -> ϕ_interp(x[1],x[2])

        func = [T_start, u_start, ϕ_start]
        println("------------------- SigVaR Method -------------------")
        println("μ = ", μ)
        println("τ = ", τ)
        w = 1.2
        μ *= w
        τ = ((μ + 1)/2)*Γ

        println("iteration # ", ł)
        println("offset from last iteration = ", offset)
    end
    return last_feasible_model
end
function complimentarity(α, solver, ϵ, func, linear) 
    model = base_model(solver, func, linear)
    @variable(model, 0 ≤ y1 ≤ 1, Infinite(model[:x]), start = 0.5)
    @variable(model, 0 ≤ y2 ≤ 1, Infinite(model[:x]), start = 0.5)
    @constraint(model, (1/2)*(y1 + y2 - √((y1-y2)^2 + ϵ^2)) ≤ ϵ) #max upper bound
    @constraint(model, (1/2)*(y1 + y2 - √((y1-y2)^2 + ϵ^2)) ≥ -ϵ) #max lower bound
    @constraint(model, y1 + y2 ≤ 1 + ϵ) #linear upper bound
    @constraint(model, y1 + y2 ≥ 1 - ϵ) #linear lower bound
    tol = 1e-9
    @constraint(model, model[:h] ≤ (1 - y1) * (1 - Tmax + 5*tol)) 
    @constraint(model, 𝔼(𝔼(y1, model[:x][1]), model[:x][2]) ≥ α)
    execute_model(model)
    return model
end
comp_obj_vals = Vector{Any}()
comp_status_vals = Vector{Any}()
comp_runtime_vals = Vector{Any}() 
ϵ_vals = Vector{Any}()  
function compliment_iterate(α, solver, linear)
    ϵs = [1.00000000e+00, 9.55000000e-01, 9.10000000e-01, 8.65000000e-01,
    8.20000000e-01, 7.75000000e-01, 7.30000000e-01, 6.85000000e-01,
    6.40000000e-01, 5.95000000e-01, 5.50000000e-01, 5.05000000e-01,
    4.60000000e-01, 4.15000000e-01, 3.70000000e-01, 3.25000000e-01,
    2.80000000e-01, 2.35000000e-01, 1.90000000e-01, 1.45000000e-01,
    1.00000000e-01, 4.46683592e-02, 1.99526231e-02, 8.91250938e-03,
    3.98107171e-03, 1.77827941e-03, 7.94328235e-04, 3.54813389e-04,
    1.58489319e-04, 7.07945784e-05, 3.16227766e-05, 1.41253754e-05,
    6.30957344e-06, 2.81838293e-06, 1.25892541e-06, 5.62341325e-07,
    2.51188643e-07, 1.12201845e-07, 5.01187234e-08, 2.23872114e-08,
    1.00000000e-08]

    # Use hard constraint as initialization
    hard = interpolate_data(linear, "hard")
    initial_solution = hard[1]
    func = hard[2]

    last_feasible_model = initial_solution
    for ϵ in ϵs
        model = complimentarity(α, solver, ϵ, func, linear)
        status = termination_status(model)
        if status != MathOptInterface.LOCALLY_SOLVED
            break
        end
        last_feasible_model = model
        xs1 = value.(model[:x][1])
        xs2 = value.(model[:x][2])
        T_vals_iter = value.(model[:T])
        u_vals_iter = value.(model[:u])

        T_interp = interpolate((xs1, xs2), T_vals_iter, Gridded(Linear()))
        T_start = x -> T_interp(x[1], x[2])
        u_interp = interpolate((xs1, xs2), u_vals_iter, Gridded(Linear()))
        u_start = x -> u_interp(x[1], x[2])
        func = [T_start, u_start]

        push!(comp_status_vals, status)
        push!(comp_obj_vals, objective_value(model))
        push!(comp_runtime_vals, JuMP.solve_time(model))
        push!(ϵ_vals, ϵ)

        println("epsilon = ", ϵ)
        println("------------------- Complimentarity Method -------------------")
    end
    return last_feasible_model
end
function runmodel(method, linear, α = 1)
    if method == "Hard"
        model = hard_constraint(linear)
    elseif method == "bigM"
        model = big_M(α, linear)
    elseif method == "bigM_juniper"
        model = big_M(α, linear, "Juniper-Ipopt")
    elseif method == "bigM_dicopt"
        model = big_M(α, linear, "dicopt")
    elseif method == "bigM_scip"
        model = big_M(α, linear, "scip")
    elseif method == "bigM_baron"
        model = big_M(α, linear, "baron")
    elseif method == "CVaR"
        model = CVaR(α, linear)
    elseif method == "SigVaR"
        num_iterations = 30
        model = SigVaR_algorithm(α, num_iterations, linear)
    elseif method == "complimentarity_ipopt"
        model = compliment_iterate(α, "Ipopt", linear)
    elseif method == "complimentarity_conopt4"
        model = compliment_iterate(α, "conopt4", linear)
    end
    return model
end
function clear_global_vectors()
    global ϕ_vals = Vector{Any}()
    global T_vals = Vector{Any}()
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
function master(α, method, linear)
    if method == "Hard"
        #---------------- Hard Constraint Method -------------------#
        model = runmodel(method, linear)
        Xs = value.(model[:x])
        T_opt = vec(value(model[:T]))
        h_opt = vec(value(model[:h]))
        u_opt = vec(value(model[:u]))
        α_min = mean(h_opt .< 1e-5)
        solvetime = JuMP.solve_time(model)
        optimizer = solver_name(model)
        println(optimizer)

        data = Dict(
            # Save Hard Constraint Results:
            "T_opt" => T_opt, "h_opt" => h_opt, 
            "u_opt" => u_opt, "α_min" => α_min,  "xs" => Xs,
            "solver time" => round(solvetime, digits = 5),
            "solver name" => optimizer, "termination_status" => termination_status(model),
            "objective value" => objective_value(model),
        )
        json_data = JSON.json(data)
        file_path = joinpath("data_folder", linear ? "L_hard.json" : "NL_hard.json")
        open(file_path, "w") do file
            println(file, json_data)
        end
    elseif method == "CVaR"
        #--------------------- CVaR Method --------------------------#
        model = runmodel(method, linear, α)
        Xs = value.(model[:x])
        ϕ_opt = vec(value(model[:ϕ]))
        T_opt = vec(value(model[:T]))
        u_opt = vec(value(model[:u]))
        h_opt = vec(value(model[:h]))
        solve_time = JuMP.solve_time(model)
        optimizer = solver_name(model)
        println(optimizer)

        data = Dict(    
            # Save CVaR Results:
            "T_opt" => T_opt, "h_opt" => h_opt, 
            "ϕ_opt" => ϕ_opt, "u_opt" => u_opt, 
            "λ values" => λ_vals, "xs" => Xs,
            "solver time" => round(solve_time, digits = 5),
            "objective value" => objective_value(model),
            "solver name" => linear ? "Gurobi" : "GAMS - dicopt",
            "termination status" => termination_status(model),
        )
        json_data = JSON.json(data)
        file_path = joinpath("data_folder", linear ? "L_CVaR_alpha_$α.json" : "NL_CVaR_alpha_$α.json")
        open(file_path, "w") do file
            println(file, json_data)
        end
    elseif method == "SigVaR"

        #--------------------- SigVaR Method --------------------------#
        model = runmodel(method, linear, α)
        Xs = value.(model[:x])
        ϕ_opt = vec(value(model[:ϕ]))
        T_opt = vec(value(model[:T]))
        u_opt = vec(value(model[:u]))
        h_opt = vec(value(model[:h]))
        optimizer = solver_name(model)
        println(optimizer)

        data = Dict(            
            # Save SigVaR Results:
            "τ values" => τ_vals, "μ values" => μ_vals, 
            "ϕ values" => ϕ_vals, "h values" => h_vals,
            "T values" => T_vals, "u values" => u_vals,
            "iterations" => iters, "termination status" => status_vals, 
    
            "T_opt" => T_opt, "h_opt" => h_opt, 
            "ϕ_opt" => ϕ_opt, "u_opt" => u_opt, 
            "xs" => Xs, 
            "solver times" => round.(runtime_vals, digits = 5),
            "objective values" => obj_vals,
            "objective value" => obj_vals[end],
            "solver name" => "GAMS - conopt4",
        )
        json_data = JSON.json(data)
        file_path = joinpath("data_folder", linear ? "L_SigVaR_alpha_$α.json" : "NL_SigVaR_alpha_$α.json")
        open(file_path, "w") do file
            println(file, json_data)
        end
    elseif method == "bigM"
        model = runmodel(method, linear, α)
        Xs = value.(model[:x])
        T_opt = vec(value(model[:T]))
        u_opt = vec(value(model[:u]))
        h_opt = vec(value(model[:h]))
        q_opt = vec(value(model[:q]))
        solve_time = JuMP.solve_time(model)
        optimizer = solver_name(model)
        println(optimizer)

        data = Dict(
            "xs" => Xs, "T_opt" => T_opt, "h_opt" => h_opt, 
            "q_opt" => q_opt, "u_opt" => u_opt, 
            "solver time" => round(solve_time, digits = 5),
            "objective value" => objective_value(model),
            "solver name" => "Gurobi",
        )
        json_data = JSON.json(data)
        file_path = joinpath("data_folder", linear ? "L_bigM_gurobi_alpha_$α.json" : "NL_bigM_gurobi_alpha_$α.json")
        open(file_path, "w") do file
            println(file, json_data)
        end
    elseif method == "bigM_juniper"
        #--------------------- BigM Juniper Method --------------------------#
        model = runmodel(method, linear, α)
        Xs = value.(model[:x])
        T_opt = vec(value(model[:T]))
        u_opt = vec(value(model[:u]))
        h_opt = vec(value(model[:h]))
        q_opt = vec(value(model[:q]))
        solve_time = JuMP.solve_time(model)
        optimizer = solver_name(model)
        println(optimizer)

        data = Dict(
            # Save BigM Juniper Results:
            "xs" => Xs, "T_opt" => T_opt, "h_opt" => h_opt, 
            "q_opt" => q_opt, "u_opt" => u_opt, 
            "solver time" => round(solve_time, digits = 5),
            "objective value" => objective_value(model),
            "solver name" => "Juniper-Ipopt",
        )
        json_data = JSON.json(data)
        file_path = joinpath("data_folder", linear ? "L_bigM_juniper_alpha_$α.json" : "NL_bigM_juniper_alpha_$α.json")
        open(file_path, "w") do file
            println(file, json_data)
        end
    elseif method == "bigM_dicopt"
        #--------------------- BigM Dicopt Method --------------------------#
        model = runmodel(method, linear, α)
        Xs = value.(model[:x])
        T_opt = vec(value(model[:T]))
        u_opt = vec(value(model[:u]))
        h_opt = vec(value(model[:h]))
        q_opt = vec(value(model[:q]))
        solve_time = JuMP.solve_time(model)
        optimizer = solver_name(model)
        println(optimizer)

        data = Dict(
            # Save BigM Dicopt Results:
            "xs" => Xs, "T_opt" => T_opt, "h_opt" => h_opt, 
            "q_opt" => q_opt, "u_opt" => u_opt, 
            "solver time" => round(solve_time, digits = 5),
            "objective value" => objective_value(model),
            "solver name" => "GAMS - dicopt",
        )
        json_data = JSON.json(data)
        file_path = joinpath("data_folder", linear ? "L_bigM_dicopt_alpha_$α.json" : "NL_bigM_dicopt_alpha_$α.json")
        open(file_path, "w") do file
            println(file, json_data)
        end
    elseif method == "bigM_scip"
        #--------------------- BigM SCIP Method --------------------------#
        model = runmodel(method, linear, α)
        Xs = value.(model[:x])
        T_opt = vec(value(model[:T]))
        u_opt = vec(value(model[:u]))
        h_opt = vec(value(model[:h]))
        q_opt = vec(value(model[:q]))
        solve_time = JuMP.solve_time(model)
        optimizer = solver_name(model)
        println(optimizer)

        data = Dict(
            # Save BigM SCIP Results:
            "xs" => Xs, "T_opt" => T_opt, "h_opt" => h_opt, 
            "q_opt" => q_opt, "u_opt" => u_opt, 
            "solver time" => round(solve_time, digits = 5),
            "objective value" => objective_value(model),
            "solver name" => "SCIP",
        )
        json_data = JSON.json(data)
        file_path = joinpath("data_folder", linear ? "L_bigM_scip_alpha_$α.json" : "NL_bigM_scip_alpha_$α.json")
        open(file_path, "w") do file
            println(file, json_data)
        end
    elseif method == "bigM_baron"
        #--------------------- BigM Baron Method --------------------------#
        model = runmodel(method, linear, α)
        Xs = value.(model[:x])
        T_opt = vec(value(model[:T]))
        u_opt = vec(value(model[:u]))
        h_opt = vec(value(model[:h]))
        q_opt = vec(value(model[:q]))
        solve_time = JuMP.solve_time(model)
        optimizer = solver_name(model)
        println(optimizer)

        data = Dict(
            # Save BigM Baron Results:
            "xs" => Xs, "T_opt" => T_opt, "h_opt" => h_opt, 
            "q_opt" => q_opt, "u_opt" => u_opt, 
            "solver time" => round(solve_time, digits = 5),
            "objective value" => objective_value(model),
            "solver name" => "baron",
        )
        json_data = JSON.json(data)
        file_path = joinpath("data_folder", linear ? "L_bigM_baron_alpha_$α.json" : "NL_bigM_baron_alpha_$α.json")
        open(file_path, "w") do file
            println(file, json_data)
        end
        #--------------------- BigM Method --------------------------#
        model = runmodel(method, linear, α)
        Xs = value.(model[:x])
        T_opt = vec(value(model[:T]))
        u_opt = vec(value(model[:u]))
        h_opt = vec(value(model[:h]))
        q_opt = vec(value(model[:q]))
        solve_time = JuMP.solve_time(model)
        optimizer = solver_name(model)
        println(optimizer)

        data = Dict(
            # Save BigM Results:
            "xs" => Xs, "T_opt" => T_opt, "h_opt" => h_opt, 
            "q_opt" => q_opt, "u_opt" => u_opt, 
            "solver time" => round(solve_time, digits = 5),
            "objective value" => objective_value(model),
            "solver name" => linear ? "Gurobi" : "GAMS - dicopt",
        )
        json_data = JSON.json(data)
        file_path = joinpath("data_folder", linear ? "L_bigM_alpha_$α.json" : "NL_bigM_alpha_$α.json")
        open(file_path, "w") do file
            println(file, json_data)
        end
    elseif method == "complimentarity_ipopt"
        # -------------- Complimentarity Method -------------- #
        model = runmodel(method, linear, α)
        Xs = value.(model[:x])
        y1_opt = vec(value(model[:y1]))
        y2_opt = vec(value(model[:y2]))
        T_opt = vec(value(model[:T]))
        u_opt = vec(value(model[:u]))
        h_opt = vec(value(model[:h]))
        optimizer = solver_name(model)
        println(optimizer)

        data = Dict(
            # Save Complimentarity Results:
            "xs" => Xs, "T_opt" => T_opt, "h_opt" => h_opt, 
            "y1_opt" => y1_opt, "y2_opt" => y2_opt, 
            "u_opt" => u_opt, 
            "solver name" => "ipopt",
            "termination status" => comp_status_vals, 
            "solver times" => round.(comp_runtime_vals, digits = 5),
            "objective values" => comp_obj_vals,
            "objective value" => comp_obj_vals[end],
            "ϵ values" => ϵ_vals[end],
        )
        json_data = JSON.json(data)
        file_path = joinpath("data_folder", linear ? "L_compliment_alpha_$(α)_ipopt.json" : "NL_compliment_alpha_$(α)_ipopt.json")
        open(file_path, "w") do file
            println(file, json_data)
        end
    elseif method == "complimentarity_conopt4"
        # -------------- Complimentarity Method -------------- #
        model = runmodel(method, linear, α)
        Xs = value.(model[:x])
        y1_opt = vec(value(model[:y1]))
        y2_opt = vec(value(model[:y2]))
        T_opt = vec(value(model[:T]))
        u_opt = vec(value(model[:u]))
        h_opt = vec(value(model[:h]))
        optimizer = solver_name(model)
        println(optimizer)

        data = Dict(
            # Save Complimentarity Results:
            "xs" => Xs, "T_opt" => T_opt, "h_opt" => h_opt, 
            "y1_opt" => y1_opt, "y2_opt" => y2_opt, 
            "u_opt" => u_opt, 
            "solver name" => "conopt4",
            "termination status" => comp_status_vals, 
            "solver times" => round.(comp_runtime_vals, digits = 5),
            "objective values" => comp_obj_vals,
            "objective value" => comp_obj_vals[end],
            "ϵ values" => ϵ_vals[end],
        )
        json_data = JSON.json(data)
        file_path = joinpath("data_folder", linear ? "L_compliment_alpha_$(α)_conopt4.json" : "NL_compliment_alpha_$(α)_conopt4.json")
        open(file_path, "w") do file
            println(file, json_data)
        end
    end
    clear_global_vectors()
end
αs = [0.9, 0.95, 0.96, 0.97, 0.99, 0.999]

for α in αs
    #------------- Linear Solution --------------#
    linear = true
    master(α, "Hard", linear) # Hard Constraint Method
    master(α, "CVaR", linear) # CVaR Method
    master(α, "SigVaR", linear) # SigVaR Method
    master(α, "bigM", linear) # Big-M Method
    master(α, "complimentarity_ipopt", linear) # Complimentarity Method Ipopt
    master(α, "complimentarity_conopt4", linear) # Complimentarity Method Conopt4

    #------------- Nonlinear Solution --------------#
    linear = false
    master(α, "Hard", linear) # Hard Constraint Method
    master(α, "CVaR", linear) # CVaR Method
    master(α, "SigVaR", linear) # SigVaR Method
    # master(α, "bigM_juniper", linear) # Big-M Method Juniper
    # master(α, "bigM_dicopt", linear) # Big-M Method Dicopt
    # master(α, "bigM_scip", linear) # Big-M Method SCIP
    # master(α, "bigM_baron", linear) # Big-M Method Baron
    master(α, "complimentarity_ipopt", linear) # Complimentarity Method Ipopt
    master(α, "complimentarity_conopt4", linear) # Complimentarity Method Conopt4
end

# Only run the big-M methods for the nonlinear case:
for α in αs
    linear = false
    master(α, "bigM_juniper", linear) # Big-M Method Juniper
    master(α, "bigM_dicopt", linear) # Big-M Method Dicopt
    master(α, "bigM_scip", linear) # Big-M Method SCIP
    master(α, "bigM_baron", linear) # Big-M Method Baron
end