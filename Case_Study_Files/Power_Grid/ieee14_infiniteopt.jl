using InfiniteOpt, LinearAlgebra, Gurobi, Distributions, LinearAlgebra, Random, DisjunctiveProgramming, JSON

"""
    ieee14_gdp

Run the ieee14 case study using GDP over a vector of alpha values `αs` with a '
`num_samples` as the number random scenarios. This returns a dictionary with 
all the results.
"""
function ieee14_gdp(
    αs, 
    num_samples, 
    method = BigM(); # the method used to reformulate the model --> see DisjunctiveProgramming.jl 
    atleast = true, # should atleast logic be used instead of intersection logic?
    min_gen = 5, # if atleast logic is used, how many generators are needed?
    min_line = 20, # if atleast logic is used, how many lines are needed?
    time_limit = 3600, # what solver time limit should be used (in seconds)
    silent_solver = true, # should we silence the raw solver output?
    mip_gap = 0.00,
    mip_gap_abs = 0.0,
    )

    # Set the covariance matrix for the uncertain parameters
    ξ_nom = [87.3; 50.0; 25.0; 28.8; 50.0; 25.0; 0; 0; 0; 0; 0]
    β = 240.0
    covar = Matrix(I, length(ξ_nom), length(ξ_nom)) * 1200.
    covar[covar .== 0] .= β

    # Specify the network details
    line_cap = 50
    gen_cap = [332; 140; 100; 100; 100]

    # Set the problem parameters
    M = 10000
    num_lines = 20
    num_gens = 5
    sig_dig = 4

    # Obtain samples
    Random.seed!(42)
    d = MvNormal(ξ_nom, covar)
    ξs = rand(d, num_samples)
    ξs[ξs .<= 0] .= 0.0
    ξs = round.(ξs, digits=5)

    # Set up the Model
    m = InfiniteGDPModel(Gurobi.Optimizer)
    if silent_solver
        set_silent(m)
    end
    set_time_limit_sec(m, time_limit)
    set_optimizer_attribute(m, "MIPGap", mip_gap)
    set_optimizer_attribute(m, "MIPGapAbs", mip_gap_abs)

    # Set parameters
    @infinite_parameter(m, ξ[1:length(ξ_nom)] ~ d)
    set_supports(ξ, ξs, label = WeightedSample) # lets it know the samples are from the distribution
    @finite_parameter(m, α == 0) # makes it so we can change it later and resolve if wanted

    # Set up the variables
    @variable(m, -1000 <= z_line[1:num_lines] <= 1000, Infinite(ξ))
    @variable(m, 0 <= z_gen[1:num_gens] <= 1000, Infinite(ξ))
    @variable(m, 0 <= d_line[1:num_lines] <= 100)
    @variable(m, 0 <= d_gen[1:num_gens] <= 300)

    # Set the infinite logical variables
    @variable(m, Y_gen[1:num_gens], InfiniteLogical(ξ)) # generator constaints hold
    @variable(m, notY_gen[i=1:num_gens], InfiniteLogical(ξ), logical_complement=Y_gen[i])
    @variable(m, Y_line_low[1:num_lines], InfiniteLogical(ξ))
    @variable(m, Y_line[1:num_lines], InfiniteLogical(ξ)) # line constraints hold
    @variable(m, Y_line_high[1:num_lines], InfiniteLogical(ξ))
    @variable(m, W, InfiniteLogical(ξ))

    # Set objective function
    @objective(m, Min, sum(d_gen) + sum(d_line) ) # + 1

    # Set the node balance constraints
    @constraints(m, begin 
        z_gen[1] - z_line[1] - z_line[6] == 0
        z_gen[2] + z_line[1] - sum(z_line[i] for i = [2; 4; 5]) - ξ[1] == 0
        z_gen[3] + z_line[2] - z_line[3] - ξ[2] == 0
        sum(z_line[i] for i = [3; 4; 8]) - sum(z_line[i] for i = [7; 11]) - ξ[3] == 0
        sum(z_line[i] for i = [5; 6; 7; 12]) - ξ[4] == 0
        z_gen[4] + sum(z_line[i] for i = [16; 18]) - sum(z_line[i] for i = [12; 19]) - ξ[5] == 0
        z_line[9] - sum(z_line[i] for i = [8; 10]) == 0
        z_gen[5] - z_line[9] == 0
        sum(z_line[i] for i = [10; 11]) - sum(z_line[i] for i = [13; 14]) - ξ[6] == 0
        sum(z_line[i] for i = [13; 20]) - ξ[7] == 0
        z_line[19] - z_line[20] - ξ[8] == 0
        z_line[17] - z_line[18] - ξ[9] == 0
        z_line[15] - sum(z_line[i] for i = [16; 17]) - ξ[10] == 0
        z_line[14] - z_line[15] - ξ[11] == 0
    end)

    # Disjunctions for generator capacity constraints
    @constraint(m, [g ∈ 1:num_gens], z_gen[g] - gen_cap[g] - d_gen[g] <= 0, Disjunct(Y_gen[g])) 
    @constraint(m, [g ∈ 1:num_gens], z_gen[g] - gen_cap[g] - d_gen[g] >= 0, Disjunct(notY_gen[g]))
    @disjunction(m, [g ∈ 1:num_gens], [Y_gen[g], notY_gen[g]])
    # @disjunction(m, [g ∈ 1:num_gens], [Y_gen[g]], exactly1=false)

    # Disjunctions for line capacity constraints
    @constraint(m, [l ∈ 1:num_lines], -z_line[l] - line_cap - d_line[l] >= 0, Disjunct(Y_line_low[l]))
    @constraint(m, [l ∈ 1:num_lines], -line_cap - d_line[l] <= z_line[l], Disjunct(Y_line[l]))
    @constraint(m, [l ∈ 1:num_lines], z_line[l] <= line_cap + d_line[l], Disjunct(Y_line[l]))
    @constraint(m, [l ∈ 1:num_lines], z_line[l] - line_cap - d_line[l] >= 0, Disjunct(Y_line_high[l]))
    @disjunction(m, [l ∈ 1:num_lines], [Y_line_low[l], Y_line[l], Y_line_high[l]])
    # @disjunction(m, [l ∈ 1:num_lines], [Y_line[l]] , exactly1=false)

    # Set the event constraint aggregation logic
    if !atleast
        @constraint(m, (logical_and(Y_gen...) ∧ logical_and(Y_line...)) ⇔ W := true)
    else
        @variable(m, Y_atleast_gen, InfiniteLogical(ξ))
        @variable(m, Y_atleast_line, InfiniteLogical(ξ))

        @variable(m, notY_atleast_gen, InfiniteLogical(ξ), logical_complement=Y_atleast_gen)
        @variable(m, dummy_gens >= 0, Infinite(ξ))
        @constraint(m, dummy_gens == sum(binary_variable.(Y_gen)))
        @constraint(m, dummy_gens >= min_gen, Disjunct(Y_atleast_gen)) 
        @constraint(m, dummy_gens <= min_gen - 1, Disjunct(notY_atleast_gen))
        @disjunction(m, [Y_atleast_gen, notY_atleast_gen])
        # @disjunction(m, [Y_atleast_gen], exactly1=false)

        @variable(m, notY_atleast_line, InfiniteLogical(ξ), logical_complement=Y_atleast_line)
        @variable(m, dummy_lines >= 0, Infinite(ξ))
        @constraint(m, dummy_lines == sum(binary_variable.(Y_line)))
        @constraint(m, dummy_lines >= min_line, Disjunct(Y_atleast_line)) 
        @constraint(m, dummy_lines <= min_line - 1, Disjunct(notY_atleast_line))
        @disjunction(m, [Y_atleast_line, notY_atleast_line])
        # @disjunction(m, [Y_atleast_line], exactly1=false)

        # Manueal bigM
        # @constraint(m, min_gen - sum(binary_variable.(Y_gen)) ≤ (1 - binary_variable(Y_atleast_gen)) * min_gen)
        # @constraint(m, min_line - sum(binary_variable.(Y_line)) ≤ (1 - binary_variable(Y_atleast_line)) * min_line)

        @constraint(m, (Y_atleast_gen ∧ Y_atleast_line) ⇔ W := true)
    end

    # Enforce the event constraint
    @expression(m, prob, 𝔼(binary_variable(W), ξ)) # makes it easy to query the value later on
    @constraint(m, event, prob >= α)

    # Prepare the data dict
    data = Dict()

    # Solve and store the results
    for α_val in αs
        set_value(α, α_val)
        optimize!(m)

        if has_values(m)
            opt_d_line = value.(d_line)
            opt_d_gen = value.(d_gen)
            opt_obj = objective_value(m)
            opt_time = solve_time(m)
            opt_prob = value(prob)

            # Print the results
            println("------------------RESULTS------------------")
            println("Solution Method:       ", method)
            println("Optimal Objective:     ", round(opt_obj, sigdigits = sig_dig))
            println("Solution Time:         ", round(opt_time, sigdigits = sig_dig), "s")
            println("Minimum Probability:   ", value(α))
            println("Predicted Probability: ", 100 * round(opt_prob, sigdigits = sig_dig), "%")
            println("Optimal Line Design Values: ", opt_d_line)
            println("Optimal Generator Design Values: ", opt_d_gen, "\n")
            # println("My print", value.(Y_atleast_gen))

            # save the data
            data[α_val] = Dict("time" => opt_time, "status" => termination_status(m), 
                               "objective" => opt_obj, "d_line" => opt_d_line,
                               "d_gen" => opt_d_gen, "w" => value(W))
        else
            println("Failed to converge with α = $α")
            data[α_val] = Dict("time" => opt_time, "status" => termination_status(m))
        end
    end

    # save data
    logic = atleast ? "atleast_" : "intersection_"
    event = atleast ? "$(min_gen)gens_$(min_line)lines_" : ""
    method = string(method)[1:findfirst(c -> c == '{' || c == '(', string(method))-1] * "_"
    open("./data/ieee14_disjlogcompl_" * logic * event * method * string(num_samples, ".json"), "w") do f
        write(f, JSON.json(data))
    end
    return data
end

# Set the run parameters
alphas = [0.5, 0.6, 0.7, 0.8, 0.85, 0.87, 0.89, 0.9, 0.92, 0.94, 0.95, 0.97, 0.98, 0.99, 0.9999]
num_samples = 1000
methods = [BigM(), Indicator(), Hull()]
active_gens_list = [5,4,3]
active_lines_list = [20,19,18,17,16,15]
atleast_logic_flag = [true, false]

# Run the models
for method in methods
    for atleast_logic in atleast_logic_flag
        if atleast_logic
            for active_gens in active_gens_list
                for active_lines in active_lines_list
                    ieee14_gdp(alphas, num_samples, method, min_gen = active_gens, min_line = active_lines, time_limit = 10000, silent_solver=true)
                end
            end
        else
            ieee14_gdp(alphas, num_samples, method, atleast = atleast_logic, time_limit = 10000, silent_solver=false) # uses intersection logic
        end
    end
end