from idaes import *
from pyomo.environ import *
from pyomo.gdp import *
from time import time
import numpy as np
import json

def ieee14(alphas, num_samples, method, active_gens=4, active_lines=20, time_limit=-1, tee=False):
    # Set the covariance matrix for the uncertain parameters
    theta_nom = np.array([87.3, 50.0, 25.0, 28.8, 50.0, 25.0, 0, 0, 0, 0, 0])
    covar = np.identity(len(theta_nom)) * 1200.
    covar[covar == 0] = 240.0

    # Specify the network details
    line_cap = 50.
    gen_cap = [332., 140., 100., 100., 100.]

    # Set the problem parameters
    num_lines = 20
    num_gens = 5

    # Obtain samples
    np.random.seed(42)
    thetas = np.random.multivariate_normal(theta_nom, covar, num_samples)
    thetas[thetas <= 0.0] = 0.0

    # np.savetxt('thetas.txt', thetas)

    start = time()

    # Define the model
    m = ConcreteModel()

    # Create the indexing sets
    m.K = RangeSet(0, num_samples-1)
    m.G = RangeSet(1, num_gens)
    m.L = RangeSet(1, num_lines)

    # Create the variables
    m.z_line = Var(m.L, m.K, bounds=(-150, 150)) # These are extra bounds needed for the reformulations
    m.z_gen = Var(m.G, m.K, bounds=(0, 632)) # These are extra bounds needed for the reformulations
    m.d_line = Var(m.L, bounds=(0, 100))
    m.d_gen = Var(m.G, bounds=(0, 300))

    # Create the objectivve function
    m.obj = Objective(sense=minimize, expr = sum(m.d_gen[g] for g in m.G) + sum(m.d_line[l] for l in m.L))

    # Create the balance constraints
    @m.Constraint(m.K)
    def h1(m, k):
        return m.z_gen[1, k] - m.z_line[1, k] - m.z_line[6, k] == 0
    @m.Constraint(m.K)
    def h2(m, k):
        return m.z_gen[2, k] + m.z_line[1, k] - sum(m.z_line[i, k] for i in [2, 4, 5]) - thetas[k, 0] == 0
    @m.Constraint(m.K)
    def h3(m, k):
        return m.z_gen[3, k] + m.z_line[2, k] - m.z_line[3, k] - thetas[k, 1] == 0
    @m.Constraint(m.K)
    def h4(m, k):
        return sum(m.z_line[i, k] for i in [3, 4, 8]) - sum(m.z_line[i, k] for i in [7, 11]) - thetas[k, 2] == 0
    @m.Constraint(m.K)
    def h5(m, k):
        return sum(m.z_line[i, k] for i in [5, 6, 7, 12]) - thetas[k, 3] == 0
    @m.Constraint(m.K)
    def h6(m, k):
        return m.z_gen[4, k] + sum(m.z_line[i, k] for i in [16, 18]) - sum(m.z_line[i, k] for i in [12, 19]) - thetas[k, 4] == 0
    @m.Constraint(m.K)
    def h7(m, k):
        return m.z_line[9, k] - sum(m.z_line[i, k] for i in [8, 10]) == 0
    @m.Constraint(m.K)
    def h8(m, k):
        return m.z_gen[5, k] - m.z_line[9, k] == 0
    @m.Constraint(m.K)
    def h9(m, k):
        return sum(m.z_line[i, k] for i in [10, 11]) - sum(m.z_line[i, k] for i in [13, 14]) - thetas[k, 5] == 0
    @m.Constraint(m.K)
    def h10(m, k):
        return sum(m.z_line[i, k] for i in [13, 20]) - thetas[k, 6] == 0
    @m.Constraint(m.K) 
    def h11(m, k):
        return m.z_line[19, k] - m.z_line[20, k] - thetas[k, 7] == 0
    @m.Constraint(m.K)
    def h12(m, k):
        return m.z_line[17, k] - m.z_line[18, k] - thetas[k, 8] == 0
    @m.Constraint(m.K)
    def h13(m, k):
        return m.z_line[15, k] - sum(m.z_line[i, k] for i in [16, 17]) - thetas[k, 9] == 0
    @m.Constraint(m.K)
    def h14(m, k):
        return m.z_line[14, k] - m.z_line[15, k] - thetas[k, 10] == 0

    # Build the generator disjunctions
    def build_satisfy_gen_constraints(disjunct, g, k):
        m = disjunct.model()
        @disjunct.Constraint()
        def gen_upper(disjunct):
            return m.z_gen[g, k] - gen_cap[g-1] - m.d_gen[g] <= 0

    def build_not_satisfy_gen_constraints(disjunct, g, k):
        m = disjunct.model()
        @disjunct.Constraint()
        def gen_upper(disjunct):
            return m.z_gen[g, k] - gen_cap[g-1] - m.d_gen[g] >= 0

    m.gen_disjunct1 = Disjunct(m.G, m.K, rule=build_satisfy_gen_constraints)
    m.gen_disjunct2 = Disjunct(m.G, m.K, rule=build_not_satisfy_gen_constraints)

    @m.Disjunction(m.G, m.K)
    def gen_constrs_satisfy_or_not(m, g, k):
        return [m.gen_disjunct1[g, k], m.gen_disjunct2[g, k]]

    # Build the line disjunctions
    def build_satisfy_line_constraints(disjunct, l, k):
        m = disjunct.model()

        @disjunct.Constraint()
        def line_lower(disjunct):
            return -m.z_line[l, k] - line_cap - m.d_line[l] <= 0
        
        @disjunct.Constraint()
        def line_upper(disjunct):
            return m.z_line[l, k] - line_cap - m.d_line[l] <= 0

    def build_not_satisfy_lower_line_constraints(disjunct, l, k):
        m = disjunct.model()
        @disjunct.Constraint()
        def line_lower(disjunct):
            return -m.z_line[l, k] - line_cap - m.d_line[l] >= 0

    def build_not_satisfy_upper_line_constraints(disjunct, l, k):
        m = disjunct.model()
        @disjunct.Constraint()
        def line_upper(disjunct):
            return m.z_line[l, k] - line_cap - m.d_line[l] >= 0

    m.line_disjunct1 = Disjunct(m.L, m.K, rule=build_satisfy_line_constraints)
    m.line_disjunct2 = Disjunct(m.L, m.K, rule=build_not_satisfy_lower_line_constraints)
    m.line_disjunct3 = Disjunct(m.L, m.K, rule=build_not_satisfy_upper_line_constraints)

    @m.Disjunction(m.L, m.K)
    def line_constrs_satisfy_or_not(m, l, k):
        return [m.line_disjunct1[l, k], m.line_disjunct2[l, k], m.line_disjunct3[l, k]]

    # Set the event constraint (this is were we play with the logic) --> if we change then update the JSON file name below
    m.w = BooleanVar(m.K)
    # @m.LogicalConstraint(m.K)
    # def event_logic(m, k):
    #     return equivalent(land(land(*m.gen_disjunct1[:, k].indicator_var), land(*m.line_disjunct1[:, k].indicator_var)), m.w[k])

    @m.LogicalConstraint(m.K)
    def event_logic(m, k):
        return equivalent(land(atleast(active_gens,*m.gen_disjunct1[:, k].indicator_var), atleast(active_lines, *m.line_disjunct1[:, k].indicator_var)), m.w[k])


    # Enforce the event constraint over scenarios
    m.min_constrs = Param(mutable=True, default=0)
    m.event_constr = LogicalConstraint(expr = atleast(m.min_constrs, m.w))

    # Transform the formulation
    if method in ('gdp.bigm', 'gdp.hull', 'gdp.cuttingplane'):
        TransformationFactory(method).apply_to(m)
    build_time = time() - start

    # Create Dict to store stuff in 
    data = {}

    # Solve the model for each alpha
    for alpha in alphas:
        # Solve the model and get results
        m.min_constrs.value = ceil(alpha * num_samples)
        if method in ('LOA', 'LBB', 'GLOA', 'RIC'):
            solver = SolverFactory('gdpopt')
            if time_limit != -1:
                results = solver.solve(m, tee=tee, strategy=method, mip_solver='gurobi_direct', time_limit=time_limit)
            else:
                results = solver.solve(m, tee=tee, strategy=method, mip_solver='gurobi_direct')
        elif method in ('gdp.bigm', 'gdp.hull', 'gdp.cuttingplane'):
            solver = SolverFactory('gurobi_direct')
            if time_limit != -1:
                solver.options['TimeLimit'] = time_limit
            results = solver.solve(m, tee=tee, report_timing=False)
        else:
            raise ValueError('Unrecognized method ' + method)
        opt_time = results.solver.wallclock_time
        status = results.solver.termination_condition

        # Get the results (hit time limit)
        if str(status) != 'optimal':
            data[alpha] = {'time' : opt_time, 'status' : status}

            print("\n--------------------------------------------------")
            print("---------------------RESULTS----------------------")
            print("--------------------------------------------------")
            print('Build & Transform Time: ', build_time, 's')
            print('Solution Method:        ', method)
            print('Status:                 ', status)
            print("~Solution Time:         ", opt_time, 's')
            print("Probability Level:      ", alpha)
            print("--------------------------------------------------\n")

            continue

        # Get the results for normal solve
        opt_obj = value(m.obj)
        data[alpha] = {'time' : opt_time, 'status' : status, 'objective' : opt_obj, 
                    'd_line' : [value(m.d_line[idx]) for idx in m.d_line],
                    'd_gen' : [value(m.d_gen[idx]) for idx in m.d_gen]}

        # Print the results
        print("\n--------------------------------------------------")
        print("---------------------RESULTS----------------------")
        print("--------------------------------------------------")
        print('Build & Transform Time: ', build_time, 's')
        print('Solution Method:        ', method)
        print('Status:                 ', status)
        print("~Solution Time:         ", opt_time, 's')
        print("Probability Level:      ", alpha)
        print("Optimal Objective:      ", opt_obj)
        print("Optimal Line Design Values: ")
        for idx in m.d_line:
            print('\tLine ', idx, ': ', value(m.d_line[idx]))
        print("Optimal Generator Design Values: ")
        for idx in m.d_gen:
            print('\tGenerator ', idx, ': ', value(m.d_gen[idx]))
        print("--------------------------------------------------\n")

    # Create the JSON file (change the name to for different logic)
    if method[:3] == 'gdp':
        method_name = method[4:]
    else:
        method_name = method

    event_name = str(active_gens)+'gens_'+str(active_lines)+'lines'
    with open('./data/ieee14_atleast_' + event_name + '_' + method_name + '_' + str(num_samples) + '.json', 'w') as fp:
        json.dump(data, fp)

    # Return the data
    return data

# Setup data for runs
alphas = [0.5, 0.6, 0.7, 0.8, 0.85, 0.87, 0.89, 0.9, 0.92, 0.94, 0.95, 0.97, 0.98, 0.99, 0.9999]
num_samples = 1000
methods = ['gdp.bigm']
active_gens_list = [5]
active_lines_list = [20,19,18,17,16,15]

# Run the models
for method in methods:
    for active_gens in active_gens_list:
        for active_lines in active_lines_list:
            ieee14(alphas, num_samples, method, active_gens=active_gens, active_lines=active_lines, time_limit=3600)
