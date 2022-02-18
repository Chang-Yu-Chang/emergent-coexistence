import re
import pandas as pd
import numpy as np
import os
os.chdir('simulation/')
# from community_simulator import *
# from community_simulator.usertools import *
from user_tools import *

input_csv = 'input_independent.csv' # Input file name
row_number = 1 # Which row of experiment to run
input_row = pd.read_csv(input_csv).loc[row_number]
seed = int(input_row['seed'])
exp_id = int(input_row['exp_id'])
input_row['ma'] = 2
assumptions = load_assumptions(input_row)
np.random.seed(seed)
# D, c, l = sample_matrices(a)
# write_matrices(input_row, D, c, l)
M, T, S, F, type_names, resource_index, consumer_index = extract_shapes(assumptions)
M_waste = assumptions['MA'][assumptions['waste_type']]
waste_name = type_names[assumptions['waste_type']]

# 
# # Set family specific secretion
# assumptions['fs'] = 0.1
# assumptions['fw'] = 0.7
assumptions['fss'] = 0.1 # sugar to sugar
assumptions['fsa'] = 0.8 # sugar to acid
assumptions['fsw'] = 0.1 # sugar to waste
assumptions['fas'] = 0.01   # acid to sugar
assumptions['faa'] = 0.10 # acid to acid
assumptions['faw'] = 0.79 # acid to waste
assumptions['fws'] = 0.01   # waste to sugar
assumptions['fwa'] = 0.01 # waste to acid
assumptions['fww'] = 0.98 # waste to waste
# 

DT = pd.DataFrame(np.zeros((M,M)),index=resource_index,columns=resource_index)
for type_name in type_names:
    MA = len(DT.loc[type_name])
    if type_name == 'T0': #sugar
        p = pd.Series(np.zeros(M), index = DT.keys())
        p.loc['T0'] = assumptions['fss']/assumptions['MA'][0] # to sugar
        p.loc['T1'] = assumptions['fsa']/assumptions['MA'][1] # to acid
        p.loc['T2'] = assumptions['fsw']/assumptions['MA'][2] # to waste
        DT.loc[type_name] = dirichlet(p/assumptions['sparsity'],size=assumptions['MA'][0])
    elif type_name == 'T1': # acid
        p = pd.Series(np.zeros(M), index = DT.keys())
        p.loc['T0'] = assumptions['fas']/assumptions['MA'][0] # to sugar
        p.loc['T1'] = assumptions['faa']/assumptions['MA'][1] # to acid
        p.loc['T2'] = assumptions['faw']/assumptions['MA'][2] # to waste
        DT.loc[type_name] = dirichlet(p/assumptions['sparsity'],size=assumptions['MA'][1])
    elif type_name == 'T2': # waste
        p = pd.Series(np.zeros(M), index = DT.keys())
        p.loc['T0'] = assumptions['fws']/assumptions['MA'][0] # to sugar
        p.loc['T1'] = assumptions['fwa']/assumptions['MA'][1] # to acid
        p.loc['T2'] = assumptions['fww']/assumptions['MA'][2] # to waste
        DT.loc[type_name] = dirichlet(p/assumptions['sparsity'],size=assumptions['MA'][2])
    else:
        print("D matrices are hard coded. The number of resource class must be exactly 3.")
D = DT.T

D.round(2)

# Sample matrices
np.random.seed(seed)
D, c, l = sample_matrices(a)
write_matrices(input_row, D, c, l)

# Update params 
params = MakeParams(a)
params['c'] = c
params['D'] = D
params['l'] = l

# Generate equations
def dNdt(N,R,params):
    return make_consumer_dynamics(a)(N,R,params) # Use the updated version to include varied l
def dRdt(N,R,params):
    return make_resource_dynamics(a)(N,R,params) # Use the updated version to include varied D and l
dynamics = [dNdt,dRdt]

# Set initial state by reading in plate conditions
N0, R0 = make_initial_state(input_row, a)
init_state = [N0,R0]

# Make plate object
Plate = Community(init_state, dynamics, params, parallel = False)
Plate.RunExperiment(np.eye(a['n_wells'])/10, T = 24, npass = 5, refresh_resource=True)
Plate.N.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "end.csv", input_row["init_N0"]))










