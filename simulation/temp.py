import sys
import pandas as pd
import numpy as np
from community_simulator import *
from community_simulator.usertools import *


input_row = pd.read_csv("input_communityPairs.csv").loc[0]

### Load Parameters from the input csv
output_dir = input_row['output_dir']
q = input_row['q']
q2 = input_row['q2']
seed = int(input_row['seed'])
exp_id = int(input_row['exp_id'])
S = int(input_row['S'])
np.random.seed(seed)

# 
N0_init = pd.read_csv(output_dir + 'communityPairs_' + str(exp_id) + '_init.csv')
n_pairs = int(N0_init.shape[1])

#Set paramaters in model note some of these paramaters are irrelevant to simulation (c1,c0 etc).

a = {'sampling':'Gamma', #Sampling method
	'SA': np.ones(3)*100, #Number of species in each family
	'MA': np.ones(3)*10, #Number of resources of each type
	'Sgen': 0, #Number of generalist species
	'muc': 10, #Mean sum of consumption rates in Gaussian model
	'q': q, #Preference strength (0 for generalist and 1 for specialist) on sugars
	'q2': q2, #Preference strength (0 for generalist and 1 for specialist) on acids
	'c0':0, #Background consumption rate in binary model
	'c1':1., #Specific consumption rate in binary model
	'fs':0.45, #Fraction of secretion flux with same resource type
	'fw':0.45, #Fraction of secretion flux to 'waste' resource
	'sparsity':0.3, #Variability in secretion fluxes among resources (must be less than 1)
	'regulation':'independent',
	'supply':'external',
	'response':'type I',
	'waste_type':1,
	'R0_food':1000, #unperturbed fixed point for supplied food
	'n_wells': n_pairs, #Number of independent wells
	'S':S, #Number of species per well
	'food':0, #index of food source
	'w':1, #energy content of resource
	'g':1, # conversion factor from energy to growth rate
	'l':0.5,#Leackage rate
	'tau':1, #timescale for esource renewal
	'm' : 1,
	'sigc':3
	}

#Generate equations
def dNdt(N,R,params):
    return MakeConsumerDynamics(a)(N,R,params)
def dRdt(N,R,params):
    return MakeResourceDynamics(a)(N,R,params)
dynamics = [dNdt,dRdt]

np.random.seed(seed)    
N0,R0 = MakeInitialState(a)
init_state = [N0,R0]
params = MakeParams(a)
N0
N0_init.index = N0.index

















M = np.sum(a['MA'])
T = len(a['MA'])
S = np.sum(a['SA'])+a['Sgen']
F = len(a['SA'])
M_waste = a['MA'][a['waste_type']]

l_row = np.concatenate((np.round(np.random.uniform(0.1, 1, a['SA'][0]), 2), np.ones(a['SA'][1])*0.01, np.zeros(a['SA'][2])))
l = pd.DataFrame(np.tile(l_row, (M,1))).transpose()
l.index = N0.index
l.columns = R0.index






l = pd.DataFrame.transpose(np.tile(l_row, (M, 1)))
l.rename(columns = R0.index)

np.tile(l_row, (M, 1))
params['l'] = np.random.uniform(0, 1, np.sum(a["SA"]))
Plate = Community(init_state, dynamics, params)
# Run to steady state. Note tol 1e-3 is the same as in Marsland et al. 
Plate.SteadyState(plot = False, tol = 1e-3, verbose = False)



params['c']
#Construct resource environments and update initial state
R0 = np.zeros(np.shape(R0))
N_init = np.zeros(np.shape(N0))
#p = [None] * n_pairs
n_pairs



# Initial N0 state
for i in range(n_pairs):
    p = numpy.random.choice(S, size=2, replace=False, p=None)
    #print(p)
    N_init[p[0], [(i*3), (i*3+1), (i*3+2)]] = [0.05, 0.5, 0.95]
    N_init[p[1], [(i*3), (i*3+1), (i*3+2)]] = [0.95, 0.5, 0.05]

N0 = pd.DataFrame(N_init, index = N0.index, columns = N0.keys())

#Construct resource environments and update initial state
init_state=[N0,R0]














