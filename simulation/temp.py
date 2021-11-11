import sys
import pandas as pd
import numpy as np
from community_simulator import *
from community_simulator.usertools import *
from community_simulator.essentialtools import *


input_row = pd.read_csv("input_poolPairs.csv").loc[0]

### Load Parameters from the input csv
output_dir = input_row['output_dir']
q = input_row['q']
q2 = input_row['q2']
seed = int(input_row['seed'])
exp_id = int(input_row['exp_id'])
vamp = float(input_row['vamp'])
S = int(input_row['S'])
n_pairs = int(input_row["n_pairs"])
np.random.seed(seed)

# 
#N0_init = pd.read_csv(output_dir + 'communityPairs_' + str(exp_id) + '_init.csv')
#n_pairs = int(N0_init.shape[1])

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
	'sigc':3 * vamp
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

#params["l"] = params["c"].loc[:,("T0","R0")]
params["l"] = params["c"]

Plate = Community(init_state, dynamics, params, parallel = False)
Plate.Propagate(24)




def sample_l_matrix(assumptions):
    #Force number of species to be an array:
    if isinstance(assumptions['MA'],numbers.Number):
        assumptions['MA'] = [assumptions['MA']]
    if isinstance(assumptions['SA'],numbers.Number):
        assumptions['SA'] = [assumptions['SA']]
    #Force numbers of species to be integers:
    assumptions['MA'] = np.asarray(assumptions['MA'],dtype=int)
    assumptions['SA'] = np.asarray(assumptions['SA'],dtype=int)
    assumptions['Sgen'] = int(assumptions['Sgen'])
    #Default waste type is last type in list:
    if 'waste_type' not in assumptions.keys():
        assumptions['waste_type']=len(assumptions['MA'])-1

    #Extract total numbers of resources, consumers, resource types, and consumer families:
    M = np.sum(assumptions['MA'])
    T = len(assumptions['MA'])
    S = np.sum(assumptions['SA'])+assumptions['Sgen']
    F = len(assumptions['SA'])
    M_waste = assumptions['MA'][assumptions['waste_type']]
    #Construct lists of names of resources, consumers, resource types, and consumer families:
    resource_names = ['R'+str(k) for k in range(M)]
    type_names = ['T'+str(k) for k in range(T)]
    family_names = ['F'+str(k) for k in range(F)]
    consumer_names = ['S'+str(k) for k in range(S)]
    waste_name = type_names[assumptions['waste_type']]
    resource_index = [[type_names[m] for m in range(T) for k in range(assumptions['MA'][m])],
                      resource_names]
    consumer_index = [[family_names[m] for m in range(F) for k in range(assumptions['SA'][m])]
                      +['GEN' for k in range(assumptions['Sgen'])],consumer_names]
    #
    l_row = np.concatenate((np.round(np.random.uniform(0.1, 1, assumptions['SA'][0]), 2), np.ones(assumptions['SA'][1])*0.01, np.zeros(assumptions['SA'][2])))
    l = pd.DataFrame(np.tile(l_row, (M,1))).transpose()
    l.index = consumer_index
    l.columns = resource_index
    
    return(l)

sample_l_matrix(a)














assumptions = a.copy()
 #Extract total numbers of resources, consumers, resource types, and consumer families:
M = np.sum(assumptions['MA'])
T = len(assumptions['MA'])
S = np.sum(assumptions['SA'])+assumptions['Sgen']
F = len(assumptions['SA'])
M_waste = assumptions['MA'][assumptions['waste_type']]
#Construct lists of names of resources, consumers, resource types, and consumer families:
resource_names = ['R'+str(k) for k in range(M)]
type_names = ['T'+str(k) for k in range(T)]
family_names = ['F'+str(k) for k in range(F)]
consumer_names = ['S'+str(k) for k in range(S)]
waste_name = type_names[assumptions['waste_type']]
resource_index = [[type_names[m] for m in range(T) for k in range(assumptions['MA'][m])],
                  resource_names]
consumer_index = [[family_names[m] for m in range(F) for k in range(assumptions['SA'][m])]
                  +['GEN' for k in range(assumptions['Sgen'])],consumer_names]
#
l_row = np.concatenate((np.round(np.random.uniform(0.1, 1, assumptions['SA'][0]), 2), np.ones(assumptions['SA'][1])*0.01, np.zeros(assumptions['SA'][2])))
l = pd.DataFrame(np.tile(l_row, (M,1))).transpose()
l.index = consumer_index
l.columns = resource_index


y_in = Plate.N.append(Plate.R).values

#PACKAGE SYSTEM STATE AND PARAMETERS IN LIST OF DICTIONARIES
if isinstance(Plate.params,list):
    well_info = [{'y0':y_in[:,k],'params':Plate.params[k]} for k in range(Plate.n_wells)]
else:
    well_info = [{'y0':y_in[:,k],'params':Plate.params} for k in range(Plate.n_wells)]

well_info[0]["y0"]

T0=0;T=1;ns=2
t = np.linspace(T0,T0+T,ns)
y0 = well_info[0]["y0"]
S = well_info[0]['params']['S']
M = len(y0)-S
not_extinct = y0>0
compress_species = True
compress_resources = True
if not compress_species:
    not_extinct[:S] = True
if not compress_resources:  #only compress resources if we're running non-renewable dynamics
    not_extinct[S:] = True
S_comp = np.sum(not_extinct[:S]) #record the new point dividing species from resources
not_extinct_idx = np.where(not_extinct)[0]
y0_comp = y0[not_extinct]
not_extinct_consumers = not_extinct[:S]
not_extinct_resources = not_extinct[S:]

"Think the probkem is that the SXM dim l is not compressed"

if 'SxM' in dimensions.keys():
    for item in dimensions['SxM']:
        if item in params_comp.keys():
            assert np.shape(params_comp[item])==(S,M), 'Invalid shape for ' + item + '. Please update dimensions dictionary with correct dimensions.'
            params_comp[item]=params_comp[item][not_extinct_consumers,:]
            params_comp[item]=params_comp[item][:,not_extinct_resources]

params_comp = CompressParams(not_extinct_consumers,not_extinct_resources,well_info[0]['params'],Plate.dimensions,S,M)
out = integrate.odeint(Plate.dydt,y0_comp,t,args=(params_comp,S_comp),mxstep=10000,atol=1e-4)
traj = np.zeros((np.shape(out)[0],S+M))
traj[:,not_extinct_idx] = out


#PREPARE INTEGRATOR FOR PARALLEL PROCESSING
IntegrateTheseWells = partial(IntegrateWell,Plate,T=T,compress_resources=True,
                              compress_species=True)

#INITIALIZE PARALLEL POOL AND SEND EACH WELL TO ITS OWN WORKER
y_out = np.asarray(list(map(IntegrateTheseWells, well_info))).squeeze().T

#HANDLE CASE OF A SINGLE-WELL PLATE
if len(np.shape(y_out)) == 1:
    y_out = y_out[:,np.newaxis]

#UPDATE STATE VARIABLES WITH RESULTS OF INTEGRATION
self.N = pd.DataFrame(y_out[:self.S,:],
                      index = self.N.index, columns = self.N.keys())
self.R = pd.DataFrame(y_out[self.S:,:],
                      index = self.R.index, columns = self.R.keys())


















range_expand.RunExperiment(np.eye(20)/600, T = 24, npass = 30,refresh_resource=True)
Plate.N
#Plate.RunExperiment():
# Run to steady state. Note tol 1e-3 is the same as in Marsland et al. 
Plate.SteadyState(plot = False, tol = 1e-3, verbose = False)





if False:
    M = np.sum(a['MA'])
    T = len(a['MA'])
    S = np.sum(a['SA'])+a['Sgen']
    F = len(a['SA'])
    M_waste = a['MA'][a['waste_type']]
    l_row = np.concatenate((np.round(np.random.uniform(0.1, 1, a['SA'][0]), 2), np.ones(a['SA'][1])*0.01, np.zeros(a['SA'][2])))
    l = pd.DataFrame(np.tile(l_row, (M,1))).transpose()
    l.index = N0.index
    l.columns = R0.index
    #l = pd.DataFrame.transpose(np.tile(l_row, (M, 1)))
    #l.rename(columns = R0.index)
    #np.tile(l_row, (M, 1))
    #params['l'] = np.random.uniform(0, 1, np.sum(a["SA"]))
    params["l"] = l

#w = 1
#params['c']*params['w']*(1-params['l'])/w 

f0 = 1
well_name = Plate.N.keys()[0]
N_well = Plate.N.copy()[well_name] * f0
R_well = Plate.R.copy()[well_name]

y0 = well_info['y0']
S = well_info['params']['S']
M = len(y0)-S
not_extinct = y0>0


y_in = Plate.N.append(Plate.R).values
len(y_in[:2])
y_in[:2]
y_in[2:]

Plate.N.values
R = np.array(R0)
R = R0["W0"]
N = N0["W0"]
params["c"]* R

sigma = {'type I': lambda R,params: params['c']*R,
         'type II': lambda R,params: params['c']*R/(1+params['c']*R/params['sigma_max']),
         'type III': lambda R,params: (params['c']*R)**params['n']/(1+(params['c']*R)**params['n']/params['sigma_max'])
    }

u = {'independent': lambda x,params: 1.,
     'energy': lambda x,params: (((params['w']*x)**params['nreg']).T
                                  /np.sum((params['w']*x)**params['nreg'],axis=1)).T,
     'mass': lambda x,params: ((x**params['nreg']).T/np.sum(x**params['nreg'],axis=1)).T
    }

h = {'off': lambda R,params: 0.,
     'external': lambda R,params: (params['R0']-R)/params['tau'],
     'self-renewing': lambda R,params: params['r']*R*(params['R0']-R),
     'predator': lambda R,params: params['r']*R*(params['R0']-R)-params['u']*R
}
assumptions = a.copy()
J_in = lambda R,params: (u[assumptions['regulation']](params['c']*R,params)
                         *params['w']*sigma[assumptions['response']](R,params))
J_out = lambda R,params: (params['l']*J_in(R,params)).dot(params['D'].T)

params["l"] = params["c"].iloc[:,0]
J_in(R,params)

lambda N,R,params: (h[assumptions['supply']](R,params) - (J_in(R,params)/params['w']).T.dot(N) + (J_out(R,params)/params['w']).T.dot(N))


R = init_state[1]
R = Plate.R.values
params['c']*R


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














