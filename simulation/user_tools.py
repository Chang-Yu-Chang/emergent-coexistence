import sys
import re
# ignore the future warning from pandas future updates
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import numpy as np
from community_simulator import *
from community_simulator.usertools import *

# Extract essential shape values
def extract_shapes(assumptions):
    #Force number of species to be an array:
    if isinstance(assumptions['MA'],numbers.Number):
        assumptions['MA'] = [assumptions['MA']]
    if isinstance(assumptions['SA'],numbers.Number):
        assumptions['SA'] = [assumptions['SA']]
    #Force numbers of species to be integers:
    assumptions['MA'] = np.asarray(assumptions['MA'], dtype=int)
    assumptions['SA'] = np.asarray(assumptions['SA'], dtype=int)
    assumptions['Sgen'] = int(assumptions['Sgen'])
    #Default waste type is last type in list:
    if 'waste_type' not in assumptions.keys():
        assumptions['waste_type']=len(assumptions['MA'])-1
    M = np.sum(assumptions['MA'])
    T = len(assumptions['MA'])
    S = np.sum(assumptions['SA'])+assumptions['Sgen']
    F = len(assumptions['SA'])
    resource_names = ['R'+str(k) for k in range(M)]
    type_names = ['T'+str(k) for k in range(T)]
    family_names = ['F'+str(k) for k in range(F)]
    consumer_names = ['S'+str(k) for k in range(S)]
    resource_index = [[type_names[m] for m in range(T) for k in range(assumptions['MA'][m])], resource_names]
    consumer_index = [[family_names[m] for m in range(F) for k in range(assumptions['SA'][m])]+['GEN' for k in range(assumptions['Sgen'])],consumer_names]
    return M, T, S, F, type_names, resource_index, family_names, consumer_index

# Overwrite the community-simulator sampling matrices
def sample_matrices(assumptions):
    # Extract total numbers of resources, consumers, resource types, and consumer families
    # Construct lists of names of resources, consumers, resource types, and consumer families
    M, T, S, F, type_names, resource_index, family_names, consumer_index = extract_shapes(assumptions)
    
    # Sample D matrix using empirical data, Specific to each functional group
    if assumptions['metabolism'] == 'empirical':
        def generalized_dirichlet(alpha,size=None): # custom Dirichlet sampling function that accepts alpha=0 (always samples 0 in those instances); use this function in case some uptake rate is 0 (e.g. in binary sampling of matrix c) because we are going to weigh secretions for each species using uptake rates
            alpha = np.asarray(alpha)
            if size == None:
                size = 1
            alpha_nonzero = np.where(alpha != 0)[0]
            out = np.zeros((size,len(alpha)))
            
            if len(alpha_nonzero)>0:
                out[:,alpha_nonzero] = dirichlet(alpha[alpha_nonzero],size=size)
            return out
        
        # Create empty list of matrices D
        D = ['NA' for i in range(S)]
        
        for s in range(S):
            DT = pd.DataFrame(np.zeros((M,M)),index=resource_index,columns=resource_index)
            # Fermenter
            if consumer_index[0][s] == 'F0':  
                for type_name in type_names:
                    MA = len(DT.loc[type_name])
                    if type_name == 'T0': #sugar
                        p = pd.Series(np.ones(M)*0/assumptions['MA'][0], index = DT.keys())
                        p.loc['T0'] = assumptions['ffss']/assumptions['MA'][0] # to sugar
                        p.loc['T1'] = assumptions['ffsa']/assumptions['MA'][1] # to acid
                        DT.loc[type_name] = generalized_dirichlet(p/assumptions['sparsity'],size=assumptions['MA'][0])
                    elif type_name == 'T1': # acid
                        p = pd.Series(np.ones(M)*0/assumptions['MA'][1], index = DT.keys())
                        p.loc['T0'] = assumptions['ffas']/assumptions['MA'][0] # to sugar
                        p.loc['T1'] = assumptions['ffaa']/assumptions['MA'][1] # to acid
                        DT.loc[type_name] = generalized_dirichlet(p/assumptions['sparsity'],size=assumptions['MA'][1])
                    else:
                        print("D matrices are hard coded. The number of resource class must be exactly 2.")
            # Respirator
            elif consumer_index[0][s] == 'F1':
                for type_name in type_names:
                    MA = len(DT.loc[type_name])
                    if type_name == 'T0': #sugar
                        p = pd.Series(np.ones(M)*0/assumptions['MA'][0], index = DT.keys())
                        p.loc['T0'] = assumptions['frss']/assumptions['MA'][0] # to sugar
                        p.loc['T1'] = assumptions['frsa']/assumptions['MA'][1] # to acid
                        DT.loc[type_name] = generalized_dirichlet(p/assumptions['sparsity'],size=assumptions['MA'][0])
                    elif type_name == 'T1': # acid
                        p = pd.Series(np.ones(M)*0/assumptions['MA'][1], index = DT.keys())
                        p.loc['T0'] = assumptions['fras']/assumptions['MA'][0] # to sugar
                        p.loc['T1'] = assumptions['fraa']/assumptions['MA'][1] # to acid
                        DT.loc[type_name] = generalized_dirichlet(p/assumptions['sparsity'],size=assumptions['MA'][1])
                    else:
                        print("D matrices are hard coded. The number of resource class must be exactly 2.")
            D[s] = DT.T

    # Sample c matrix; default by gamma
    assert assumptions['muc'] < M*assumptions['c1'], 'muc not attainable with given M and c1.'
    c = pd.DataFrame(np.zeros((S,M)), columns = resource_index, index = consumer_index)
    ## Use default. Fermenters do not have advantage on sugar, nor do respirators on acids
    if assumptions['c_symmetry'] == 'symmetry':
        for k in range(F):
            for j in range(T):
                if k==0 and j==0: # fermenter on sugar
                    c_mean = (assumptions['muc']/M)*(1+assumptions['q1'])
                    c_var = (assumptions['sigc']**2/M)*(1+assumptions['q1'])
                    thetac = c_var/c_mean
                    kc = c_mean**2/c_var
                    c.loc['F'+str(k)]['T'+str(j)] = np.random.gamma(kc,scale=thetac,size=(assumptions['SA'][k],assumptions['MA'][j]))
                elif k==1 and j==1: # respirator on acid
                    c_mean = (assumptions['muc']/M)*(1+assumptions['q2'])
                    c_var = (assumptions['sigc']**2/M)*(1+assumptions['q2'])
                    thetac = c_var/c_mean
                    kc = c_mean**2/c_var
                    c.loc['F'+str(k)]['T'+str(j)] = np.random.gamma(kc,scale=thetac,size=(assumptions['SA'][k],assumptions['MA'][j]))
                elif k==1 and j==0: # respirator on sugar
                    c_mean = (assumptions['muc']/M)*(1-assumptions['q1'])
                    c_var = (assumptions['sigc']**2/M)*(1-assumptions['q1'])
                    thetac = c_var/c_mean
                    kc = c_mean**2/c_var
                    c.loc['F'+str(k)]['T'+str(j)] = np.random.gamma(kc,scale=thetac,size=(assumptions['SA'][k],assumptions['MA'][j]))
                elif k==0 and j==1: # fermenter on acid
                    c_mean = (assumptions['muc']/M)*(1-assumptions['q2'])
                    c_var = (assumptions['sigc']**2/M)*(1-assumptions['q2'])
                    thetac = c_var/c_mean
                    kc = c_mean**2/c_var
                    c.loc['F'+str(k)]['T'+str(j)] = np.random.gamma(kc,scale=thetac,size=(assumptions['SA'][k],assumptions['MA'][j]))
    # Parameterize the c matrix using empirical data
    elif assumptions['c_symmetry'] == 'empirical':
        for k in range(F):
            for j in range(T):
                if k==0 and j==0: # fermenter on sugar
                    c_mean = assumptions['c_fs']
                    c_var = assumptions['sigc_fs']**2
                    thetac = c_var/c_mean
                    kc = c_mean**2/c_var
                    #c.loc['F'+str(k)]['T'+str(j)] = 1
                    #c.loc[('F0'), ('T0', 'R0')] = 1
                    c.loc['F'+str(k), 'T'+str(j)] = np.random.gamma(kc, scale = thetac, size = (assumptions['SA'][k], assumptions['MA'][j]))
                elif k==1 and j==1: # respirator on acid
                    c_mean = assumptions['c_ra']
                    c_var = assumptions['sigc_ra']**2
                    thetac = c_var/c_mean
                    kc = c_mean**2/c_var
                    c.loc['F'+str(k), 'T'+str(j)] = np.random.gamma(kc,scale=thetac,size=(assumptions['SA'][k],assumptions['MA'][j]))
                elif k==1 and j==0: # respirator on sugar
                    c_mean = assumptions['c_rs']
                    c_var = assumptions['sigc_rs']**2
                    thetac = c_var/c_mean
                    kc = c_mean**2/c_var
                    c.loc['F'+str(k), 'T'+str(j)] = np.random.gamma(kc,scale=thetac,size=(assumptions['SA'][k],assumptions['MA'][j]))
                elif k==0 and j==1: # fermenter on acid
                    c_mean = assumptions['c_fa']
                    c_var = assumptions['sigc_fa']**2
                    thetac = c_var/c_mean
                    kc = c_mean**2/c_var
                    c.loc['F'+str(k), 'T'+str(j)] = np.random.gamma(kc,scale=thetac,size=(assumptions['SA'][k],assumptions['MA'][j]))



    # Sample l matrix. Default by uniform
    l = pd.DataFrame(np.zeros((S,M)), columns = resource_index, index = consumer_index)
    for k in range(F):
        for j in range(T):
            if k==0 and j==0: # fermenter on sugar
                l_mean = assumptions['l1']
                l_var = assumptions['l1_sd']**2
                a_l = max(l_mean - np.sqrt(3*l_var), 0)
                b_l = min(l_mean + np.sqrt(3*l_var), 1)
                l.loc['F'+str(k), 'T'+str(j)] = np.random.uniform(low = a_l, high = b_l, size = (assumptions['SA'][k], assumptions['MA'][j]))
            elif k==1 and j==1: # respirator on acid
                l_mean = assumptions['l2']
                l_var = assumptions['l2_sd']**2
                a_l = max(l_mean - np.sqrt(3*l_var), 0)
                b_l = min(l_mean + np.sqrt(3*l_var), 1)
                l.loc['F'+str(k), 'T'+str(j)] = np.random.uniform(low = a_l, high = b_l, size = (assumptions['SA'][k], assumptions['MA'][j]))
            elif k==1 and j==0: # respirator on sugar
                l_mean = assumptions['l2']
                l_var = assumptions['l2_sd']**2
                a_l = max(l_mean - np.sqrt(3*l_var), 0)
                b_l = min(l_mean + np.sqrt(3*l_var), 1)
                l.loc['F'+str(k), 'T'+str(j)] = np.random.uniform(low = a_l, high = b_l, size = (assumptions['SA'][k], assumptions['MA'][j]))
            elif k==0 and j==1: # fermenter on acids
                l_mean = assumptions['l1']
                l_var = assumptions['l1_sd']**2
                a_l = max(l_mean - np.sqrt(3*l_var), 0)
                b_l = min(l_mean + np.sqrt(3*l_var), 1)
                l.loc['F'+str(k), 'T'+str(j)] = np.random.uniform(low = a_l, high = b_l, size = (assumptions['SA'][k], assumptions['MA'][j]))
            elif type_names[j] == waste_name:
                l.loc['F'+str(k), 'T'+str(j)] = np.ones((assumptions['SA'][k], assumptions['MA'][j]))
                
    # if 'GEN' in l.index:
    #         l_mean = assumptions['l1']
    #         l_var = assumptions['l1_sd']**2
    #         a_l = max(l_mean - np.sqrt(3*l_var), 0)
    #         b_l = min(l_mean + np.sqrt(3*l_var), 1)
    #         l.loc['GEN'] = np.random.uniform(low = a_l, high = b_l, size = (assumptions['Sgen'], M))
    
    return D, c, l

# Update the make params function to remove c, D, and l
def make_params(assumptions):
    N0,R0 = MakeInitialState(assumptions)
    
    # if not isinstance(assumptions['food'], int) or not isinstance(assumptions['R0_food'], int):
    #     params=[{'m':1,
    #             'w':1,
    #             'g':1,
    #             'R0': R0.values[:,k],
    #             'tau':1,
    #             'r':1,
    #             'sigma_max':1,
    #             'nreg':10,
    #             'n':2
    #             } for k in range(assumptions['n_wells'])]
    #     for item in ['m','w','g','tau','r','sigma_max','n','nreg']:
    #         if item in assumptions.keys():
    #             for k in range(assumptions['n_wells']):
    #                 params[k][item] = assumptions[item]
    # 
    # else:
    params={'m':1,
            'w':1,
            'g':1,
            'R0':R0.values[:,0],
            'tau':1,
            'r':1,
            'sigma_max':1,
            'nreg':10,
            'n':2
            }
        
    for item in ['m','w','g','tau','r','sigma_max','n','nreg']:
        if item in assumptions.keys():
            params[item] = assumptions[item]

    return params

# Update resource equation to include species-specific D and resource specific l
def make_resource_dynamics(assumptions):
    """
    This function is adapted so that it returns a function (dRdt) which,
    in turn, can expect params['D'] to be a list of metabolic matrices D_i
    of length equal to the total number of species (Stot).
    The i-th element of the list should the metabolic matrix of species i.
    
    Original description from the Community Simulator package:
        
    Construct resource dynamics. 'assumptions' must be a dictionary containing at least
    three entries:
    
    response = {'type I', 'type II', 'type III'} specifies nonlinearity of growth law
    
    regulation = {'independent','energy','mass'} allows microbes to adjust uptake
        rates to favor the most abundant accessible resources (measured either by
        energy or mass)
    
    supply = {'off','external','self-renewing'} sets choice of
        intrinsic resource dynamics
        
    Returns a function of N, R, and the model parameters, which itself returns the
        vector of resource rates of change dR/dt
    """
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
    
    J_in = lambda R,params: (u[assumptions['regulation']](params['c']*R,params)
                             *params['w']*sigma[assumptions['response']](R,params))
    
    """
    ### original in community-simulator package (single D matrix)
    J_out = lambda R,params: (params['l']*J_in(R,params)).dot(params['D'].T)
    ### new version --accepts D as a list of matrices (using list comprehension and trimming J_in before .dot() to speed up)
    J_out = lambda R,params: params['l']*np.array([ J_in(R,params)[i,:].dot(params['D'][i].T) for i in range(len(params['D'])) ])
    ###
    """
    
    ### new version incorporating both scenarios (D can be a matrix or a list, depending on the choice of 'metabolism')
    J_out = {'common': lambda R,params: (params['l']*J_in(R,params)).dot(params['D'].T),
             'specific': lambda R,params: params['l']*np.array([ J_in(R,params)[i,:].dot(params['D'][i].T) for i in range(len(params['D'])) ]),
             'empirical': lambda R,params: params['l']*np.array([ J_in(R,params)[i,:].dot(params['D'][i].T) for i in range(len(params['D'])) ]),
             'two-families': lambda R,params: (params['l']*J_in(R,params)).dot(params['D'].T)
            }
    
    return lambda N,R,params: (h[assumptions['supply']](R,params)
                               -(J_in(R,params)/params['w']).T.dot(N)
                               +(J_out[assumptions['metabolism']](R,params)/params['w']).T.dot(N))

# Update consumer equation to include resource specific l
def make_consumer_dynamics(assumptions):
    """
    Construct resource dynamics. 'assumptions' must be a dictionary containing at least
    three entries:
    
    response = {'type I', 'type II', 'type III'} specifies nonlinearity of growth law
    
    regulation = {'independent','energy','mass'} allows microbes to adjust uptake
        rates to favor the most abundant accessible resources (measured either by
        energy or mass)
    
    supply = {'off','external','self-renewing','predator'} sets choice of
        intrinsic resource dynamics
        
    Returns a function of N, R, and the model parameters, which itself returns the
        vector of consumer rates of change dN/dt
    """
    sigma = {'type I': lambda R,params: params['c']*R,
             'type II': lambda R,params: params['c']*R/(1+params['c']*R/params['sigma_max']),
             'type III': lambda R,params: (params['c']*R)**params['n']/(1+(params['c']*R)**params['n']/params['sigma_max'])
            }
    
    u = {'independent': lambda x,params: 1.,
         'energy': lambda x,params: (((params['w']*x)**params['nreg']).T
                                      /np.sum((params['w']*x)**params['nreg'],axis=1)).T,
         'mass': lambda x,params: ((x**params['nreg']).T/np.sum(x**params['nreg'],axis=1)).T
        }
    
    J_in = lambda R,params: (u[assumptions['regulation']](params['c']*R,params)
                             *params['w']*sigma[assumptions['response']](R,params))
    J_growth = lambda R,params: (1-params['l'])*J_in(R,params)
    
    return lambda N,R,params: params['g']*N*(np.sum(J_growth(R,params),axis=1)-params['m'])

# Load model parameters from the input csv
def load_assumptions(input_row):
    # Set parameters in model
    assumptions = a_default.copy()
    assumptions['SA'] = np.ones(2) * int(input_row['sa']) # Number of species per specialist family
    assumptions['MA'] = np.ones(2) * int(input_row['ma']) # Number of resource per resource class
    assumptions['Sgen'] = 0 # No generalist
    assumptions['n_wells'] = int(input_row['n_wells'])
    assumptions['n_communities'] = int(input_row['n_communities'])
    assumptions['m'] = 0 # Turn off maintanence cost. i.e., no cell dies
    assumptions['S'] = int(input_row['S']) # Number of initial per-well species sampled from the species pool
    
    ## Consumer consumption rates
    assumptions['muc'] = 10 # Mean sum of consumption rates (used in all models)
    assumptions['sigc'] = 3 # Standard deviation of sum of consumption rates for Gaussian and Gamma models
    assumptions['c_symmetry'] = str(input_row['c_symmetry'])
    assumptions['c_fs'] = float(input_row['c_fs'])
    assumptions['sigc_fs'] = float(input_row['sigc_fs'])
    assumptions['c_fa'] = float(input_row['c_fa'])
    assumptions['sigc_fa'] = float(input_row['sigc_fa'])
    assumptions['c_rs'] = float(input_row['c_rs'])
    assumptions['sigc_rs'] = float(input_row['sigc_rs'])
    assumptions['c_ra'] = float(input_row['c_ra'])
    assumptions['sigc_ra'] = float(input_row['sigc_ra'])

    
    ## Consumer metabolism
    assumptions['fs'] = 0.3 # Fraction of secretion flux with same resource type
    assumptions['fw'] = 0.01 # Fraction of secretion flux to 'waste' resource
    assumptions['ffss'] = float(input_row['ffss'])
    assumptions['ffsa'] = float(input_row['ffsa'])
    assumptions['ffas'] = float(input_row['ffas'])
    assumptions['ffaa'] = float(input_row['ffaa'])
    assumptions['frss'] = float(input_row['frss'])
    assumptions['frsa'] = float(input_row['frsa'])
    assumptions['fras'] = float(input_row['fras'])
    assumptions['fraa'] = float(input_row['fraa'])
    assumptions['sparsity'] = 0.1 # Effective sparsity of metabolic matrix (between 0 and 1)
    assumptions['metabolism'] = str(input_row['metabolism']) #{'common','specific'} determines whether to use a common metabolic matrix or each species having its own
    assumptions['rs'] = float(input_row['rs']) # control parameter (only used if 'metabolism' is 'specific'): if 1, each species secretes only resources that it can consume (or waste resources), preferentially those that it can consume more efficiently; if 0 secretions are randomized (default behavior of the original community-simulator package) 
    assumptions['l1'] = float(input_row['l1']) # Mean leakage rate of specialist family 1
    assumptions['l2'] = float(input_row['l2']) # Meab leakage rate of specialist family 2
    assumptions['l1_sd'] = float(input_row['l1_sd']) # SD of l1
    assumptions['l2_sd'] = float(input_row['l2_sd']) # SD of l2
    
    ## Resource
    assumptions['supply'] = 'off' # 'off' for batch culture. 'external' and 'self-renewing' for constant supply
    assumptions['tau'] = 1 # timescale for esource renewal
    assumptions['food'] = 0 # Index of the food source when a single resource is externally supplied
    assumptions['R0_food'] = float(input_row['R0_food'])
    assumptions['response'] = 'type I'
    assumptions['regulation'] = 'independent' # 'independent' is standard, 'energy' or 'mass' enable preferential consumption of resources
    
    return assumptions

# Make initial state N0 and R0 
def make_initial_state(input_row, assumptions):
    M, T, S, F, type_names, resource_index, family_names, consumer_index = extract_shapes(assumptions)
    
    # Default
    #if not isinstance(input_row['init_N0'], str) and not isinstance(input_row['init_R0'], str):
    N0,R0 = MakeInitialState(assumptions)
    # Read both N0 and R0 when provided
    #elif isinstance(input_row['init_N0'], str) and isinstance(input_row['init_R0'], str):
    N0 = pd.read_csv(input_row['output_dir'] + input_row['init_N0'])
    N0 = N0.loc[:,~N0.columns.isin(['Family', 'Species'])].astype(float)
    N0 = N0.set_index(consumer_index)
    
    R0 = pd.read_csv(input_row['output_dir'] + input_row['init_R0'])
    R0 = R0.loc[:,~R0.columns.isin(['Class', 'Resource'])].astype(float)
    R0 = R0.set_index(resource_index)
    assert N0.shape[1] == R0.shape[1], "Well numbers in the input N0 and R0 are different."
    # Update n_wells in assumptions
    #assumptions["n_wells"] = N0.shape[1]
    # # Read N0 from external csv 
    # elif isinstance(input_row['init_N0'], str) and not isinstance(input_row['init_R0'], str): 
    #     N0 = pd.read_csv(input_row['output_dir'] + input_row['init_N0'])
    #     N0 = N0.loc[:,~N0.columns.isin(['Family', 'Species'])].astype(float)
    #     #N0 = pd.read_csv(input_row['output_dir'] + input_row['init_N0'], dtype = 'float', usecols = ['W' + str(i) for i in range(n_wells)])
    #     N0 = N0.set_index(consumer_index)
    #     # Update n_wells in assumptions and create a corresponding R0
    #     #assumptions["n_wells"] = N0.shape[1]
    #     #temp, R0 = MakeInitialState(assumptions)
    # # Read R0 from external csv 
    # elif not isinstance(input_row['init_N0'], str) and isinstance(input_row['init_R0'], str): 
    #     R0 = pd.read_csv(input_row['output_dir'] + input_row['init_R0'])
    #     R0 = R0.loc[:,~R0.columns.isin(['Class', 'Resource'])].astype(float)
    #     #R0 = pd.read_csv(input_row['output_dir'] + input_row['init_R0'], dtype = 'float', usecols = ['W' + str(i) for i in range(n_wells)])
    #     R0 = R0.set_index(resource_index)
    #     # Update n_wells in assumptions and create a corresponding N0
    #     assumptions["n_wells"] = R0.shape[1]
    #     N0, temp = MakeInitialState(assumptions)
    return N0, R0

# Output matrices
def write_matrices(input_row, D, c, l):
    if input_row['metabolism'] == "common":
        D.to_csv(input_row['output_dir'] + 'D_seed' + str(input_row['seed']) + '.csv')
    elif input_row['metabolism'] == "two-families":
        D.to_csv(input_row['output_dir'] + 'D_seed' + str(input_row['seed']) + '.csv')
    elif input_row['metabolism'] == "specific":
        D[0].to_csv(input_row['output_dir'] + 'D_seed' + str(input_row['seed']) + '.csv')
    elif input_row['metabolism'] == "empirical":
        S = np.sum(input_row['sa'] * 2)
        # Load only one species'D
        D[0].to_csv(input_row['output_dir'] + 'D_S' + str(0) +  '_seed' + str(input_row['seed']) + '.csv')
        # Uncomment to load all species
        # for s in range(S):
        #     D[s].to_csv(input_row['output_dir'] + 'D_S' + str(s) +  '_seed' + str(input_row['seed']) + '.csv')
    
    c.to_csv(input_row['output_dir'] + 'c_seed' + str(input_row['seed']) + '.csv')
    l.to_csv(input_row['output_dir'] + 'l_seed' + str(input_row['seed']) + '.csv')

# Run the simulation using the customized parameter sets
def run_simulations(input_row):
    seed = int(input_row['seed'])
    exp_id = int(input_row['exp_id'])

    # Update model parameters
    a = load_assumptions(input_row)
    assumptions = a
    
    # Sample matrices
    np.random.seed(seed)
    D, c, l = sample_matrices(a)
    write_matrices(input_row, D, c, l)
    print("Matrices written")

    # Update params 
    params = make_params(a)
    params['c'] = c
    params['D'] = D
    params['l'] = l
    #print(params)
    # print(type(params))
    # for k in range(assumptions['n_wells']):
    #     params[k]['c'] = c
    #     params[k]['D'] = D
    #     params[k]['l'] = l

    
    # Generate equations
    def dNdt(N,R,params):
        return make_consumer_dynamics(a)(N,R,params) # Use the updated version to include varied l
    def dRdt(N,R,params):
        return make_resource_dynamics(a)(N,R,params) # Use the updated version to include varied D and l
    dynamics = [dNdt,dRdt]

    # Set initial state by reading in plate conditions
    N0, R0 = make_initial_state(input_row, a)
    init_state = [N0,R0]
    #print(N0.shape)
    
    # Make plate object
    Plate = Community(init_state, dynamics, params, parallel = False)
#    if input_row['save_timepoint']:
    for i in range(12): # number of passages
        for j in range(10): # time of propagation 
            Plate.Propagate(T = 1)
            print("T" + str(i+1) + "t" + str(j+1))
            # if j == 0:
            #     Plate.N.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "T" + str(i+1) + "t" + str(0) + ".csv", input_row["init_N0"]))
            #     Plate.R.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "T" + str(i+1) + "t" + str(0) + ".csv", input_row["init_R0"]))
        #Plate.N.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "T" + str(i+1) + "t" + str(j+1) + ".csv", input_row["init_N0"]))
        Plate.N.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "T" + str(i+1) + ".csv", input_row["init_N0"]))
        Plate.R.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "T" + str(i+1) + ".csv", input_row["init_R0"]))
        if i == 12: # last transfer
            Plate.N.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "end.csv", input_row["init_N0"]))
            Plate.R.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "end.csv", input_row["init_R0"]))

        Plate.Passage(f = np.eye(a['n_wells'])/10, refresh_resource = True)

    # elif input_row['save_timepoint'] == False:
    #     Plate.RunExperiment(np.eye(a['n_wells'])/10, T = 1, npass = 1, refresh_resource = True)
    #     Plate.N.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "end.csv", input_row["init_N0"]))
    #     Plate.R.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "end.csv", input_row["init_R0"]))
    #     Plate.R.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "Rend.csv", input_row["init_N0"]))

# For common metabolism, use T = 10, npass=5, muc1 = 10, muc2 = 20
# For two-families metabolism, use T = 10, npass=5, muc1 = 10, muc2 = 10





