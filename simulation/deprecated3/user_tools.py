import sys
import re
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
    return M, T, S, F, type_names, resource_index, consumer_index

# Overwrite the community-simulator sampling matrices
def sample_matrices(assumptions):
    #PREPARE VARIABLES
    #Extract total numbers of resources, consumers, resource types, and consumer families:
    #Construct lists of names of resources, consumers, resource types, and consumer families:
    M, T, S, F, type_names, resource_index, consumer_index = extract_shapes(assumptions)
    M_waste = assumptions['MA'][assumptions['waste_type']]
    waste_name = type_names[assumptions['waste_type']]
    
    # Sample D matrix
    if assumptions['metabolism'] == 'common': # if the metabolic matrix is common for all species, continue normally as in the original community-simulator package
        # SAMPLE METABOLIC MATRIX FROM DIRICHLET DISTRIBUTION
        DT = pd.DataFrame(np.zeros((M,M)),index=resource_index,columns=resource_index)
        for type_name in type_names:
            MA = len(DT.loc[type_name])
            if type_name is not waste_name:
                #Set background secretion levels
                p = pd.Series(np.ones(M)*(1-assumptions['fs']-assumptions['fw'])/(M-MA-M_waste),index = DT.keys())
                #Set self-secretion level
                p.loc[type_name] = assumptions['fs']/MA
                #Set waste secretion level
                p.loc[waste_name] = assumptions['fw']/M_waste
                #Sample from dirichlet
                DT.loc[type_name] = dirichlet(p/assumptions['sparsity'],size=MA)
            else:
                if M > MA:
                    #Set background secretion levels
                    p = pd.Series(np.ones(M)*(1-assumptions['fw']-assumptions['fs'])/(M-MA),index = DT.keys())
                    #Set self-secretion level
                    p.loc[type_name] = (assumptions['fw']+assumptions['fs'])/MA
                else:
                    p = pd.Series(np.ones(M)/M,index = DT.keys())
                #Sample from dirichlet
                DT.loc[type_name] = dirichlet(p/assumptions['sparsity'],size=MA)
        D = DT.T
    elif assumptions['metabolism'] == 'two-families':
        DT = pd.DataFrame(np.zeros((M,M)),index=resource_index,columns=resource_index)
        for type_name in type_names:
            MA = len(DT.loc[type_name])
            if type_name == 'T0': #sugar
                p = pd.Series(np.ones(M)*assumptions['fsw']/assumptions['MA'][2], index = DT.keys())
                p.loc['T0'] = assumptions['fss']/assumptions['MA'][0] # to sugar
                p.loc['T1'] = assumptions['fsa']/assumptions['MA'][1] # to acid
                p.loc['T2'] = assumptions['fsw']/assumptions['MA'][2] # to waste
                DT.loc[type_name] = dirichlet(p/assumptions['sparsity'],size=assumptions['MA'][0])
            elif type_name == 'T1': # acid
                p = pd.Series(np.ones(M)*assumptions['fsw']/assumptions['MA'][2], index = DT.keys())
                p.loc['T0'] = assumptions['fas']/assumptions['MA'][0] # to sugar
                p.loc['T1'] = assumptions['faa']/assumptions['MA'][1] # to acid
                p.loc['T2'] = assumptions['faw']/assumptions['MA'][2] # to waste
                DT.loc[type_name] = dirichlet(p/assumptions['sparsity'],size=assumptions['MA'][1])
            elif type_name == 'T2': # waste
                p = pd.Series(np.ones(M)*assumptions['fsw']/assumptions['MA'][2], index = DT.keys())
                p.loc['T0'] = assumptions['fws']/assumptions['MA'][0] # to sugar
                p.loc['T1'] = assumptions['fwa']/assumptions['MA'][1] # to acid
                p.loc['T2'] = assumptions['fww']/assumptions['MA'][2] # to waste
                DT.loc[type_name] = dirichlet(p/assumptions['sparsity'],size=assumptions['MA'][2])
            else:
                print("D matrices are hard coded. The number of resource class must be exactly 3.")
        D = DT.T
    elif assumptions['metabolism'] == 'specific': # if we want a different metabolic matrix for each species, generate a list
        
        def generalized_dirichlet(alpha,size=None): # custom Dirichlet sampling function that accepts alpha=0 (always samples 0 in those instances); use this function in case some uptake rate is 0 (e.g. in binary sampling of matrix c) because we are going to weigh secretions for each species using uptake rates
            alpha = np.asarray(alpha)
            if size == None:
                size = 1
                
            alpha_nonzero = np.where(alpha != 0)[0]
            out = np.zeros((size,len(alpha)))
            
            if len(alpha_nonzero)>0:
                out[:,alpha_nonzero] = dirichlet(alpha[alpha_nonzero],size=size)
            
            return out
    
        # create empty list of matrices D
        D = ['NA' for i in range(S)]
        
        # loop through species and generate matrices
        for s in range(S):
            
            # initialize matrix for species s
            DT = pd.DataFrame(np.zeros((M,M)),index=resource_index,columns=resource_index)
            
            # weights
            cs = c.iloc[s,:] # uptake rates of species s (will be used to weigh secretions based on assumptions['rs'])
            weight = (cs - 1)*assumptions['rs'] + 1 # each species secretion levels will depend on the affinity (uptake rate) of that species for the secreted resource (controlled by parameter 'rs' in the assumptions)
            weight.loc[waste_name] = 1 # weighing does not apply to waste resource (waste secretions remain random)
            # normalize weights so that final fluxes (determined by fs and fw) are not altered
            wsum = weight.copy()
            for type_name in type_names:
                MA = len(DT.loc[type_name])
                if weight.loc[type_name].sum()>0:
                    wsum.loc[type_name] = weight.loc[type_name].sum()/MA
                else:
                    wsum.loc[type_name] = 1 # if there is no flux to an entire resource type, this would give a 0 and then a NaN when dividing. Instead, we set it to 1 so it gives 0 when dividing (this will modify energy fluxes defined by fs and fw but there is no way around it)
            weight = weight.div(wsum)
            
            # make matrix
            for type_name in type_names:
                MA = len(DT.loc[type_name])
                if type_name is not waste_name:
                    #Set background secretion levels
                    p = pd.Series(np.ones(M)*(1-assumptions['fs']-assumptions['fw'])/(M-MA-M_waste),index = DT.keys())
                    #Set self-secretion level
                    p.loc[type_name] = assumptions['fs']/MA
                    #Set waste secretion level
                    p.loc[waste_name] = assumptions['fw']/M_waste
                    # new: add weights
                    p = p*weight
                    #Sample from dirichlet
                    DT.loc[type_name] = generalized_dirichlet(p/assumptions['sparsity'],size=MA) # use generalized Dirichlet in case some weighs are zero (in the limit case where an uptake rate is zero, e.g. if matrix c is binary)
                else: # the waste resource is treated as before, i.e. there is no weighing: all species can produce arbitrary amounts of waste resources regardless of their ability to consume them
                    if M > MA:
                        #Set background secretion levels
                        p = pd.Series(np.ones(M)*(1-assumptions['fw']-assumptions['fs'])/(M-MA),index = DT.keys())
                        #Set self-secretion level
                        p.loc[type_name] = (assumptions['fw']+assumptions['fs'])/MA
                    else:
                        p = pd.Series(np.ones(M)/M,index = DT.keys())
                    #Sample from dirichlet
                    DT.loc[type_name] = dirichlet(p/assumptions['sparsity'],size=MA)
            D[s] = DT.T


    
    # Sample c matrix; default by gamma
    assert assumptions['muc'] < M*assumptions['c1'], 'muc not attainable with given M and c1.'
    c = pd.DataFrame(np.zeros((S,M)),columns=resource_index,index=consumer_index)
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
    ## When fermenters have advantage on sugar, and so do respirators on acids
    ## Now only coded so that respirator overall is a good user
    elif assumptions['c_symmetry'] == 'asymmetry':
        for k in range(F):
            for j in range(T):
                if k==0 and j==0: # fermenter on sugar
                    c_mean = (assumptions['muc1']/M)*(1+assumptions['q1'])
                    c_var = (assumptions['sigc']**2/M)*(1+assumptions['q1'])
                    thetac = c_var/c_mean
                    kc = c_mean**2/c_var
                    c.loc['F'+str(k)]['T'+str(j)] = np.random.gamma(kc,scale=thetac,size=(assumptions['SA'][k],assumptions['MA'][j]))
                elif k==1 and j==1: # respirator on acid
                    c_mean = (assumptions['muc2']/M)*(1+assumptions['q2'])
                    c_var = (assumptions['sigc']**2/M)*(1+assumptions['q2'])
                    thetac = c_var/c_mean
                    kc = c_mean**2/c_var
                    c.loc['F'+str(k)]['T'+str(j)] = np.random.gamma(kc,scale=thetac,size=(assumptions['SA'][k],assumptions['MA'][j]))
                elif k==1 and j==0: # respirator on sugar
                    c_mean = (assumptions['muc2']/M)*(1-assumptions['q1'])
                    c_var = (assumptions['sigc']**2/M)*(1-assumptions['q1'])
                    thetac = c_var/c_mean
                    kc = c_mean**2/c_var
                    c.loc['F'+str(k)]['T'+str(j)] = np.random.gamma(kc,scale=thetac,size=(assumptions['SA'][k],assumptions['MA'][j]))
                elif k==0 and j==1: # fermenter on acid
                    c_mean = (assumptions['muc1']/M)*(1-assumptions['q2'])
                    c_var = (assumptions['sigc']**2/M)*(1-assumptions['q2'])
                    thetac = c_var/c_mean
                    kc = c_mean**2/c_var
                    c.loc['F'+str(k)]['T'+str(j)] = np.random.gamma(kc,scale=thetac,size=(assumptions['SA'][k],assumptions['MA'][j]))

    # Sample l matrix. Default by uniform
    l = pd.DataFrame(np.zeros((S,M)),columns=resource_index,index=consumer_index)
    for k in range(F):
        for j in range(T):
            if k==0 and j==0: # fermenter on sugar
                l_mean = assumptions['l1']
                l_var = assumptions['l1_var']
                a_l = max(l_mean - np.sqrt(3*l_var), 0)
                b_l = min(l_mean + np.sqrt(3*l_var), 1)
                l.loc['F'+str(k)]['T'+str(j)] = np.random.uniform(low = a_l, high = b_l, size = (assumptions['SA'][k], assumptions['MA'][j]))
            elif k==1 and j==1: # respirator on acid
                l_mean = assumptions['l2']
                l_var = assumptions['l2_var']
                a_l = max(l_mean - np.sqrt(3*l_var), 0)
                b_l = min(l_mean + np.sqrt(3*l_var), 1)
                l.loc['F'+str(k)]['T'+str(j)] = np.random.uniform(low = a_l, high = b_l, size = (assumptions['SA'][k], assumptions['MA'][j]))
            elif k==1 and j==0: # respirator on sugar
                l_mean = assumptions['l2']
                l_var = assumptions['l2_var']
                a_l = max(l_mean - np.sqrt(3*l_var), 0)
                b_l = min(l_mean + np.sqrt(3*l_var), 1)
                l.loc['F'+str(k)]['T'+str(j)] = np.random.uniform(low = a_l, high = b_l, size = (assumptions['SA'][k], assumptions['MA'][j]))
            elif k==0 and j==1: # fermenter on acids
                l_mean = assumptions['l1']
                l_var = assumptions['l1_var']
                a_l = max(l_mean - np.sqrt(3*l_var), 0)
                b_l = min(l_mean + np.sqrt(3*l_var), 1)
                l.loc['F'+str(k)]['T'+str(j)] = np.random.uniform(low = a_l, high = b_l, size = (assumptions['SA'][k], assumptions['MA'][j]))
            elif type_names[j] == waste_name:
                l.loc['F'+str(k)]['T'+str(j)] = np.ones((assumptions['SA'][k], assumptions['MA'][j]))
    if 'GEN' in l.index:
            l_mean = assumptions['l1']
            l_var = assumptions['l1_var']
            a_l = max(l_mean - np.sqrt(3*l_var), 0)
            b_l = min(l_mean + np.sqrt(3*l_var), 1)
            l.loc['GEN'] = np.random.uniform(low = a_l, high = b_l, size = (assumptions['Sgen'], M))
    
    return D, c, l

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

# Load model parameters
def load_assumptions(input_row):
    # Load Parameters from the input csv
    output_dir = input_row['output_dir']
    save_timepoint = input_row['save_timepoint']
    sa = int(input_row['sa'])
    ma = int(input_row['ma'])
    S = int(input_row['S'])
    q1 = float(input_row['q1'])
    q2 = float(input_row['q2'])
    l1 = float(input_row['l1'])
    l2 = float(input_row['l2'])
    l1_var = float(input_row['l1_var'])
    l2_var = float(input_row['l2_var'])
    c_symmetry = str(input_row['c_symmetry'])
    muc1 = int(input_row['muc1'])
    muc2 = int(input_row['muc2'])
    n_communities = int(input_row['n_communities'])
    n_wells = int(input_row['n_wells'])
    metabolism = str(input_row['metabolism'])
    rs = float(input_row['rs'])

    # Set parameters in model
    assumptions = a_default.copy()
    assumptions['SA'] = np.ones(2) * sa # Number of species per specialist family
    assumptions['MA'] = np.ones(3) * ma # Number of resource per resource class
    assumptions['Sgen'] = 0 # No generalist
    assumptions['n_wells'] = n_wells
    assumptions['n_communities'] = n_communities
    assumptions['m'] = 0 # Turn off maintanence cost. i.e., no cell dies
    assumptions['S'] = S # Number of initial per-well species sampled from the species pool
    
    ## Consumer consumption rates
    #assumptions['sampling'] = 'Gamma' # 'Binary', 'Gaussian', 'Uniform' or 'Gamma' (sampling of the matrix c)
    assumptions['muc'] = 10 # Mean sum of consumption rates (used in all models)
    assumptions['sigc'] = 3 # Standard deviation of sum of consumption rates for Gaussian and Gamma models
    #assumptions['c0'] = 0 # Sum of background consumption rates in binary model
    #assumptions['c1'] = 1 # Specific consumption rate in binary model
    assumptions['q1'] = q1 # Preference strengeth of specialist family 1 (set to 0 for generalist and 1 for specialist)
    assumptions['q2'] = q2 # Preference strengeth of specialist family 2 
    assumptions['c_symmetry'] = c_symmetry
    assumptions['muc1'] = muc1
    assumptions['muc2'] = muc2

    
    ## Consumer metabolism
    assumptions['fs'] = 0.3 # Fraction of secretion flux with same resource type
    assumptions['fw'] = 0.01 # Fraction of secretion flux to 'waste' resource
    assumptions['fss'] = 0.10 # sugar to sugar
    assumptions['fsa'] = 0.80 # sugar to acid
    assumptions['fsw'] = 0.10 # sugar to waste
    assumptions['fas'] = 0.001 # acid to sugar
    assumptions['faa'] = 0.10 # acid to acid
    assumptions['faw'] = 0.899 # acid to waste
    assumptions['fws'] = 0.001 # waste to sugar
    assumptions['fwa'] = 0.001 # waste to acid
    assumptions['fww'] = 0.998 # waste to waste
    assumptions['sparsity'] = 0.1 # Effective sparsity of metabolic matrix (between 0 and 1)
    assumptions['metabolism'] = metabolism #{'common','specific'} determines whether to use a common metabolic matrix or each species having its own
    assumptions['rs'] = rs # control parameter (only used if 'metabolism' is 'specific'): if 1, each species secretes only resources that it can consume (or waste resources), preferentially those that it can consume more efficiently; if 0 secretions are randomized (default behavior of the original community-simulator package) 
    assumptions['l1'] = l1 # Mean leakage rate of specialist family 1
    assumptions['l2'] = l2 # Meab leakage rate of specialist family 2
    assumptions['l1_var'] = l1_var # Variance of l1
    assumptions['l2_var'] = l2_var # Variance of l2
    
    ## Resource
    assumptions['supply'] = 'off' # 'off' for batch culture. 'external' and 'self-renewing' for constant supply
    assumptions['tau'] = 1 # timescale for esource renewal
    assumptions['food'] = 0 # Index of the food source when a single resource is externally supplied
    assumptions['R0_food'] = 1000
    assumptions['response'] = 'type I'
    assumptions['regulation'] = 'independent' # 'independent' is standard, 'energy' or 'mass' enable preferential consumption of resources
    
    return assumptions

# Make initial state N0 and R0 
def make_initial_state(input_row, assumptions):
    M, T, S, F, type_names, resource_index, consumer_index = extract_shapes(assumptions)
    # Default
    if not isinstance(input_row['init_N0'], str) and not isinstance(input_row['init_R0'], str):
        N0,R0 = MakeInitialState(assumptions)
    # Read N0 from external csv 
    elif isinstance(input_row['init_N0'], str) and not isinstance(input_row['init_R0'], str): 
        N0 = pd.read_csv(input_row['output_dir'] + input_row['init_N0'])
        N0 = N0.loc[:,~N0.columns.isin(['Family', 'Species'])].astype(float)
        #N0 = pd.read_csv(input_row['output_dir'] + input_row['init_N0'], dtype = 'float', usecols = ['W' + str(i) for i in range(n_wells)])
        N0 = N0.set_index(consumer_index)
        # Update n_wells in assumptions and create a corresponding R0
        assumptions["n_wells"] = N0.shape[1]
        temp, R0 = MakeInitialState(assumptions)
    # Read R0 from external csv 
    elif not isinstance(input_row['init_N0'], str) and isinstance(input_row['init_R0'], str): 
        R0 = pd.read_csv(input_row['output_dir'] + input_row['init_N0'])
        R0 = R0.loc[:,~R0.columns.isin(['Class', 'Resource'])].astype(float)
        #R0 = pd.read_csv(input_row['output_dir'] + input_row['init_R0'], dtype = 'float', usecols = ['W' + str(i) for i in range(n_wells)])
        R0 = R0.set_index(resource_index)
        # Update n_wells in assumptions and create a corresponding N0
        assumptions["n_wells"] = R0.shape[1]
        N0, temp = MakeInitialState(assumptions)
    # Read both N0 and R0 when provided
    elif isinstance(input_row['init_N0'], str) and isinstance(input_row['init_R0'], str):
        N0 = pd.read_csv(input_row['output_dir'] + input_row['init_N0'])
        N0 = N0.loc[:,~N0.columns.isin(['Family', 'Species'])].astype(float)
        N0 = N0.set_index(consumer_index)
        R0 = pd.read_csv(input_row['output_dir'] + input_row['init_N0'])
        R0 = R0.loc[:,~R0.columns.isin(['Class', 'Resource'])].astype(float)
        R0 = R0.set_index(resource_index)
        assert N0.shape[1] == R0.shape[1], "Well numbers in the input N0 and R0 are different."
        # Update n_wells in assumptions 
        assumptions["n_wells"] = N0.shape[1]
    return N0, R0

# Output matrices
def write_matrices(input_row, D, c, l):
    if input_row['metabolism'] == "common":
        D.to_csv(input_row['output_dir'] + 'D_seed' + str(input_row['seed']) + '.csv')
    elif input_row['metabolism'] == "two-families":
        D.to_csv(input_row['output_dir'] + 'D_seed' + str(input_row['seed']) + '.csv')
    elif input_row['metabolism'] == "specific":
        D[0].to_csv(input_row['output_dir'] + 'D_seed' + str(input_row['seed']) + '.csv')
    
    c.to_csv(input_row['output_dir'] + 'c_seed' + str(input_row['seed']) + '.csv')
    l.to_csv(input_row['output_dir'] + 'l_seed' + str(input_row['seed']) + '.csv')

# Run the simulation using the customized parameter sets
def run_simulations(input_row):
    seed = int(input_row['seed'])
    exp_id = int(input_row['exp_id'])

    # Update model parameters
    a = load_assumptions(input_row)
    
    # Sample matrices
    np.random.seed(seed)
    D, c, l = sample_matrices(a)
    write_matrices(input_row, D, c, l)
    print("Matrices written")

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
    if input_row['save_timepoint']:
        for i in range(5): # number of passages
            for j in range(10): # time of propagation 
                Plate.Propagate(T = 1)
                Plate.N.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "T" + str(i+1) + "t" + str(j+1) + ".csv", input_row["init_N0"]))
                Plate.Passage(f = np.eye(a['n_wells'])/10)
                print("T" + str(i+1) + "t" + str(j+1))
        
        Plate.N.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "end.csv", input_row["init_N0"]))
    
    elif input_row['save_timepoint'] == False:
        Plate.RunExperiment(np.eye(a['n_wells'])/10, T = 10, npass = 5, refresh_resource=True)
        Plate.N.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "end.csv", input_row["init_N0"]))


# For common metabolism, use T = 10, npass=5, muc1 = 10, muc2 = 20
# For two-families metabolism, use T = 10, npass=5, muc1 = 10, muc2 = 10












