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
    #Extract total numbers of resources, consumers, resource types, and consumer families:
    M = np.sum(assumptions['MA'])                        # Total number of resources
    T = len(assumptions['MA'])                           # Number of resource types 
    S = np.sum(assumptions['SA'])+assumptions['Sgen']    # Total number of species
    F = len(assumptions['SA'])                           # Number of species families
    M_waste = assumptions['MA'][assumptions['waste_type']]
    #Construct lists of names of resources, consumers, resource types, and consumer families:
    resource_names = ['R'+str(k) for k in range(M)]
    type_names = ['T'+str(k) for k in range(T)]
    family_names = ['F'+str(k) for k in range(F)]
    consumer_names = ['S'+str(k) for k in range(S)]
    waste_name = type_names[assumptions['waste_type']]
    resource_index = [[type_names[m] for m in range(T) for k in range(assumptions['MA'][m])], resource_names]
    consumer_index = [[family_names[m] for m in range(F) for k in range(assumptions['SA'][m])]+['GEN' for k in range(assumptions['Sgen'])],consumer_names]
    return M, T, S, F, M_waste, type_names, family_names, waste_name, resource_index, family_names, consumer_index

# Load model parameters from the input csv
def load_assumptions(input_row):
    # Set parameters in model
    assumptions = a_default.copy()
    assumptions['SA'] = np.ones(int(input_row['fa']), dtype = int) * int(input_row['sa']) # Number of species per specialist family
    assumptions['MA'] = np.ones(int(input_row['fa']), dtype = int) * int(input_row['ma']) # Number of resource per resource class
    assumptions['Sgen'] = int(input_row['Sgen'])
    assumptions['n_wells'] = int(input_row['n_wells'])
    
    #
    for item in ['ma','sa']:
        if item in input_row.keys():
            assumptions[item] = int(input_row[item])
            
    # Goldford et al 2018 parameters
    for item in ['mu_f','sigma_f', 'omega_f','n_timesteps']:
        if item in input_row.keys():
            assumptions[item] = float(input_row[item])
    for item in ['n_timepoints']:
        if item in input_row.keys():
            assumptions[item] = int(input_row[item])
    
    
    # Community-simulator assumption parameters
    for item in ['sampling_c','sampling_D','regulation','response','supply']:
        if item in input_row.keys():
            assumptions[item] = str(input_row[item])
    for item in ['muc','sigc','q','c0','c1','fs','fw','sparsity','R0_food', 'tau']:
        if item in input_row.keys():
            assumptions[item] = float(input_row[item])
    assumptions['food'] = int(input_row['food'])
            
    # Communit-simulation params parameters. Implicit to the package
    for item in ['g','w','l','m','tau','r','n','sigma_max','nreg']:
        if item in input_row.keys():
            assumptions[item] = input_row[item]
    
    # Customized paramters 
    for item in ['n_pass', 't_propagation', 'dilution_factor', 'save_timepoint', 'n_timepoint']:
        if item in input_row.keys():
            assumptions[item] = input_row[item]

    return assumptions

# Make matrices. Modified to update the pandas syntax

def make_matrices(assumptions):
    """
    Construct consumer matrix and metabolic matrix.
    
    Returns:
    c = consumer matrix
    D = metabolic matrix
    """
    # Prepare variables 
    M, T, S, F, M_waste, type_names, family_names, waste_name, resource_index, family_names, consumer_index = extract_shapes(assumptions)
    
    # Consumer uptake rate matrix c
    if assumptions['sampling_c'] == 'Gaussian':
        c = pd.DataFrame(np.zeros((S,M)),columns=resource_index,index=consumer_index)
        #Add Gaussian-sampled values, biasing consumption of each family towards its preferred resource:
        for k in range(F):
            for j in range(T):
                if k==j:
                    c_mean = (assumptions['muc']/M)*(1+assumptions['q']*(M-assumptions['MA'][j])/assumptions['MA'][j])
                    c_var = (assumptions['sigc']**2/M)*(1+assumptions['q']*(M-assumptions['MA'][j])/assumptions['MA'][j])
                else:
                    c_mean = (assumptions['muc']/M)*(1-assumptions['q'])
                    c_var = (assumptions['sigc']**2/M)*(1-assumptions['q'])
                c.loc['F'+str(k), 'T'+str(j)] = c_mean + np.random.randn(assumptions['SA'][k],assumptions['MA'][j])*np.sqrt(c_var)
        if 'GEN' in c.index:
            c_mean = assumptions['muc']/M
            c_var = assumptions['sigc']**2/M
            c.loc['GEN'] = c_mean + np.random.randn(assumptions['Sgen'],M)*np.sqrt(c_var)
                    
    elif assumptions['sampling_c'] == 'Binary':
        assert assumptions['muc'] < M*assumptions['c1'], 'muc not attainable with given M and c1.'
        #Construct uniform matrix at total background consumption rate c0:
        c = pd.DataFrame(np.ones((S,M))*assumptions['c0']/M,columns=resource_index,index=consumer_index)
        #Sample binary random matrix blocks for each pair of family/resource type:
        for k in range(F):
            for j in range(T):
                if k==j:
                    p = (assumptions['muc']/(M*assumptions['c1']))*(1+assumptions['q']*(M-assumptions['MA'][j])/assumptions['MA'][j])
                else:
                    p = (assumptions['muc']/(M*assumptions['c1']))*(1-assumptions['q'])
                    
                c.loc['F'+str(k), 'T'+str(j)] = (c.loc['F'+str(k), 'T'+str(j)].values 
                                                + assumptions['c1']*BinaryRandomMatrix(assumptions['SA'][k],assumptions['MA'][j],p))
        #Sample uniform binary random matrix for generalists:
        if 'GEN' in c.index:
            p = assumptions['muc']/(M*assumptions['c1'])
            c.loc['GEN'] = c.loc['GEN'].values + assumptions['c1']*BinaryRandomMatrix(assumptions['Sgen'],M,p)

    elif assumptions['sampling_c'] == 'Gamma':
        c = pd.DataFrame(np.zeros((S,M)),columns=resource_index,index=consumer_index)
        #Add Gamma-sampled values, biasing consumption of each family towards its preferred resource
        for k in range(F):
            for j in range(T):
                if k==j:
                    c_mean = (assumptions['muc']/M)*(1+assumptions['q']*(M-assumptions['MA'][j])/assumptions['MA'][j])
                    c_var = (assumptions['sigc']**2/M)*(1+assumptions['q']*(M-assumptions['MA'][j])/assumptions['MA'][j])
                    thetac = c_var/c_mean
                    kc = c_mean**2/c_var
                    c.loc['F'+str(k), 'T'+str(j)] = np.random.gamma(kc,scale=thetac,size=(assumptions['SA'][k],assumptions['MA'][j]))
                else:
                    c_mean = (assumptions['muc']/M)*(1-assumptions['q'])
                    c_var = (assumptions['sigc']**2/M)*(1-assumptions['q'])
                    thetac = c_var/c_mean
                    kc = c_mean**2/c_var
                    c.loc['F'+str(k), 'T'+str(j)] = np.random.gamma(kc,scale=thetac,size=(assumptions['SA'][k],assumptions['MA'][j]))
        if 'GEN' in c.index:
            c_mean = assumptions['muc']/M
            c_var = assumptions['sigc']**2/M
            thetac = c_var/c_mean
            kc = c_mean**2/c_var
            c.loc['GEN'] = np.random.gamma(kc,scale=thetac,size=(assumptions['Sgen'],M))
    
    elif assumptions['sampling_c'] == 'Uniform':
        c = pd.DataFrame(np.zeros((S,M)),columns=resource_index,index=consumer_index)
        #Add uniformly sampled values, biasing consumption of each family towards its preferred resource:
        for k in range(F):
            for j in range(T):
                if k==j:
                    c_mean = (assumptions['muc']/M)*(1+assumptions['q']*(M-assumptions['MA'][j])/assumptions['MA'][j])
                else:
                    c_mean = (assumptions['muc']/M)*(1-assumptions['q'])
                c.loc['F'+str(k), 'T'+str(j)] = c_mean + (np.random.rand(assumptions['SA'][k],assumptions['MA'][j])-0.5)*assumptions['b']
        if 'GEN' in c.index:
            c_mean = assumptions['muc']/M
            c.loc['GEN'] = c_mean + (np.random.rand(assumptions['Sgen'],M)-0.5)*assumptions['b']
    
    elif assumptions['sampling_c'] == 'Dirichlet':
        c = pd.DataFrame(np.zeros((S,M)),columns=resource_index,index=consumer_index)
        # Sample update rate proportion from Dirichlet
        for k in range(F):
            theta = pd.Series(np.zeros(M), index = c.keys())
            theta_f = np.random.normal(assumptions['mu_f'], assumptions['sigma_f'])
            # Specialized resource
            theta['T' + str(k)] = theta_f
            # Other resources
            theta_a = np.random.uniform(0, 1, size = M-1)
            theta_a = (1 - theta_f) * theta_a / sum(theta_a)
            theta.loc[~theta.index.isin([('T' + str(k), 'R' + str(k))])] = theta_a
            c.loc['F'+str(k), :] = dirichlet(theta * assumptions['omega_f'], size = assumptions['SA'][k])
        # Sample total uptake capacity from Normal
        Ti = np.random.normal(1, 0.01, size = S)
        c = c.mul(Ti, axis = 0)
    
    else:
        print('Invalid distribution choice. Valid choices are kind=Gaussian, kind=Binary, kind=Gamma, kind=Uniform, kind = Dirichlet')
        return 'Error'

    # Metabolic matrix D
    if assumptions['sampling_D'] == 'Dirichlet':
        #SAMPLE METABOLIC MATRIX FROM DIRICHLET DISTRIBUTION
        DT = pd.DataFrame(np.zeros((M,M)),index=c.keys(),columns=c.keys())
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
    elif assumptions['sampling_D'] == 'Uniform': 
        DT = pd.DataFrame(np.zeros((M,M)),index=c.keys(),columns=c.keys())
        for i in range(T):
            for j in range(T):
                DT.loc['T'+str(i), 'T'+str(j)] = np.random.uniform(0, 1/M, size = (assumptions['MA'][i], assumptions['MA'][j]))
    
    elif assumptions['sampling_D'] == 'Zeros': 
        DT = pd.DataFrame(np.zeros((M,M)),index=c.keys(),columns=c.keys())
    else:
        print('Invalid distribution choice. Valid choices are kind = Dirichlet, kind=Uniform, kind = Zeros')
        return 'Error'
     
    # Leakiness. legacy
    l = pd.DataFrame(assumptions['l'] * np.ones((S,M)),columns=resource_index,index=consumer_index)
       
        
    return c, DT.T, l

def make_params(assumptions):
    """
    Makes a dictionary of parameters, using MakeMatrices for the matrices, MakeInitialState
    for the resource supply point, and setting everything else to 1, except l which is zero.
    
    Parameter values can be modified from 1 (or zero for l) by adding their name-value pairs
    to the assumptions dictionary.
    """

    c, D, l = make_matrices(assumptions)
    N0,R0 = MakeInitialState(assumptions)
    
    if not isinstance(assumptions['food'],int) or not isinstance(assumptions['R0_food'],int):
        params=[{'c':c,
                'm':1,
                'w':1,
                'D':D,
                'g':1,
                'l':l,
                'R0':R0.values[:,k],
                'tau':1,
                'r':1,
                'sigma_max':1,
                'nreg':10,
                'n':2
                } for k in range(assumptions['n_wells'])]
        for item in ['m','w','g','tau','r','sigma_max','n','nreg']:
            if item in assumptions.keys():
                for k in range(assumptions['n_wells']):
                    params[k][item] = assumptions[item]

    else:
        params={'c':c,
                'm':1,
                'w':1,
                'D':D,
                'g':1,
                'l':l,
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

def make_resource_dynamics(assumptions):
    """
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
    #J_out = lambda R,params: (params['l']*J_in(R,params)).dot(params['D'].T)
    J_out = lambda R,params: (J_in(R,params)).dot(params['D'].T) # remove leakiness
    
    return lambda N,R,params: (h[assumptions['supply']](R,params)
                               -(J_in(R,params)/params['w']).T.dot(N)
                               +(J_out(R,params)/params['w']).T.dot(N))

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
                             *params['w']
                             *sigma[assumptions['response']](R,params))
    #J_growth = lambda R,params: (1-params['l'])*J_in(R,params)
    J_growth = lambda R,params: J_in(R,params) # Remove leakiness
    
    return lambda N,R,params: params['g']*N*(np.sum(J_growth(R,params), axis = 1) - params['m'])

# Make initial state N0 and R0 
def make_initial_state(input_row, assumptions):
    M, T, S, F, M_waste, type_names, family_names, waste_name, resource_index, family_names, consumer_index = extract_shapes(assumptions)
    
    # Default
    N0,R0 = MakeInitialState(assumptions)
    
    # Read both N0 and R0 when provided
    N0 = pd.read_csv(input_row['output_dir'] + input_row['init_N0'])
    N0 = N0.loc[:,~N0.columns.isin(['Family', 'Species'])].astype(float)
    N0 = N0.set_index(consumer_index)
    
    R0 = pd.read_csv(input_row['output_dir'] + input_row['init_R0'])
    R0 = R0.loc[:,~R0.columns.isin(['Class', 'Resource'])].astype(float)
    R0 = R0.set_index(resource_index)
    assert N0.shape[1] == R0.shape[1], "Well numbers in the input N0 and R0 are different."
    
    return N0, R0

# Write matrices
def write_matrices(input_row, D, c, l):
    D.to_csv(input_row['output_dir'] + '00-D.csv')
    c.to_csv(input_row['output_dir'] + '00-c.csv')
    l.to_csv(input_row['output_dir'] + '00-l.csv')

# Run the simulation using the customized parameter sets
def run_simulations(input_row):
    # input_csv = 'simulation/01-input_parameters.csv' # Input file name
    # input_csv = 'simulation/02b-input_communities.csv' # Input file name
    # row_number = 0 # Which row of experiment to run
    # input_row = pd.read_csv(input_csv).loc[row_number]

    seed = int(input_row['seed'])
    exp_id = int(input_row['exp_id'])

    # Update model parameters
    assumptions = load_assumptions(input_row)

    # Update params 
    np.random.seed(seed)
    params = make_params(assumptions)[0]
    write_matrices(input_row, params['D'], params['c'], params['l'])

    # Generate equations
    def dNdt(N,R,params):
        return make_consumer_dynamics(assumptions)(N,R,params) 
    def dRdt(N,R,params):
        return make_resource_dynamics(assumptions)(N,R,params) 
    dynamics = [dNdt,dRdt]

    # Set initial state by reading in plate conditions
    N0, R0 = make_initial_state(input_row, assumptions)
    init_state = [N0,R0]
    
    # Make plate object
    Plate = Community(init_state, dynamics, params, parallel = False)
    
    # Batch culture
    if assumptions['supply'] == 'off': 
        
        for i in range(int(assumptions['n_timepoints'])):
            Plate.Propagate(T = assumptions['n_timesteps_batch'], compress_resources = False, compress_species = True)
            Plate.N.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "T" + str(i+1) + ".csv", input_row["init_N0"]))
            Plate.R.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "T" + str(i+1) + ".csv", input_row["init_R0"]))
            print("T" + str(i+1))
            if i == (assumptions['n_timepoints']-1): # last transfer
                Plate.N.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "end.csv", input_row["init_N0"]))
                Plate.R.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "end.csv", input_row["init_R0"]))
    
    # Chemostat
    elif assumptions['supply'] == 'external':
        if assumptions['save_timepoint'] == True:
            for i in range(int(assumptions['n_timepoints'])):
                Plate.Propagate(T = assumptions['n_timesteps'] / assumptions['n_timepoints'], compress_resources = False, compress_species = True)
                Plate.N.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "T" + str(i+1) + ".csv", input_row["init_N0"]))
                Plate.R.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "T" + str(i+1) + ".csv", input_row["init_R0"]))
                print("T" + str(i+1))
                if i == (assumptions['n_timepoints']-1): # last transfer
                    Plate.N.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "end.csv", input_row["init_N0"]))
                    Plate.R.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "end.csv", input_row["init_R0"]))
    
        elif assumptions['save_timepoint'] == False:
            Plate.Propagate(T = assumptions['n_timesteps'], compress_resources = False, compress_species = True)
            Plate.N.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "end.csv", input_row["init_N0"]))
            Plate.R.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "end.csv", input_row["init_R0"]))
















