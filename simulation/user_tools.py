import sys
import re
# ignore the future warning from pandas future updates
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import numpy as np
from community_simulator import *
from community_simulator.usertools import *

# input_csv = 'simulation/02b-input_communities.csv' # Input file name
# row_number = 0 # Which row of experiment to run
# input_row = pd.read_csv(input_csv).loc[row_number]

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
    waste_name = type_names[assumptions['waste_type']]
    resource_index = [[type_names[m] for m in range(T) for k in range(assumptions['MA'][m])], resource_names]
    consumer_index = [[family_names[m] for m in range(F) for k in range(assumptions['SA'][m])]+['GEN' for k in range(assumptions['Sgen'])],consumer_names]
    return M, T, S, F, type_names, family_names, waste_name, resource_index, family_names, consumer_index

# Load model parameters from the input csv
def load_assumptions(input_row):
    
    # Set parameters in model
    assumptions = a_default.copy()
    assumptions['SA'] = np.ones(int(input_row['fa'])) * int(input_row['sa']) # Number of species per specialist family
    assumptions['MA'] = np.ones(int(input_row['fa'])) * int(input_row['ma']) # Number of resource per resource class
    assumptions['Sgen'] = int(input_row['Sgen'])
    assumptions['n_wells'] = int(input_row['n_wells'])
    assumptions['n_communities'] = int(input_row['n_communities'])
    assumptions['m'] = 0 # Set m=0 to turn off maintanence cost. i.e., no cell dies
    assumptions['S'] = 2 # Number of initial per-well species sampled from the species pool
    
    # Community-simulator assumption parameters
    for item in ['sampling','regulation','response','supply']:
        if item in input_row.keys():
            assumptions[item] = str(input_row[item])
    for item in ['muc','sigc','q','c0','c1','fs','fw','sparsity','R0_food']:
        if item in input_row.keys():
            assumptions[item] = float(input_row[item])
    assumptions['food'] = int(input_row['food'])
            
    # Communit-simulation params paramters
    for item in ['m', 'w', 'g', 'l','tau','r','sigma_max','nreg','n']:
        if item in input_row.keys():
            assumptions[item] = input_row[item]
    
    # Customized paramters 
    for item in ['n_pass', 't_propagation', 'dilution_factor']:
        if item in input_row.keys():
            assumptions[item] = input_row[item]

    return assumptions

# Make matrices. It's modified to update the pandas syntax
def make_matrices(assumptions):
    """
    Construct consumer matrix and metabolic matrix.
    
    assumptions = dictionary of metaparameters
        'sampling' = {'Gaussian','Binary','Gamma'} specifies choice of sampling algorithm
        'SA' = number of species in each family
        'MA' = number of resources of each type
        'Sgen' = number of generalist species
        'muc' = mean sum of consumption rates
        'sigc' = standard deviation for Gaussian sampling of consumer matrix
        'q' = family preference strength (from 0 to 1)
        'c0' = row sum of background consumption rates for Binary sampling
        'c1' = specific consumption rate for Binary sampling
        'fs' = fraction of secretion flux into same resource type
        'fw' = fraction of secretion flux into waste resource type
        'sparsity' = effective sparsity of metabolic matrix (from 0 to 1)
        'wate_type' = index of resource type to designate as "waste"
    
    Returns:
    c = consumer matrix
    D = metabolic matrix
    """
    #PREPARE VARIABLES
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
    
    #PERFORM GAUSSIAN SAMPLING
    if assumptions['sampling'] == 'Gaussian':
        #Initialize dataframe:
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
                    
    #PERFORM BINARY SAMPLING
    elif assumptions['sampling'] == 'Binary':
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

    elif assumptions['sampling'] == 'Gamma':
        #Initialize dataframe
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
        #PERFORM GAUSSIAN SAMPLING
    elif assumptions['sampling'] == 'Uniform':
        #Initialize dataframe:
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
    
    else:
        print('Invalid distribution choice. Valid choices are kind=Gaussian, kind=Binary, kind=Gamma, kind=Uniform.')
        return 'Error'

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
        
    return c, DT.T

# Make l matrix
def make_l_matrix(assumptions):
    M = np.sum(assumptions['MA'])
    T = len(assumptions['MA'])
    S = np.sum(assumptions['SA'])+assumptions['Sgen']
    F = len(assumptions['SA'])
    resource_names = ['R'+str(k) for k in range(M)]
    type_names = ['T'+str(k) for k in range(T)]
    family_names = ['F'+str(k) for k in range(F)]
    consumer_names = ['S'+str(k) for k in range(S)]
    waste_name = type_names[assumptions['waste_type']]
    resource_index = [[type_names[m] for m in range(T) for k in range(assumptions['MA'][m])],
                      resource_names]
    consumer_index = [[family_names[m] for m in range(F) for k in range(assumptions['SA'][m])]
                      +['GEN' for k in range(assumptions['Sgen'])],consumer_names]
    l = pd.DataFrame(assumptions['l'] * np.ones((S,M)),columns=resource_index,index=consumer_index)
    
    return l

# Make initial state N0 and R0 
def make_initial_state(input_row, assumptions):
    M, T, S, F, type_names, family_names, waste_name, resource_index, family_names, consumer_index = extract_shapes(assumptions)
    
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
    
    # Load only one species' D
    # Uncomment to load all species
    # for s in range(S):
    #     D[s].to_csv(input_row['output_dir'] + 'D_S' + str(s) +  '_seed' + str(input_row['seed']) + '.csv')
    #S = np.sum(input_row['sa'] * input_row['fa'])

    D.to_csv(input_row['output_dir'] + '00-D.csv')
    c.to_csv(input_row['output_dir'] + '00-c.csv')
    l.to_csv(input_row['output_dir'] + '00-l.csv')

# Run the simulation using the customized parameter sets
def run_simulations(input_row):
    seed = int(input_row['seed'])
    exp_id = int(input_row['exp_id'])

    # Update model parameters
    assumptions = load_assumptions(input_row)
    # for i in range(len(assumptions)):
    #     print(str(list(assumptions)[i]) + " = " + str(list(assumptions.values())[i]))
    
    # Sample matrices
    np.random.seed(seed)
    c, D = make_matrices(assumptions)
    l = make_l_matrix(assumptions)
    write_matrices(input_row, D, c, l)
    
    # Update params 
    params = MakeParams(assumptions)[0]
    params['l'] = l
    params['c'] = c
    params['D'] = D
    
    # Generate equations
    def dNdt(N,R,params):
        return MakeConsumerDynamics(assumptions)(N,R,params) 
    def dRdt(N,R,params):
        return MakeResourceDynamics(assumptions)(N,R,params) 
    dynamics = [dNdt,dRdt]

    # Set initial state by reading in plate conditions
    #init_state = MakeInitialState(assumptions)
    N0, R0 = make_initial_state(input_row, assumptions)
    init_state = [N0,R0]
    
    # Make plate object
    Plate = Community(init_state, dynamics, params, parallel = False)
    print("Start passaging")
    for i in range(assumptions['n_pass']): # number of passages
        Plate.Propagate(T = assumptions['t_propagation'], compress_resources = False, compress_species = True)
        print("T" + str(i+1))
        Plate.N.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "T" + str(i+1) + ".csv", input_row["init_N0"]))
        Plate.R.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "T" + str(i+1) + ".csv", input_row["init_R0"]))
        if i == (assumptions['n_pass']-1): # last transfer
            Plate.N.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "end.csv", input_row["init_N0"]))
            Plate.R.round(2).to_csv(input_row['output_dir'] + re.sub("init.csv", "end.csv", input_row["init_R0"]))
        Plate.Passage(f = np.eye(assumptions['n_wells'])*assumptions['dilution_factor'], refresh_resource = True)






