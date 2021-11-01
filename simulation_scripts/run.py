import pandas as pd
import numpy as np
from community_simulator import *
from community_simulator.usertools import *
from community_simulator.visualization import *
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import gc
from scipy import stats
import sys

def sample_d_matrix(assumptions,dtype):
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
    
    #SAMPLE METABOLIC MATRIX FROM DIRICHLET DISTRIBUTION
    DT = pd.DataFrame(np.zeros((M,M)),index=resource_index,columns=resource_index)
    for type_name in type_names:
        MA = len(DT.loc[type_name])
        if dtype == 1:
            #No metabolic structure
            p = pd.Series(np.ones(M)/(M),index = DT.keys())
            DT.loc[type_name,p.index[p!=0]] = dirichlet(p[p!=0]/assumptions['sparsity'],size=len(DT.loc[type_name]))
        if dtype == 2:
			#Self Metabolic Structure
            p = pd.Series((np.ones(M)*0.1)/M,index = DT.keys()) #Background
            p2 =  pd.Series((np.ones(M)*0.9)/MA,index = DT.keys()) #Self
            p.loc[type_name] = (p2.loc[type_name] + p.loc[type_name]).values
            DT.loc[type_name,p.index[p!=0]] = dirichlet(p[p!=0]/assumptions['sparsity'],size=len(DT.loc[type_name]))
        if dtype ==3:
			#Waste Metabolic Structure
            p = pd.Series((np.ones(M)*0.1)/M,index = DT.keys()) #Background
            p2 =  pd.Series((np.ones(M)*0.9)/M_waste,index = DT.keys()) #Self
            p.loc[waste_name] = (p2.loc[waste_name] + p.loc[waste_name]).values
            DT.loc[type_name,p.index[p!=0]] = dirichlet(p[p!=0]/assumptions['sparsity'],size=len(DT.loc[type_name]))
    return (DT.T)


def sample_c_matrix(assumptions):
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
    assert assumptions['muc'] < M*assumptions['c1'], 'muc not attainable with given M and c1.'
    c = pd.DataFrame(np.zeros((S,M)),columns=resource_index,index=consumer_index)
    #Add Gamma-sampled values, biasing consumption of each family towards its preferred resource
    for k in range(F):
        for j in range(T):
            if k==0 and j ==0:
                c_mean = (assumptions['muc']/M)*(1+assumptions['q'])
                c_var = (assumptions['sigc']**2/M)*(1+assumptions['q'])
                thetac = c_var/c_mean
                kc = c_mean**2/c_var
                c.loc['F'+str(k)]['T'+str(j)] = np.random.gamma(kc,scale=thetac,size=(assumptions['SA'][k],assumptions['MA'][j]))
            elif k==1 and j== 1:
                c_mean = (assumptions['muc']/M)*(1+assumptions['q2'])
                c_var = (assumptions['sigc']**2/M)*(1+assumptions['q2'])
                thetac = c_var/c_mean
                kc = c_mean**2/c_var
                c.loc['F'+str(k)]['T'+str(j)] = np.random.gamma(kc,scale=thetac,size=(assumptions['SA'][k],assumptions['MA'][j]))
            elif k==1 and j ==0:
                c_mean = (assumptions['muc']/M)*(1-assumptions['q'])
                c_var = (assumptions['sigc']**2/M)*(1-assumptions['q'])
                thetac = c_var/c_mean
                kc = c_mean**2/c_var
                c.loc['F'+str(k)]['T'+str(j)] = np.random.gamma(kc,scale=thetac,size=(assumptions['SA'][k],assumptions['MA'][j]))
            elif k==0 and j ==1:
                c_mean = (assumptions['muc']/M)*(1-assumptions['q2'])
                c_var = (assumptions['sigc']**2/M)*(1-assumptions['q2'])
                thetac = c_var/c_mean
                kc = c_mean**2/c_var
                c.loc['F'+str(k)]['T'+str(j)] = np.random.gamma(kc,scale=thetac,size=(assumptions['SA'][k],assumptions['MA'][j]))
    if 'GEN' in c.index:
        c_mean = assumptions['muc']/M
        c_var = assumptions['sigc']**2/M
        thetac = c_var/c_mean
        kc = c_mean**2/c_var
        c.loc['GEN'] = np.random.gamma(kc,scale=thetac,size=(assumptions['Sgen'],M))
    return(c)

def run(input_row):
	### Load Paramaters from File
	q = input_row['q']
	q2 = input_row['q2']
	seed = input_row['seed']
	exp_id = input_row['exp_id']
	S = input_row['S']
	np.random.seed(seed)
	
	#Set paramaters in model note some of these paramaters are irrelevant to simulation (c1,c0 etc).
	rep = 1 #Number of replicates per run
	a= {'sampling':'Gamma', #Sampling method
		'SA': np.ones(2)*100, #Number of species in each family
		'MA': np.ones(2)*10, #Number of resources of each type
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
		'n_wells':7*rep, #Number of independent wells
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

	#Set Random Seed
	np.random.seed(seed)	
	#Sample D matrix
	D = sample_d_matrix(a,1)
	#Sample C matrix
	c = sample_c_matrix(a)
	#Set initial state
	N0,R0 = MakeInitialState(a)
	#Construct resource environments and update initial state
	R0 = np.zeros(np.shape(R0))
	for i in range(0,rep):
		R0[0+i,(0+i*7)] = a['R0_food'] #Sugar A
		R0[int(a['MA'][:1].sum()/2)+i,(1+i*7)] = a['R0_food'] #Sugar B
		R0[a['MA'][:1].sum()+i,(2+ i*7)] = a['R0_food'] #Acid A
		R0[a['MA'][:1].sum()+i +int(a['MA'][1:].sum()/2),3+ (i*7)] = a['R0_food'] #Acid B
		R0[:,4 + (i*7)] = 0.5*(R0[:,0 + (i*7)]+R0[:,2 + (i*7)]) #Sugar + Acid
		R0[:,5 + (i*7)] = 0.5*(R0[:,0 +  (i*7)]+R0[:,1 + (i*7)]) #2 Sugars
		R0[:,6 + (i*7)] = 0.5*(R0[:,2 +  (i*7)]+R0[:,3 +(i*7)]) #2 Acids
	R0 = pd.DataFrame(R0,index=D.index,columns=N0.keys())	
	init_state=[N0,R0]
	#Create list of paramaters (differet R0 values so each well has to have different paramaters
	params=[{**a,**{'R0':R0.values[:,k],
			}} for k in range(len(N0.T))]
	#Set c and D matrix to be same for 7 environments (same species)
	for k in range(len(params)):
		params[k]['c'] = c
		params[k]['D'] = D
	#Load plate object
	Plate = Community(init_state,dynamics,params)
	#Run to steady state. Note tol 1e-3 is the same as in Marsland et al. 
	Plate.SteadyState(plot=False,tol=1e-3,verbose=False)
	#Record Error. 

	err_list = pd.Series(dtype='float64')
	N = Plate.N.copy()
	R = Plate.R.copy()
	params = Plate.params[0]
	
	#Record Error. Code taken from Bobby's jupyter notebook
	err_list = pd.Series(dtype = 'float64')
	thresh = 1e-4
	for well in N.keys():
		params = Plate.params[np.where(N.keys()==well)[0][0]]
		params['R0'] = R0[well]
		Nwell = N[well].values
		Rwell = R[well].values
		Nwell[Nwell<thresh] = 0
		err = dNdt(N[well].values,R[well].values,params)
		err = err[err != 0]
		err_list.loc[well] = np.max(np.abs(err/Nwell[Nwell>thresh]))
		
	#Save Plate data and error file.	
	Plate.N.to_csv('../Data/Raw/' + str(exp_id) + '.csv')
	err_list.to_csv('../Data/Raw/Error_' + str(exp_id) + '.csv')
	return 0

if __name__ == "__main__":
	input_csv = str(sys.argv[1]) # Input file name
	row_number = int(sys.argv[2]) # Which row of experiment to run
	input_row = pd.read_csv(input_csv).loc[row_number]
	run(input_row)
