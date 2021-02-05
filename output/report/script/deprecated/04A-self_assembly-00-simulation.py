"""
This script simulates the self-assembly of microbial communities:
"""

# Import modules
import datetime
import time
from tqdm import tqdm
from IPython.display import Image
from community_simulator import *
from community_simulator.usertools import *
from community_simulator.visualization import *
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends import backend_pdf as bpdf

# Call the essenteial functions for self-assembly
exec(open("04A-self_assembly-01-experiment_functions.py").read()) 


# Parameters

# Assumptions
assumptions = a_default.copy() # Start with default parameters
assumptions.update({'n_wells':96, 'c1' :0.01, 'muc':0.1, 'm':0, 'response':'type II'})

# Prepare experiment
params, species_pool = prepare_experiment(assumptions)

## Simulation parameters
params_simulation = {
    "n_propagation": 12, # Length of propagation, or hours within a growth cycle
    "n_transfer": 10, # Number of transfer, or number of passage
    "dilution": 1/1000, # Dilution factor for transfer
}


# 96 self-assembled communities 
print("Self-assembly")
# Make initial state
init_state = MakeInitialState(assumptions)
# Make plate
plate = Community(init_state, dynamics, params, scale = 10**6, parallel = True) 
# Simulation
simulate_community(plate, assumptions, params_simulation, file_name = "../data/self_assembly_TypeII-community", write_composition = True)

# Monocultures
print("monocultures")
# Make initial state
N0 = make_synthetic_mono(assumptions)
init_state = make_initial_state(assumptions, N0)
# Make plate
plate = Community(init_state, dynamics, params, scale = 10**6, parallel = True) 
# Simulation
simulate_community(plate, assumptions, params_simulation, file_name = "../data/self_assembly_TypeII-mono", write_composition = True)


















