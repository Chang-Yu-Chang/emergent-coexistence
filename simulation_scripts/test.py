#!/usr/local/bin/python3

import sys
import os
from community_selection.usertools import *
import timeit
import time
import inspect


#input_csv = "~/Desktop/Lab/community-selection/Data/Mapping_Files/input_independent_f1d_additive_medium2.csv" 
#input_csv = "~/Desktop/Lab/community-selection/Data/Mapping_Files/input_independent_f5_invader_suppression.csv"
#input_csv = "~/Desktop/Lab/community-selection/Data/Mapping_Files/input_iteration_f5_invader_suppression.csv" 
input_csv = "/Users/cychang/Desktop/Lab/invasion-network/data/raw/simulation/mapping_files/input_set_simple_medium.csv" 
row_number = 191

assumptions = make_assumptions(input_csv, row_number)
{i: assumptions[i] for i in ("l", "dilution", "n_propagation", "response", "sampling", "invader_strength")}
assumptions.update({
    "rn": 3,
    "rf": 3,
    "SA": 4 * np.ones(3),
    "MA": 3 * np.ones(3),
    "muc": 1,
#    "c1": 1,
    "q": 0.8, 
#    "sf": 1,
#   # "l": 0.5,
    "sampling": "Binary"
    })
assumptions["MA"]
assumptions["SA"]
np.random.seed(assumptions['seed']) 

#plate = make_plate(assumptions,params)
#plate.N.iloc[assumptions["invader_index"],:] = 0
#plate = add_community_function(plate, assumptions, params)
params, params_simulation , params_algorithm, plate = prepare_experiment(assumptions)
params["c"].iloc[1,:]

#plate.N.iloc[assumptions["invader_index"],:]

i=0
# community_function = globals()[params_algorithm["community_phenotype"][0]](plate, params_simulation = params_simulation) # Community phenotype
# print(community_function)
phenotype_algorithm = params_algorithm["community_phenotype"][i]
selection_algorithm = params_algorithm["selection_algorithm"][i]
plate.Propagate(params_simulation["n_propagation"])
community_function = globals()[params_algorithm["community_phenotype"][0]](plate, params_simulation = params_simulation) # Community phenotype
print(community_function)

richness = np.sum(plate.N >= 1/params_simulation["scale"], axis = 0) # Richness
biomass = list(np.sum(plate.N, axis = 0)) # Biomass
temp = pd.DataFrame([biomass, richness,community_function])
temp
#print(temp)
#print([biomass[i] for i in range(len(biomass))])
#print([richness[i] for i in range(len(richness))])

#Store prior state before passaging (For coalescence)
setattr(plate, "prior_N", plate.N)
setattr(plate, "prior_R", plate.R)
setattr(plate, "prior_R0", plate.R0)

# Passage and transfer matrix
transfer_matrix = globals()[selection_algorithm](community_function)
if params_simulation['monoculture']:
    plate = passage_monoculture(plate, params_simulation["dilution"])
else:
    plate.Passage(transfer_matrix * params_simulation["dilution"])



















simulate_community(params = params, params_simulation = params_simulation, params_algorithm = params_algorithm, plate = plate)

print(str(round(timeit.default_timer() - start_time, 5)) + " seconds")

#start_time = timeit.default_timer()
#{i: assumptions[i] for i in ("l", "dilution", "n_propagation", "response")}
#print(str(round(timeit.default_timer() - start_time, 5)) + " seconds")

if "invader" in assumptions["selected_function"]:
    plate.N.iloc[assumptions["invader_index"],:] = 0
