#!/usr/bin/python
from community_selection import *
from community_selection.usertools import *

# Make assumptions
assumptions = make_assumptions("../data/mapping_files/input_internal.csv", 0)
params, params_simulation , params_algorithm, plate = prepare_experiment(assumptions)

#
n_time_points = 10

#
list_communtiy_function_internal = list()
list_plate_data_internal = list()
community_function = globals()[params_algorithm["community_phenotype"][0]](plate, params_simulation = params_simulation) 
 

# Save the inocula composition
if params_simulation['save_composition']:
    plate_data = reshape_plate_data(plate, params_simulation,transfer_loop_index=0)  # Initial state
    plate_data["Time"] = 0
    list_plate_data_internal.append(plate_data)

# Save the initial community function + richness + biomass
if params_simulation['save_function']:
    richness = np.sum(plate.N >= 1/params_simulation["scale"], axis = 0) # Richness
    biomass = list(np.sum(plate.N, axis = 0)) # Biomass
    function_data = reshape_function_data(params_simulation,community_function, richness, biomass, transfer_loop_index =0)
    function_data["Time"] = 0
    list_communtiy_function_internal.append(function_data)


# Propagation
print("\nStart propogation")
# Run simulation
for i in range(0, params_simulation["n_transfer"]):
    # Algorithms used in this transfer
    phenotype_algorithm = params_algorithm["community_phenotype"][i]
    selection_algorithm = params_algorithm["selection_algorithm"][i]
    print("Transfer " + str(i+1))
    

    for j in range(n_time_points):
        # Propogation
        plate.Propagate(params_simulation["n_propagation"]/n_time_points)    

        # Function
        community_function = globals()[params_algorithm["community_phenotype"][0]](plate, params_simulation = params_simulation) # Community phenotype

        # Plate composition
        plate_data = reshape_plate_data(plate, params_simulation,transfer_loop_index=i+1)  # Initial state
        plate_data["Time"] = j
        list_plate_data_internal.append(plate_data)

        # Community function 
        richness = np.sum(plate.N >= 1/params_simulation["scale"], axis = 0) # Richness
        biomass = list(np.sum(plate.N, axis = 0)) # Biomass
        function_data = reshape_function_data(params_simulation,community_function, richness, biomass, transfer_loop_index = i+1)
        function_data["Time"] = j
        list_communtiy_function_internal.append(function_data)
        
    # Transfer
    transfer_matrix = globals()[selection_algorithm](community_function)
    plate.Passage(transfer_matrix * params_simulation["dilution"])



# Output 
composition_filename = "../Data/Raw/" + params_simulation['exp_id'] + '_composition.csv'
function_filename = "../Data/Raw/" + params_simulation['exp_id'] + '_function.csv'   
pd.concat(list_plate_data_internal).to_csv(composition_filename ,index=False)
pd.concat(list_communtiy_function_internal).to_csv(function_filename ,index=False)



