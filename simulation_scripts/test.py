#!/usr/local/bin/python3

import sys
import os
from community_selection.usertools import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

#input_csv = "~/Desktop/Lab/community-selection/Data/Mapping_Files/input_independent_f1d_additive_medium2.csv" 
#input_csv = "~/Desktop/Lab/community-selection/Data/Mapping_Files/input_independent_f5_invader_suppression.csv"
#input_csv = "~/Desktop/Lab/community-selection/Data/Mapping_Files/input_iteration_f5_invader_suppression.csv" 
input_csv = "~/Desktop/Lab/invasion-network/data/raw/simulation/mapping_files/input_set_simple_medium.csv" ; row_number = 280

assumptions = make_assumptions(input_csv, row_number)
{i: assumptions[i] for i in ("l", "q", "sampling", "SA", "sampling_D", "fww")}
assumptions.update({
#     "rn": 2,
#     "rf": 2,
     "SA": 100 * np.ones(2),
     "MA": 10 * np.ones(2),
#     "muc": 10,
#    "c1": 1,
     "q": 0.5, 
#    "sf": 1,
#   # "l": 0.5,
     "sampling": "Gamma",
#    "fs": 0.49,
#    "fw": 0.49,
#    "sparsity": 0.01,
    "sampling_D": "default"
    })
assumptions["MA"]
assumptions["SA"]
np.random.seed(assumptions['seed']) 

#plate = make_plate(assumptions,params)
#plate.N.iloc[assumptions["invader_index"],:] = 0
#plate = add_community_function(plate, assumptions, params)
params, params_simulation , params_algorithm, plate = prepare_experiment(assumptions)
def plot_matrix(X):
    fig, ax = plt.subplots()
    im = ax.imshow(X, interpolation='nearest', cmap=plt.cm.gist_gray_r)
    numrows, numcols = X.shape
    
    def format_coord(x, y):
        col = int(x + 0.5)
        row = int(y + 0.5)
        if col >= 0 and col < numcols and row >= 0 and row < numrows:
            z = X[row, col]
            return 'x=%1.4f, y=%1.4f, z=%1.4f' % (x, y, z)
        else:
            return 'x=%1.4f, y=%1.4f' % (x, y)
    ax.format_coord = format_coord
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="10%", pad=0.1)
    plt.colorbar(im, cax=cax)
    plt.show()
plot_matrix(params["D"])
plot_matrix(params["c"])
#plate.N.iloc[assumptions["invader_index"],:]





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
