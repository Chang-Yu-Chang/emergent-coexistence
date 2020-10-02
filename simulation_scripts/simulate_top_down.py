"""
Simulate top-down assembly
"""

from community_selection.usertools import *


# The updated assumptions has to be relocated in a csv file
assumptions = make_assumptions("../data/mapping_files/input_internal.csv", 0)
assumptions.update({
        "l": 0.5,
        "q": 0.8
    })


def assemble_top_down(seed = 1, assumptions = assumptions):
    
    assumptions.update({
        "n_wells": 96,
#        "SA": 50*np.ones(3),
#        "MA": 20*np.ones(3),
        "SA": 50*np.ones(3),
        "MA": np.concatenate((1*np.ones(1), 20*np.ones(2))),
        "S":1,
        "Sgen": 0,
        "muc": 5,
        "n_transfer": 20,
        "dilution": 0.01,
        "response": "type I",
        "rich_medium": False,
        "save_composition": True,
        "composition_lograte": 1,
        "save_function": False,
        "output_dir": "../data/raw/simulation/",
        "exp_id": "community-top_down-" + str(seed),
        "seed": seed,
        "synthetic_community": False,
        "init_richness": 2
    })
      
    params, params_simulation , params_algorithm, plate = prepare_experiment(assumptions)
    simulate_community(params = params, params_simulation = params_simulation, params_algorithm = params_algorithm,plate = plate)


for i in range(1,2):
    assemble_top_down(seed = i, assumptions = assumptions)

    
