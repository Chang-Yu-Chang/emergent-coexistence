"""
Simulate pairs of random isolates from the pool
"""

from community_selection.usertools import *

# The updated assumptions has to be relocated in a csv file
assumptions = make_assumptions("../data/mapping_files/input_internal.csv", 0)
assumptions.update({
    "n_wells": 960,
#        "SA": 50*np.ones(3),
#        "MA": 20*np.ones(3),
    "SA": 50*np.ones(3),
    "MA": np.concatenate((1*np.ones(1), 20*np.ones(2))),
    "S": 1,
    "Sgen": 0,
    "muc": 5,
    "l": 0.5,
    "q": 0.8,
    "n_transfer": 10,
    "dilution": 0.01,
    "response": "type I",
    "rich_medium": False,
    "save_composition": True,
    "composition_lograte": 1,
    "save_function": False,
    "output_dir": "../data/raw/simulation/",
    "exp_id": "pair-random_pairs-1",
    "seed": 1,
    "synthetic_community": True,
    "init_richness": 2
})


def compete_synthetics(seed = 1, init_richness = 2, l = 0, q = 0.8, assumptions = assumptions):
    """
    Randomly draw a systhetic community with given number of members
    """
    if init_richness == 1:
        filename_prefix = "monoculture"
    if init_richness == 2:
        filename_prefix = "pair"
    elif init_richness == 3:
        filename_prefix = "trio"
    elif init_richness == 4:
        filename_prefix = "quartet"
    elif init_richness == 5:
        filename_prefix = "quintet"
    elif init_richness > 5:
        filename_prefix = "diverse"
        
    # Update the assumptions
    assumptions.update({
        "output_dir": "../data/raw/simulation/",
        "exp_id": filename_prefix + "-random_isolates-leakage" + str(int(l*100)) + "-specialist" + str(int(q*100)) + "-" +  str(seed),
        "seed": seed,
        "synthetic_community": True,
        "init_richness": init_richness,
        "l": l,
        "q": q
    })
    
    if init_richness == 1:
        assumptions.update({"monoculture":True})
    else:
        assumptions.update({"monoculture":False})
    
    params, params_simulation , params_algorithm, plate = prepare_experiment(assumptions)
    plate.N = sample_from_pool2(plate.N, assumptions, init_richness = init_richness)
    simulate_community(params = params, params_simulation = params_simulation, params_algorithm = params_algorithm,plate = plate)


#
leakages = np.linspace(0, 0.9, 10)
specialist = [0, 0.3, 0.8]
for i in range(len(leakages)):
    for j in range(len(specialist)):
        compete_synthetics(seed = 1, init_richness = 2, l = leakages[i], q = specialist[j])
