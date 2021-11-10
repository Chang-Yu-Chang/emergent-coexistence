"""
Simulate reconstituted pairs from top-down assembled communities
"""

from community_selection.usertools import *
assumptions = make_assumptions("../data/mapping_files/input_internal.csv", 0)


#
def compete_reconstituted_pairs(seed = 1, community = "W0"):
    # Number of pairs
    df = pd.read_csv("../data/raw/simulation/reconstituted_pairs/community-top_down-" + str(seed) + "-" + community + "-reconstituted_pairs.txt")
    n_pairs = len(np.unique(df[df["Type"] == "consumer"].Well))


    # Update the assumptions
    assumptions.update({
        "n_wells": n_pairs,
        "output_dir": "../data/raw/simulation/",
        "overwrite_plate": "../data/raw/simulation/reconstituted_pairs/community-top_down-" + str(seed) + "-" + community + "-reconstituted_pairs.txt",
        "exp_id": "community-top_down-" + str(seed) + "-" + community + "-reconstituted_pairs",
        "seed": seed,
        "synthetic_community": True,
        "init_richness": 2,
        "l": 0.5,
        "q": 0.8,
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
    })


    params, params_simulation , params_algorithm, plate = prepare_experiment(assumptions)
    simulate_community(params = params, params_simulation = params_simulation, params_algorithm = params_algorithm,plate = plate) 


seed = 1
# list of community
df = pd.read_csv("../data/raw/simulation/community-top_down-" + str(seed) + "_composition.txt")
communities = np.unique(df[(df["Transfer"] == 20) & (df["Type"] == "consumer")].Well)

for i in range(len(communities)):
    compete_reconstituted_pairs(seed = seed, community = communities[i])
    
