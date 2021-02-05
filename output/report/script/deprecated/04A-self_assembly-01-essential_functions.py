#' Essential functions for simulation in self-assembly 

# Make dynanmics by default we will use the microbial consumer resource model
def dNdt(N,R,params):
    return MakeConsumerDynamics(assumptions)(N,R,params)
def dRdt(N,R,params):
    return MakeResourceDynamics(assumptions)(N,R,params)
dynamics = [dNdt,dRdt]

# Create a regional species pool
def make_regional_pool(assumptions):
    # Total number of species (specialist + generalist)
    S_tot = int(np.sum(assumptions['SA']) + assumptions['Sgen']) 
    # Assign drawn values based on power-law distribution
    pool = np.random.power(1, size  = S_tot)  # Creating an array with M rows an n_wells columns with 0 entries
    return pool/np.sum(pool)  # Relative species abundance in regional pool

# Sample communities from regional species pool
def sample_from_pool(plate_N, pool, scale=10**6, inocula=10**6):
    N0 = np.zeros((plate_N.shape))
    consumer_index = plate_N.index
    well_names = plate_N.columns
    
    for k in range(plate_N.shape[1]):
        consumer_list = np.random.choice(len(pool), size=len(pool), replace=True, p=pool)
        my_tab = pd.crosstab(index=consumer_list, columns="count") # Calculate the biomass count
        N0[my_tab.index.values,k] = np.ravel(my_tab.values / np.sum(my_tab.values) * inocula / scale) # Scale to sum

    # Make data.frame
    N0 = pd.DataFrame(N0,index=consumer_index,columns=well_names)
    return N0


def prepare_experiment(assumptions):
    
    
    
    return 



# Reshape the plate resource and consumer matrix for saving into a txt file
def reshape_plate_data(plate, transfer_loop_index):
    # Temporary function for adding variables to and melting df
    def melt_df(plate_df, data_type = "consumer"):
        # Consumers
        temp_df = pd.DataFrame(plate_df)
        total_number = temp_df.shape[0]
        
        ## Add variables
        temp_df["Type"] = np.repeat(data_type, total_number)
        temp_df["ID"] = range(total_number)
        temp_df["Transfer"] = np.repeat(str(transfer_loop_index), total_number)
        temp_df["Assembly"] = np.repeat("self-assembly", total_number)
        
        ## Melt the df
        temp_df = pd.melt(temp_df, id_vars = ["Transfer", "Assembly", "Type", "ID"], var_name = "Community", value_name = "Abundance")
        temp_df = temp_df[["Assembly", "Community", "Transfer", "Type", "ID", "Abundance"]]
        temp_df = temp_df[temp_df.Abundance != 0] # Remove zero abundances
        return temp_df
        
    # Melt the df
    temp_plate = plate.copy() # Copy the original plate 
    df_N = melt_df(temp_plate.N, data_type = "consumer")
    df_R = melt_df(temp_plate.R, data_type = "resource")
    
    # Concatenate dataframes
    merged_df = pd.concat([df_N, df_R]) 
    merged_df["Index"] = list(range(0, merged_df.shape[0]))
    merged_df.set_index("Index", inplace = True)

    return merged_df # Return concatenated dataframe
    
# Create synthetic community with all pairwise combination
def synthetic_mono(species_pool):
    """
    Make the synthetic community inocula
    
    species_pool = relative abundance of species in the pool
    
    Return:
    N0 = initial consumer populations
    """
    # Number of consumers in the species pool
    pool_richness = len(species_pool)

    # Create empty plate
    N0 = np.zeros(shape = [pool_richness, pool_richness])
    
    # Fill the plate swith 
    for i in range(pool_richness):
        N0[i, i] = 1

    # Make data.frame for community-simulator input
    N0 = pd.DataFrame(N0, columns = ["W" + str(i) for i in range(pool_richness)])
    
    return N0

# Create synthetic communities with all pairwise combination
def synthetic_pairs(consumer_list, species_pool, initial_frequency = [[0.5, 0.5], [0.05, 0.95], [0.95, 0.05]]):
    """
    Make the synthetic community inocula
    """
    # stopifnot, max(species_list) <= len(species_pool)
    from itertools import combinations
    
    # Number of consumers in the community or in the species pool
    community_richness = len(consumer_list)
    pool_richness = len(species_pool)

    # All possible pairs of consumers
    consumer_pairs = list(combinations(consumer_list, 2))
    
    # Create empty plate
    N0 = np.zeros(shape = [pool_richness, len(consumer_pairs) * len(initial_frequency)])
    
    # Fill the plate swith 
    for k in range(len(initial_frequency)):
        for i in range(len(consumer_pairs)):
            N0[consumer_pairs[i], k * len(consumer_pairs) + i] = 1 * initial_frequency[k]

    # Make data.frame for community-simulator input
    N0 = pd.DataFrame(N0, columns = ["W" + str(i) for i in range(len(consumer_pairs) * len(initial_frequency))])
    
    # Make the pair list
    pairs_list = pd.DataFrame(consumer_pairs * len(initial_frequency), columns = ["Isolate1", "Isolate2"])
    pairs_list["Community"] =  ["W" + str(i) for i in range(len(consumer_pairs) * len(initial_frequency))]
    pairs_list["Pair"] =  ["P" + str(i) for i in range(len(consumer_pairs))] * len(initial_frequency)
    pairs_list["InitialFrequency"] =  ["IF" + str(i) for i in range(len(initial_frequency))] * len(consumer_pairs)

    return N0, pairs_list



# Simulate community dynamics

params_simulation = {
    "n_propagation" = 12, # Lenght of propagation, or hours within a growth cycle
    "n_transfer" = 10, # Number of transfer, or number of passage
    "dilution" = 1/125, # Dilution factor for transfer
}




def simulate_community(plate, assumptions, params_simulation, file_name = "community", replicate_index = 1, write_composition = False):
    """
    Main function for simulating community dynamics using consumer-resource model
    
    assumptions = dictionary of metaparameters from community-simulator
    plate = plate with 
    
    """
    
    # Setup
    ## Parameters for different treatments
    np.random.seed(0) # Global random seed (i.e all participants)

    ## Initial state for the plate based on the parameters
    #assumptions = a_default.copy() # Start with default parameters
    assumptions.update({'n_wells':n, 'c1' :0.01, 'muc':0.1, 'm':0}) # Switch off mortality
    init_state = MakeInitialState(assumptions) 
    params = MakeParams(assumptions)
    species_pool = make_regional_pool(assumptions) # Generate a species pool
    
    # Make plate
    np.random.seed(0) # Set seed to get fixed initial plate
    ## Make plate by default setting. The inital community composition is updated later
    plate = Community(init_state, dynamics, params, scale = 10**6, parallel = True) 
    ## Populate the plate by sampling from species pool
    plate.N = sample_from_pool(plate.N, species_pool) 
    
    # Simulation
    ## Create empty list for recording community/resource composition and community phenotype 
    plate_data_list = list()
    
    ## Save the inocula composition
    plate_data = reshape_plate_data(plate, transfer_loop_index = 0) 
    plate_data_list.append(plate_data)

    ## Run simulation
    for i in range(0, params_simulation["n_transfer"]):
        # Propagation
        plate.Propagate(params_simulation["n_propagation"])
        
        # Save composition to an empty data
        if write_composition == True:
            # Convert commnity_function into a df and melt the df
            plate_data = reshape_plate_data(plate, transfer_loop_index = i + 1)
            plate_data_list.append(plate_data)


        # Transfer/passage by usigg transfer matrix. For simple passage, the transfer matrix is an identity matrix
        transfer_matrix = np.eye(n) * params_simulation["dilution"]
        plate.Passage(transfer_matrix)
        
        # Print the propagation progress
        print("propagation: " + str(i)) 
    
    # Output the concatenated melted data.frame
    if write_composition == True:
        plate_save = pd.concat(plate_data_list)
        plate_save.to_csv("data/self_assembly-rep" + "{:02d}".format(replicate_index) + "-" + file_name + ".txt", index = False)
        print("Species and resource composition recorded in folder data/passage/")

# Simulate pairwise dynamics
def simulate_pairs(assumptions, community, consumer_list, initial_frequency = [[0.5, 0.5], [0.05, 0.95], [0.95, 0.05]], n_propagation = 24, n_transfer = 10, dilution=1/1000, replicate_index = 1, write_composition = False):
    """
    Main function for simulating community dynamics using consumer-resource model
    """
    # Setup
    ## Parameters for different treatments
    np.random.seed(0) # Global random seed (i.e all participants)

    ## Initial state for the plate based on the parameters
    assumptions_temp = assumptions.copy()
    n = int(len(consumer_list) * (len(consumer_list) - 1) / 2 * len(initial_frequency))
    assumptions_temp.update({'n_wells': n, 'c1' :0.01, 'muc':0.1, 'm':0}) # Switch off mortality
    init_state = MakeInitialState(assumptions_temp) 
    params = MakeParams(assumptions_temp)
    species_pool = make_regional_pool(assumptions_temp) # Generate a species pool
    
    # Make plate
    np.random.seed(0) # Set seed to get fixed initial plate
    
    ## Make plate by default setting. The inital community composition is updated later
    plate = Community(init_state, dynamics, params, scale = 10**6, parallel = True) 
    
    ## Populate the plate by synthetically putting together pairs
    plate.N, pairs_list = synthetic_pairs(consumer_list = consumer_list, species_pool = species_pool, initial_frequency = initial_frequency) 
    
    # Simulation
    ## Create empty list for recording community/resource composition and community phenotype 
    plate_data_list = list()
    
    ## Save the inocula composition
    plate_data = reshape_plate_data(plate, transfer_loop_index = 0) 
    plate_data_list.append(plate_data)

    ## Run simulation
    for i in range(0, n_transfer):
        # Propagation
        plate.Propagate(n_propagation)
        
        # Save composition to an empty data
        if write_composition == True:
            # Convert commnity_function into a df and melt the df
            plate_data = reshape_plate_data(plate, transfer_loop_index = i + 1)
            plate_data_list.append(plate_data)

        # Transfer/passage by usigg transfer matrix. For simple passage, the transfer matrix is an identity matrix
        transfer_matrix = np.eye(n) * dilution
        plate.Passage(transfer_matrix)
        
        # Print the propagation progress
        print("propagation: " + str(i)) 
    
    # Output the concatenated melted data.frame
    if write_composition == True:
        plate_save = pd.concat(plate_data_list)
        plate_save = plate_save.merge(pairs_list[["Community", "Pair", "InitialFrequency" ]], on = "Community")
        plate_save = plate_save[["Assembly", "Transfer", "Pair", "InitialFrequency", "Type", "ID", "Abundance"]]

        plate_save.to_csv("data/self_assembly-rep" + "{:02d}".format(replicate_index) + "-community-" + community + "-pairs.txt", index = False)


















if False: 
    # Plot community function as a function of time    
    def plot_community_function(function_df):
        import matplotlib.pyplot as plt
        time = range(0, len(function_df))
        plt.plot(time,function_df,'ko', markersize=2)
        ax = plt.gca()
        ax.set_xlabel("transfer")
        ax.set_ylabel("Community function")
        plt.show()
    
    # Additive community function
    def additive_community_function(plate, sigma = 0.001): #Sigma is the measurement error
        # Number of communities on the plate
        N_tot = plate.N.shape[1]
        
        # Assign additive neutral traits to each of species
        ## Number of species in pools
        S_tot = plate.N.shape[0]
        ## Assign additive traits to each species. Set seed so that each species always has the same trait value
        np.random.seed(0)
        traits = np.random.normal(0, 0.1, size=S_tot)
        
        # Sum of speices traits. Sum up by columns 
        community_traits = np.sum(plate.N.values*traits[:,None],axis=0)
        
        # Impose measurement error which is drawn from a normal distribution
        measurement_error = np.random.normal(0,sigma,N_tot)
    
        return(community_traits * (1 + measurement_error))
        
    
    # Algorithms for transfer
    def pairwise_XZ(community_function):
        import itertools
        n = len(community_function)
        # Community function per transfer
        sorted_community_function = np.sort(community_function)
        # cutoff top ten 
        top_10 = sorted_community_function[86]
        winner_index = np.where(community_function >= top_10)
        pairwise = list(itertools.combinations(winner_index[0], 2))
        # Empty transfer matrix
        t = np.zeros((n, n))
        c=0
        for i in range (0,90): #populate rows one by one
            t[i,pairwise[c][0]]=1
            t[i,pairwise[c][1]]=1
            c = c+1
            if c == 45:
                c = 0
        c=0
        for i in range (90,96): 
            t[i,winner_index[0][c]]=1
            c = c+1
        return t
        
    # Compute the distances from the target resource 
    def resource_distance_community_function(plate,low=0,high=0.1,sigma = 0.01): #Sigma is the measurement error
       R_tot = plate.R.shape[0]
       well_tot = plate.R.shape[1]
       # Assign additive traits to each species
       if not hasattr(resource_distance_community_function, 'R_target'):
           np.random.seed(0)
           R_target = np.concatenate([np.zeros(1),np.random.uniform(low, high, size=R_tot-1)])
       R_dist = np.sqrt(np.sum((np.tile(R_target,(well_tot,1)) - plate.R.T)**2,axis=1))
       return np.array(R_dist.T)* -1 #(so we select for positive community function)
    
    # Migrate first 8 
    def migrate_first_half(community_function):
        # Number of wells
        n_wells = len(community_function)
        
        # Migration
        migration_factor = [1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0] * int(n_wells/16)
    
        return migration_factor
    
    # Migrate from species pool to the plate 
    def migrate_from_pool(plate, pool, migration_factor):
        # Migration plate
        migration_plate = SampleFromPool(plate.N, pool) * migration_factor
        
        # Migration
        plate_migrated = plate.N + migration_plate 
    
        return plate_migrated
    

    # Plot the transfer matrix
    def plot_transfer_matrix(t):
        t    
        fig,ax=plt.subplots()
        sns.heatmap(t,ax=ax)
        ax.set_xlabel('Old well',fontsize=14)
        ax.set_ylabel('New well',fontsize=14)
        ax.set_title(r'Transfer Matrix $f$',fontsize=14)
        plt.show()
        
    # Algorithms for transfer
    def pairwise_XZ(community_function):
        import itertools
        n = len(community_function)
        # Community function per transfer
        sorted_community_function = np.sort(community_function)
        # cutoff top ten 
        top_10 = sorted_community_function[86]
        winner_index = np.where(community_function >= top_10)
        pairwise = list(itertools.combinations(winner_index[0], 2))
        # Empty transfer matrix
        t = np.zeros((n, n))
        c=0
        for i in range (0,90): #populate rows one by one
            t[i,pairwise[c][0]]=1
            t[i,pairwise[c][1]]=1
            c = c+1
            if c == 45:
                c = 0
        c=0
        for i in range (90,96): 
            t[i,winner_index[0][c]]=1
            c = c+1
        return t
    
    def exponealing(community_function, propagation_time = 1, n_propagation=20, n_select=0.06):
        from scipy.stats import expon
        
        # Read number of wells 
        n = len(community_function)
        ngens = n_propagation #number of generations
    
        # Community function per transfer
        sorted_community_function = np.sort(community_function)
        ranked_index = np.argsort(community_function)
        #print(community_function[ranked_index])
        #print(sorted_community_function)
        
        # cutoff for selecting communities
        cut_off = sorted_community_function[int(np.round(len(community_function)*(1-n_select))) - 1]
        
        # Winner wells
        winner_index = np.where(community_function >= cut_off)
        
        # xyzzy
        if not hasattr(exponealing, 'xyzzy'): 
            exponealing.xyzzy = -1 
            
        exponealing.xyzzy = propagation_time
     #   print(exponealing.xyzzy)
        
        # Empty transfer matrix
        t = np.zeros((n, n))
        t_x = range(n)
        #t_y = np.repeat(winner_index, int(np.round(1/n_select)))
        #t_y = t_y[:n]
        #print(t_y)
        # Fill in the transfer matrix
        for i in range(len(t_x)):
            t[t_x[i], t_x[i]] = 1   
          
        # add in the random pairwise exponential tail annealing distribution 
    
        annealing_curve = np.zeros(ngens+1)
        for i in range(ngens+1):
            annealing_curve[i] = 1.5 + 28.5/((.95*i)+1) 
        
        data_expon = np.round(expon.rvs(scale=annealing_curve[exponealing.xyzzy],loc=0,size=90)) % 90
        # print(data_expon)
       
        
        for i in range(len(t_x)-6):
            t[ranked_index[t_x[i]], ranked_index[95-int(data_expon[i])]] = 1  
        if exponealing.xyzzy >= ngens-2:
            t = np.zeros((n, n))
            t_x = range(n)
             # last rond Fill in the transfer matrix with best output
            for i in range(len(t_x)):
                t[t_x[i], ranked_index[95]] = 1   
          
        # plot_transfer_matrix(t)
        return t
