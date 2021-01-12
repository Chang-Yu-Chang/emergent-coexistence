#!/usr/bin/python
"""
Read the input csv and output a shell script or txt that contains jobs 
arg1 = {independent, iteration, robustness}
arg2 = input csv file
arg3 = output file name
"""

import os
import sys
import subprocess
import csv

simulation_job_type = str(sys.argv[1]) # Input csv file type
input_csv = str(sys.argv[2]) # Input file 
output_joblist = str(sys.argv[3]) # Output file name
iteration_rounds = range(1, 24)
#iteration_rounds = range(11, 24)
cluster_envir_string = "source activate invnet; "
commandline_tool = "ecoprospector "
#list_protocols = ["simple_screening"] + ["iteration_"+str(i) for i in range(1,8)]
#list_protocols = ["simple_screening"] + ["iteration_"+str(i) for i in range(]

if simulation_job_type == "set" or simulation_job_type == "synthetic":
    # Read input csv
    with open(input_csv,"r") as f:
        reader = csv.reader(f,delimiter = ",")
        df_input = list(reader)
    df_input.pop(0)
    
    # Number of experiments = number of jobs
    with open(input_csv,"r") as f:
        reader = csv.reader(f,delimiter = ",")
        data = list(reader)
        total_experiments = len(data) - 1 # Not counting header
    
    # Total number of seeds
    n_seeds = int(df_input[len(df_input)-1][2])
    
    # Extract the list of exp_id
    list_exp_id = list()
    list_seeds = list()
    for i in range(len(df_input)): 
        exp_id = df_input[i][3] # expid 
        seed = df_input[i][2] # seed
        list_exp_id.append(exp_id)
        list_seeds.append(seed)


if simulation_job_type == "set" or simulation_job_type == "synthetic":
    fout = open(output_joblist, "wt")
    line_initial = cluster_envir_string
    
    counter = 0
    line_temp = cluster_envir_string
    for i in range(total_experiments):
        # Monoculture takes longer, run one monoculture in one job
        if "monoculture" in list_exp_id[i]:
            line_monoculture = cluster_envir_string + commandline_tool + input_csv + " " + str(i) + ";\n"
            fout.write(line_monoculture)
            continue
            
        line_temp += commandline_tool + input_csv + " " + str(i) + ";"
        counter += 1
    
        # Group 10 simulations into one job; if the number of expeirment cannot be fully divided by 10, also print the residue
        if counter == 2 or i == (total_experiments-1): 
            counter = 0
            line_temp += "\n"
            fout.write(line_temp)
            line_temp = line_initial
    
    fout.close()

elif simulation_job_type == "make_pairs":
    fout = open(output_joblist, "wt")
    line_initial = "\n" + cluster_envir_string + " source ~/project/invasion-network/wrapper/cluster_job_commands.sh; "
    for i in range(1,9):
        line_temp = line_initial + "cd ~/project/invasion-network/wrapper/set_simple_medium; " + "make_culturable_isolates input_set_simple_medium.csv " + "simple_medium" + str(i) + ";" +  "cd ~/project/invasion-network/wrapper/synthetic_simple_medium; " + "make_pair_from_random_species ../set_simple_medium/input_set_simple_medium.csv input_synthetic_simple_medium.csv simple_medium" + str(i) + ";" + "make_pair_from_top_down_community ../set_simple_medium/input_set_simple_medium.csv input_synthetic_simple_medium.csv simple_medium" +str(i)
        fout.write(line_temp)
    for i in range(1,9):
        line_temp = line_initial + "cd ~/project/invasion-network/wrapper/set_rich_medium; " + "make_culturable_isolates input_set_rich_medium.csv " + "rich_medium" + str(i) + ";" + "cd ~/project/invasion-network/wrapper/synthetic_rich_medium; " + "make_pair_from_random_species ../set_rich_medium/input_set_rich_medium.csv input_synthetic_rich_medium.csv rich_medium" + str(i) + ";" + "make_pair_from_top_down_community ../set_rich_medium/input_set_rich_medium.csv input_synthetic_rich_medium.csv rich_medium" +str(i)
        fout.write(line_temp)
    fout.close()

print("Created " + output_joblist)
