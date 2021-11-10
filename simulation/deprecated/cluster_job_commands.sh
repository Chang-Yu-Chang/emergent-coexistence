#!/usr/bin/env bash

# Upload files to cluster
upload_to_cluster() {
   for i in "${@:1:$# - 1}"
   do
       scp -r $i cc2553@transfer-grace.hpc.yale.edu:${@: -1}
   done
#    scp -r $1 cc2553@transfer-grace.hpc.yale.edu:$2
}
#upload_to_cluster ~/Desktop/test.txt project/community_selection_project/
#scp -r $i cc2553@transfer-grace.hpc.yale.edu:${@: -1}


# Download files from cluster
download_from_cluster() {
    scp cc2553@transfer-grace.hpc.yale.edu:$1 $2
}
#download_from_cluster project/community_selection_project/test.txt ~/Desktop/Lab/


# Create a wrapper folder; one wrapper folder one input csv
create_wrapper_folder() {
    mkdir ~/project/invasion-network/wrapper/$1
    echo "Created ~/project/invasion-network/wrapper/$1"
    mkdir ~/project/invasion-network/data/$1
    echo "Created ~/project/invasion-network/data/$1"
    mkdir ~/project/invasion-network/job/$1
    echo "Created ~/project/invasion-network/job/$1"
}

create_wrapper_folder2(){
    create_wrapper_folder set_$1
    create_wrapper_folder synthetic_$1
}

# Remove wrapper folder
remove_wrapper_folder() {
    rm -r ~/project/invasion-network/wrapper/$1
    echo "Deleted ~/project/invasion-network/wrapper/$1"
    rm -r ~/project/invasion-network/data/$1
    echo "Deleted ~/project/invasion-network/data/$1"
    rm -r ~/project/invasion-network/job/$1
    echo "Deleted ~/project/invasion-network/job/$1"
}

remove_wrapper_folder2() {
    remove_wrapper_folder set_$1
    remove_wrapper_folder synthetic_$1
}


# Update/refresh the prerequisite for python environment
refresh_ecoprospector() {
    conda activate invnet
    pip install -e ~/project/invasion-network/community-simulator/
    pip install -e ~/project/invasion-network/ecoprospector/
    conda deactivate
}


# Make joblist
make_joblist(){
    python ~/project/invasion-network/wrapper/make_joblist.py $1 $2 $3
    #python ~/Desktop/Lab/invasion-network/Simulations_Scripts/make_joblist.py $1 $2 $3
}
#make_joblist input_independent.csv

# Make sbatch file from joblist
make_sbatch_file() {
    module load dSQ;
    job_name=$(echo $1 | awk '{ gsub("joblist_", "") ; print $0 }' | awk '{ gsub(".txt", "") ; print $0 }')
    dsq --job-file "joblist_$job_name.txt" --mem-per-cpu "$2g" -t "$3:00:00" --mail-type ALL --job-name $job_name --batch-file  "batch_$job_name.sh" --output "/home/cc2553/project/invasion-network/job/$job_name/%A_%3a.txt"
}
#make_sbatch_file joblist_independent.txt 20 10


# Make synthetic communities
make_synthetic_communities() {
    #Rscript ~/project/invasion-network/wrapper/make_synthetic_communities.R $1
    Rscript ~/project/invasion-network/wrapper/make_culturable_isolates.R $1
    Rscript ~/project/invasion-network/wrapper/make_pair_from_random_species.R $1
    Rscript ~/project/invasion-network/wrapper/make_pair_from_top_down_community.R $1
}
make_culturable_isolates() {
    Rscript ~/project/invasion-network/wrapper/make_culturable_isolates.R $1 $2
}
make_pair_from_random_species() {
    Rscript ~/project/invasion-network/wrapper/make_pair_from_random_species.R $1 $2 $3
}
make_pair_from_top_down_community() {
    Rscript ~/project/invasion-network/wrapper/make_pair_from_top_down_community.R $1 $2 $3
}


make_synthetic_communities2() {
    conda activate invnet
    cd ~/project/invasion-network/wrapper/synthetic_$1
    make_synthetic_communities input_synthetic_$1.csv
    conda deactivate
}


# Wrapper files for creating files
submit_job() {
    cd ~/project/invasion-network/wrapper/$1_$2
    make_joblist $1 input_$1_$2.csv joblist_$1_$2.txt
    make_sbatch_file joblist_$1_$2.txt 20 24
    sbatch batch_$1_$2.sh
}


# Wrapper function for cluster command
print_cluster_commands() {
x="
# Cluster
create_wrapper_folder independent_$1
create_wrapper_folder iteration_$1
create_wrapper_folder robustness_$1

# Locally
Rscript make_mapping_files.R
upload_to_cluster ~/Desktop/Lab/invasion-network/Data/Mapping_Files/input_independent_$1.csv project/invasion-network/wrapper/independent_$1/
upload_to_cluster ~/Desktop/Lab/invasion-network/Data/Mapping_Files/input_iteration_$1.csv project/invasion-network/wrapper/iteration_$1/
upload_to_cluster ~/Desktop/Lab/invasion-network/Data/Mapping_Files/input_robustness_$1.csv project/invasion-network/wrapper/robustness_$1/

# Cluster independent
cd ~/project/invasion-network/wrapper/independent_$1
make_joblist independent input_independent_$1.csv joblist_independent_$1.txt
make_sbatch_file joblist_independent_$1.txt 20 24
sbatch batch_independent_$1.sh

# Cluster iteration
cd ~/project/invasion-network/wrapper/iteration_$1
make_joblist iteration input_iteration_$1.csv joblist_iteration_$1.txt
make_sbatch_file joblist_iteration_$1.txt 20 24
sbatch batch_iteration_$1.sh

# Make synthetic community
conda activate invnet
cd  ~/project/invasion-network/wrapper/iteration_$1
make_synthetic_communities input_iteration_$1.csv
conda deactivate

# Robustness
cd ~/project/invasion-network/wrapper/robustness_$1
make_joblist independent input_robustness_$1.csv joblist_robustness_$1.txt
make_sbatch_file joblist_robustness_$1.txt 10 10
sbatch batch_robustness_$1.sh

# Download
scp 'cc2553@transfer-grace.hpc.yale.edu:~/project/invasion-network/data/independent_$1/*' ~/Dropbox/invasion-network/Data/Raw_Rebuttal/
scp 'cc2553@transfer-grace.hpc.yale.edu:~/project/invasion-network/data/iteration_$1/*_function.txt' ~/Dropbox/invasion-network/Data/Raw_Rebuttal/
scp 'cc2553@transfer-grace.hpc.yale.edu:~/project/invasion-network/data/robustness_$1/*' ~/Dropbox/invasion-network/Data/Raw_Rebuttal/

"
echo $x
}

print_cluster_commands2() {
x="
create_wrapper_folder2 $1

upload_to_cluster ~/Desktop/Lab/invasion-network/Data/Mapping_Files/input_independent_$1.csv project/invasion-network/wrapper/independent_$1/
upload_to_cluster ~/Desktop/Lab/invasion-network/Data/Mapping_Files/input_iteration_$1.csv project/invasion-network/wrapper/iteration_$1/
upload_to_cluster ~/Desktop/Lab/invasion-network/Data/Mapping_Files/input_robustness_$1.csv project/invasion-network/wrapper/robustness_$1/

submit_job independent $1
submit_job iteration $1
make_synthetic_communities2 $1
submit_job robustness $1

scp 'cc2553@transfer-grace.hpc.yale.edu:~/project/invasion-network/data/independent_$1/*' ~/Dropbox/invasion-network/Data/Raw_Rebuttal/
scp 'cc2553@transfer-grace.hpc.yale.edu:~/project/invasion-network/data/iteration_$1/*_function.txt' ~/Dropbox/invasion-network/Data/Raw_Rebuttal/
scp 'cc2553@transfer-grace.hpc.yale.edu:~/project/invasion-network/data/robustness_$1/*' ~/Dropbox/invasion-network/Data/Raw_Rebuttal/
"
echo $x
}

# Count lines in a file
cl() {
    cat $1 | wc -l
}

# Check mapping file selected function
check_mapping_files() {
sed -n 2p ~/project/invasion-network/wrapper/$1_$2/input_$1_$2.csv | awk -F "," '{ print $1 }'
}

# check data
lldata() {
    ll ~/project/invasion-network/data/$1
}

# Check job
lljob() {
    ll ~/project/invasion-network/job/$1
}

# Print job
catjob() {
    cat ~/project/invasion-network/job/$1
}



