
Scripts for the preprint entitled [Emergent coexistence in multispecies microbial communities](https://www.biorxiv.org/content/10.1101/2022.05.20.492860v2)


This repository includes three types of scripts for:

1. Command-line tools for image segmentation and random forest classification
2. Implementing analyses including 16S sequences, bootstrapping, and networks
3. Generating the figures for main text and supplements, as well as the supplementary pdfs


# Cloning this repository

Change the current working directory to the location where you want the cloned directory

```
$ cd /directory/to/where/you/want/
$ git clone https://github.com/Chang-Yu-Chang/emergent-coexistence
```

# Downloading the raw data to Dropbox

Link to the raw data 

https://www.dropbox.com/sh/60rw7xmerown3vk/AACMzcZMR-DzF_5vG3KF0CG-a?dl=0

The data are currently stored in the directory `~/Dropbox/lab/emergent-coexistence/data/`, with three subdirectories `raw/`, `temp/`, and `output/` where the all scripts are written to read these absolute directories. For the following scripts to work, this directory has to be set up.


# Setup to reproduce the analysis

## Step 0.1 Specifying metadata

`00-metadata.R` stores all metadata used for analysis, including the folder directory, pipeline scripts, feature names, etc.

Edit this script to specify three folders for the scripts to work:

- `folder_script` is the directory of analysis scripts
- `folder_pipeline` is the directory for saving both the raw and temporary images. For this project, a minimum space of 200 GB is required.
- `folder_data` stores the raw data from sequencing and OD, as well as the processed image data. These data are csv or Rdata formats. 

For this project on my local end, I specified the directory as followed. I will abbreviate the directory `"~/Desktop/lab/emergent-coexistence/analysis/"` into `"analysis/`, otherwise the full path is specified if I refer to another directory.

```
> folder_script <- "~/Desktop/lab/emergent-coexistence/analysis/" 
> folder_pipeline <- "~/Dropbox/lab/emergent-coexistence/plate_scan_pipeline/" 
> folder_data <- "~/Dropbox/lab/emergent-coexistence/data/"
```


## Step 0.2 Generating folder structure and mapping files

Once the directories are specified, navigate to `analysis/` and execute the following scripts to set up the subfolders in `folder_pipeline` for image processing pipeline, mapping files, and ID for pairs and cocultures.

```
$ Rscript 00a-folder_structure.R
$ Rscript 00b-generate_mapping_files.R
$ Rscript 00c-generate_pairs_ID.R
```

Two groups of mapping files are generated:

- `00-list_images-BATCH-CHANNEL.csv`: is used for image processing. Each row represents one image file, and the columns specify the directory where temporary image files are stored.
- `00-list_image_mapping-BATCH.csv`: is used for matching coculture to monocultures. Each row is one coculture pair and the columns specify the batch, community, isolates, mixing frequencies, and the image file name of both isolates.


## Step 1. Command-line tools

These scripts are wrapped into command-line tools that takes the mapping files stored in `analysis/mapping_files/` as input. All temporary output images and data are stored in the subfolders under `folder_pipeline`.

For instance, implementing the image processing pipeline and random forest classification for all cocultures in the batch B2 requires executing the following scripts in order.

```
$ cd analysis/
$ Rscript 01-channel.R mapping_files/00-list_images-B2-red.csv
$ Rscript 01-channel.R mapping_files/00-list_images-B2-green.csv
$ Rscript 01-channel.R mapping_files/00-list_images-B2-blue.csv
$ python 02-rolling_ball.py mapping_files/00-list_images-B2-red.csv
$ python 02-rolling_ball.py mapping_files/00-list_images-B2-green.csv
$ python 02-rolling_ball.py mapping_files/00-list_images-B2-blue.csv
$ Rscript 03-segmentation.R mapping_files/00-list_images-B2-green.csv
$ Rscript 04-feature.R mapping_files/00-list_images-B2-red.csv
$ Rscript 04-feature.R mapping_files/00-list_images-B2-green.csv
$ Rscript 04-feature.R mapping_files/00-list_images-B2-blue.csv
$ Rscript 04a-merge_features.R mapping_files/00-list_images-B2-green.csv
$ Rscript 05-random_forest.R mapping_files/00-list_images-B2-green.csv mapping_files/00-list_image_mapping-B2.csv
```

Below is the overview for the image processing pipeline

![](plots/cartoons/image_processing.png)


## Step 2. Data wrangling and analysis

These scripts take data from either the 16S sequences, or those data generated from the command-line as described above. These data are clean up and outputed into the folder `~/labdata/temp/`

The output data files are

```
$ cd analysis/
$ Rscript 11-align_isolate_sequences.R
$ Rscript 12-assign_isolate_RDP.R
$ Rscript 13-match_community_abundance.R
$ Rscript 14-pairwise_16s_mismatch.py
$ Rscript 15-samebug_pairs.R
$ Rscript 16-match_pair_RDP.R
```


## Step 3. Generating the figures and supplementary pdfs


The main figures Fig.1-4 and supplementary figures Fig.S4-6 are genearting using the following scripts.

```
$ cd analysis/
$ Rscript 96-figures.R
$ Rscript 96a-supp_figures.R
```

Cartoons and Fig.S1-2 are generated using Adobe Illustrator.

The four supplementary pdfs that contain the images and random forest results are generated using the script. The command-line function `convert` is from `imagemagick`.

```
$ cd analysis/
$ Rscript 97-combine_images_and_random_forest.R

# Once the individual pngs are generated, merge them into multi-page pdfs
$ cd ~/Dropbox/lab/emergent-coexistence/plate_scan_pipeline/random_forest/
$ convert -quality 60 B2_*.png random_forest-B2.pdf
$ convert -quality 60 C_*.png random_forest-C.pdf
$ convert -quality 60 C2_*.png random_forest-C2.pdf
$ convert -quality 60 D_*.png random_forest-D.pdf
```

## Sum up

To execute all steps decribed above, from the raw data to ready-for-paper figures, basically run all scripts using terminal commands saved in a master shell script `analysis/00e-commands.sh`. Note that for the shell script to work, the working directory has to be the project directory (where `emergent-coexistence.Rproj` is located)

```
$ zsh analysis/00e-commands.sh
```



## Installing Python and R package dependency

The scripts are based in Mac environment

### Python

We use `pipenv` to keep track the python package dependency

```sh
$ pip install pipenv
$ brew install pipenv # Alternative for Mac users
```

Browse to the local directory and install dependency from Pipfile. Make sure that the Pipfile is in the current directory.

```sh
$ cd <your_local_directory>
$ pipenv install
```

Check dependency. ecoprospector should depend on community-selection and all other dependent packages.

```sh
$ pipenv graph
```

To activate the Pipenv shell:

```sh
$ pipenv shell
$ exit
```

### R

We use `renv` to keep track the R package dependency. This project is operated under the latest R `4.2.0`

```
> sessionInfo()
R version 4.2.0 (2022-04-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur 11.5.2

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices datasets  utils     methods   base     

loaded via a namespace (and not attached):
[1] compiler_4.2.0      BiocManager_1.30.18 tools_4.2.0         renv_0.15.5    
```

To install the package dependency, we use `renv` to record the packages used in this project. After cloning this repository into the local directory, open your code editor (e.g., Rstudio) and run the following two lines in R console. 

```
> install.packages("renv")
> renv::restore()
```

You will be prompt to confirm the installation. This will automatically install all packages on which this project depends. It may takes a few minutes. The installed packages will not be stored in your global environment but instead remain project-specific (saved in the subdirectory `renv/library/`). When you open a new R session under the R project structure (the folder that contains `emergent-coexistence.Rproj`, the directory you decided in step 1), for instance in R studio, these project-specific packages will be already installed. 










