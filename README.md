
Scripts for the experimental project on emergent coexistence

Preprint link: [Emergent coexistence in multispecies microbial communities](https://www.biorxiv.org/content/10.1101/2022.05.20.492860v1)

## Setup to repeat the analysis

### Step 1: clone this repository

Change the current working directory to the location where you want the cloned directory.

```
$ cd /directory/to/where/you/want/ # Directory to where you want
$ git clone https://github.com/Chang-Yu-Chang/emergent-coexistence
```

### Step 2: sync the raw data to Dropbox

Link to the raw data 

https://www.dropbox.com/sh/60rw7xmerown3vk/AACMzcZMR-DzF_5vG3KF0CG-a?dl=0

The data are currently saved in the directory `~/Dropbox/lab/emergent-coexistence/data/`, where the all scripts are written to read this absolute directory. For the following scripts to work, this directory has to be set up.


## Data

The data wrangling steps are conducted in order in four master Rmd along with associated R scripts in `output/report/script/`: 

- `output/report/01-isolates.Rmd`
- `output/report/02-pairs.Rmd`
- `output/report/03-networks.Rmd`
- `output/report/04-summary.Rmd`

A quick overview of these steps is reported in the generated notebooks `01-isolates.nb.html`, `02-pairs.nb.html`, `03-networks.nb.html`, and `04-summary.nb.html`. 

To generate all ready-for-figure data in `~/Dropbox/lab/emergent-coexistence/data/output/`, basically run all R scripts in `output/report/script/` in order, either by executing the code chunks in the four Rmd files, or using terminal commands saved in a master zsh.

```
$ cd output/report/script/
$ zsh run.zsh
```

## Figures

The R scripts used to generate the figures shown in the paper is in `plotting_scripts/`

## To be continuted once the figures are finalized











