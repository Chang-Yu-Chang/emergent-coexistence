# Simulating community assembly and pairwise competition using consumer-resource models

This README describes the dependency environment and commands to reproduce the simulation 

## Dependency

Install community-simulator

```
cd simulation/community-simulator/
pip install -e .

pip3 install osqp==0.6.1
```

Enter the pipenv environment

```
pipenv shell
```

In R environment

```
> sessionInfo()
R version 4.2.2 (2022-10-31)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Ventura 13.2.1

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] cowplot_1.1.1   forcats_1.0.0   stringr_1.5.0   dplyr_1.1.0     purrr_1.0.1     readr_2.1.4    
 [7] tidyr_1.3.0     tibble_3.1.8    ggplot2_3.4.1   tidyverse_1.3.2

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.0    haven_2.5.1         gargle_1.3.0        colorspace_2.1-0    vctrs_0.5.2        
 [6] generics_0.1.3      utf8_1.2.3          rlang_1.0.6         pillar_1.8.1        glue_1.6.2         
[11] withr_2.5.0         DBI_1.1.3           bit64_4.0.5         dbplyr_2.3.0        modelr_0.1.10      
[16] readxl_1.4.2        lifecycle_1.0.3     munsell_0.5.0       gtable_0.3.1        cellranger_1.1.0   
[21] ragg_1.2.5          rvest_1.0.3         labeling_0.4.2      tzdb_0.3.0          parallel_4.2.2     
[26] fansi_1.0.4         broom_1.0.3         scales_1.2.1        backports_1.4.1     googlesheets4_1.0.1
[31] vroom_1.6.1         jsonlite_1.8.4      systemfonts_1.0.4   farver_2.1.1        fs_1.6.1           
[36] bit_4.0.5           textshaping_0.3.6   hms_1.1.2           stringi_1.7.12      grid_4.2.2         
[41] rprojroot_2.0.3     here_1.0.1          cli_3.6.0           tools_4.2.2         magrittr_2.0.3     
[46] crayon_1.5.2        pkgconfig_2.0.3     ellipsis_0.3.2      xml2_1.3.3          reprex_2.0.2       
[51] googledrive_2.0.0   lubridate_1.9.2     timechange_0.2.0    assertthat_0.2.1    httr_1.4.4         
[56] rstudioapi_0.14     R6_2.5.1            compiler_4.2.2     
```

To ensure the commandline tool `Rscript` has the same installed R packages as in the R environment, include the following line to edit `PATH` in `~/.bashrc` (for bash) or in `~/.zshrc` (for zsh).

```
# Change the terminal R PATH such that it loads the Rstudio R version with the packages installed
export PATH="/Library/Frameworks/R.framework/Resources:$PATH"
```


## Simulation

The is parameterized by an universal input file `01-input_parameters.csv`

```
Rscript 01-generate_input.R
```

Generate the mapping files for 1) monocultures and 2) top-down communities. 

```
# Generate the mapping files
Rscript 02-generate_input_mono_comm.R
# Run simulations
python3 run.py 02a-input_monocultures.csv 0
python3 run.py 02b-input_communities.csv 0
python3 run.py 02c-input_communitiesWithoutCrossfeeding.csv 0
```

After the community and monoculture simulation is done

```
# Generate the mapping files
Rscript 03-generate_input_pairs.R
# Run simulations
zsh 03a-run_poolPairs.sh
zsh 03b-run_withinCommunityPairs.sh
```

Aggregate the simulation data

```
Rscript 12-aggregate_mono_comm.R
Rscript 13-aggregate_pairs.R
Rscript 14-pairwise_network.R
```


Make figures from the aggregated data

```
Rscript 21-matrices.R               # c, D. l matrices
Rscript 22-visualize_mono_comm.R    # community composition, monocultures
Rscript 23-visualize_pairs.R        # pairwise competition
Rscript 24-visualize_networks.R     # pairwise networks
```






























