#!/usr/bin/env bash

# Make synthetic communities
make_synthetics() {
    # Culturable isolates
    Rscript make_culturable_isolates.R $1
    # Pairs from culturable isolates
    Rscript make_culturable_pairs.R $1
    # Trios from culturable isolates
    Rscript make_culturable_trios.R $1
    # Pairs from trios
    Rscript make_culturable_pair_from_trio.R $1
    # Pairs from top-down communities
    Rscript make_pair_from_top_down_community.R $1
}


