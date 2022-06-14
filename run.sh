echo  '${Running script}', 01A-isolates_experiments-01-metadata.R; Rscript output/report/script/01A-isolates_experiments-01-metadata.R
echo  '${Running script}', 01B-isolate_sanger_seq-01-read_align.R; Rscript output/report/script/01B-isolate_sanger_seq-01-read_align.R
echo  '${Running script}', 01B-isolate_sanger_seq-02-RDP.R; Rscript output/report/script/01B-isolate_sanger_seq-02-RDP.R
echo  '${Running script}', 01B-isolate_sanger_seq-03-make_tree.R; Rscript output/report/script/01B-isolate_sanger_seq-03-make_tree.R
echo  '${Running script}', 01C-isolates_OD_CFU-01-read_OD.R; Rscript output/report/script/01C-isolates_OD_CFU-01-read_OD.R
echo  '${Running script}', 01C-isolates_OD_CFU-02-read_CFU.R; Rscript output/report/script/01C-isolates_OD_CFU-02-read_CFU.R
echo  '${Running script}', 01C-isolates_OD_CFU-03-epsilon.R; Rscript output/report/script/01C-isolates_OD_CFU-03-epsilon.R
echo  '${Running script}', 01D-isolates_growth_rate-01-read_data.R; Rscript output/report/script/01D-isolates_growth_rate-01-read_data.R
echo  '${Running script}', 01E-match_community_abundance-01-community_sequence.R; Rscript output/report/script/01E-match_community_abundance-01-community_sequence.R
echo  '${Running script}', 01E-match_community_abundance-02-match_isolate_16S.R; Rscript output/report/script/01E-match_community_abundance-02-match_isolate_16S.R
echo  '${Running script}', 01E-match_community_abundance-03-matched_abundance.R; Rscript output/report/script/01E-match_community_abundance-03-matched_abundance.R
echo  '${Running script}', 02A-pairs_experiment-01-pairwise_manual_key_in.R; Rscript output/report/script/02A-pairs_experiment-01-pairwise_manual_key_in.R
echo  '${Running script}', 02A-pairs_experiment-02-colony_count.R; Rscript output/report/script/02A-pairs_experiment-02-colony_count.R
echo  '${Running script}', 02A-pairs_experiment-03-pairwise_ambiguous.R; Rscript output/report/script/02A-pairs_experiment-03-pairwise_ambiguous.R
echo  '${Running script}', 02A-pairs_experiment-04-pairwise_plate_layout.R; Rscript output/report/script/02A-pairs_experiment-04-pairwise_plate_layout.R
echo  '${Running script}', 02B-CASEU_sanger_seq-00-test.R; Rscript output/report/script/02B-CASEU_sanger_seq-00-test.R
echo  '${Running script}', 02B-CASEU_sanger_seq-01-pilot1.R; Rscript output/report/script/02B-CASEU_sanger_seq-01-pilot1.R
echo  '${Running script}', 02B-CASEU_sanger_seq-02-pilot2.R; Rscript output/report/script/02B-CASEU_sanger_seq-02-pilot2.R
echo  '${Running script}', 02B-CASEU_sanger_seq-03-pilot3.R; Rscript output/report/script/02B-CASEU_sanger_seq-03-pilot3.R
echo  '${Running script}', 02B-CASEU_sanger_seq-04-pilot4.R; Rscript output/report/script/02B-CASEU_sanger_seq-04-pilot4.R
echo  '${Running script}', 02B-CASEU_sanger_seq-05-six_plates.R; Rscript output/report/script/02B-CASEU_sanger_seq-05-six_plates.R
echo  '${Running script}', 02C-pairs_OD_CFU-01-epsilon_uncertainty.R; Rscript output/report/script/02C-pairs_OD_CFU-01-epsilon_uncertainty.R
echo  '${Running script}', 02C-pairs_OD_CFU-02-CFU_frequency.R; Rscript output/report/script/02C-pairs_OD_CFU-02-CFU_frequency.R
echo  '${Running script}', 02C-pairs_OD_CFU-03-CFU_frequency_uncertainty.R; Rscript output/report/script/02C-pairs_OD_CFU-03-CFU_frequency_uncertainty.R


# Run the following scripts in bash

echo  '${Running script}', 02D-determine_pairwise_interaction-01-combine_CFU_CASEU_result.R; Rscript output/report/script/02D-determine_pairwise_interaction-01-combine_CFU_CASEU_result.R
echo  '${Running script}', 02D-determine_pairwise_interaction-02-determine_pairwise_interaction.R; Rscript output/report/script/02D-determine_pairwise_interaction-02-determine_pairwise_interaction.R
echo  '${Running script}', 02D-determine_pairwise_interaction-03-isolate_tournament.R; Rscript output/report/script/02D-determine_pairwise_interaction-03-isolate_tournament.R
echo  '${Running script}', 02E-competition_phylogeny-01-pairs_taxonomy.R; Rscript output/report/script/02E-competition_phylogeny-01-pairs_taxonomy.R
echo  '${Running script}', 02E-competition_phylogeny-02-pairs_16S.R; Rscript output/report/script/02E-competition_phylogeny-02-pairs_16S.R
echo  '${Running script}', 03A-make_network-01-make_pairwise_network.R; Rscript output/report/script/03A-make_network-01-make_pairwise_network.R
echo  '${Running script}', 03A-make_network-02-randomize_network.R; Rscript output/report/script/03A-make_network-02-randomize_network.R
echo  '${Running script}', 03B-network_measure-01-detect_motif.R; Rscript output/report/script/03B-network_measure-01-detect_motif.R
echo  '${Running script}', 03B-network_measure-02-competitive_hierarchy.R; Rscript output/report/script/03B-network_measure-02-competitive_hierarchy.R
echo  '${Running script}', 03B-network_measure-03-distance_to_diagonal.R; Rscript output/report/script/03B-network_measure-03-distance_to_diagonal.R
echo  '${Running script}', 03B-network_measure-04-node_degree.R; Rscript output/report/script/03B-network_measure-04-node_degree.R
echo  '${Running script}', 03B-network_measure-05-component_count.R; Rscript output/report/script/03B-network_measure-05-component_count.R
echo  '${Running script}', 04A-summary-01-isolates.R; Rscript output/report/script/04A-summary-01-isolates.R
echo  '${Running script}', 04B-summary-01-pairs.R; Rscript output/report/script/04B-summary-01-pairs.R
echo  '${Running script}', 04B-summary-02-outcome_types.R; Rscript output/report/script/04B-summary-02-outcome_types.R
