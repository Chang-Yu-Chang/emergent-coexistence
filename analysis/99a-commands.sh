cd ~/Desktop/lab/emergent-coexistence/analysis/
# Batch B2
Rscript 01-channel.R mapping_files/00-list_images-B2-red.csv
Rscript 01-channel.R mapping_files/00-list_images-B2-green.csv
Rscript 01-channel.R mapping_files/00-list_images-B2-blue.csv
python 02-rolling_ball.py mapping_files/00-list_images-B2-red.csv
python 02-rolling_ball.py mapping_files/00-list_images-B2-green.csv
python 02-rolling_ball.py mapping_files/00-list_images-B2-blue.csv
Rscript 03-segmentation.R mapping_files/00-list_images-B2-green.csv
Rscript 04-feature.R mapping_files/00-list_images-B2-red.csv
Rscript 04-feature.R mapping_files/00-list_images-B2-green.csv
Rscript 04-feature.R mapping_files/00-list_images-B2-blue.csv
Rscript 04a-merge_features.R mapping_files/00-list_images-B2-green.csv
Rscript 05-random_forest.R mapping_files/00-list_images-B2-green.csv mapping_files/00-list_image_mapping-B2.csv
# Batch C
Rscript 01-channel.R mapping_files/00-list_images-C-red.csv
Rscript 01-channel.R mapping_files/00-list_images-C-green.csv
Rscript 01-channel.R mapping_files/00-list_images-C-blue.csv
python 02-rolling_ball.py mapping_files/00-list_images-C-red.csv
python 02-rolling_ball.py mapping_files/00-list_images-C-green.csv
python 02-rolling_ball.py mapping_files/00-list_images-C-blue.csv
Rscript 03-segmentation.R mapping_files/00-list_images-C-green.csv
Rscript 04-feature.R mapping_files/00-list_images-C-red.csv
Rscript 04-feature.R mapping_files/00-list_images-C-green.csv
Rscript 04-feature.R mapping_files/00-list_images-C-blue.csv
Rscript 04a-merge_features.R mapping_files/00-list_images-C-green.csv
Rscript 05-random_forest.R mapping_files/00-list_images-C-green.csv mapping_files/00-list_image_mapping-C.csv
# Batch C2
Rscript 01-channel.R mapping_files/00-list_images-C2-red.csv
Rscript 01-channel.R mapping_files/00-list_images-C2-green.csv
Rscript 01-channel.R mapping_files/00-list_images-C2-blue.csv
python 02-rolling_ball.py mapping_files/00-list_images-C2-red.csv
python 02-rolling_ball.py mapping_files/00-list_images-C2-green.csv
python 02-rolling_ball.py mapping_files/00-list_images-C2-blue.csv
Rscript 03-segmentation.R mapping_files/00-list_images-C2-green.csv
Rscript 04-feature.R mapping_files/00-list_images-C2-red.csv
Rscript 04-feature.R mapping_files/00-list_images-C2-green.csv
Rscript 04-feature.R mapping_files/00-list_images-C2-blue.csv
Rscript 04a-merge_features.R mapping_files/00-list_images-C2-green.csv
Rscript 05-random_forest.R mapping_files/00-list_images-C2-green.csv mapping_files/00-list_image_mapping-C2.csv
# Batch D
Rscript 01-channel.R mapping_files/00-list_images-D-red.csv
Rscript 01-channel.R mapping_files/00-list_images-D-green.csv
Rscript 01-channel.R mapping_files/00-list_images-D-blue.csv
python 02-rolling_ball.py mapping_files/00-list_images-D-red.csv
python 02-rolling_ball.py mapping_files/00-list_images-D-green.csv
python 02-rolling_ball.py mapping_files/00-list_images-D-blue.csv
Rscript 03-segmentation.R mapping_files/00-list_images-D-green.csv
Rscript 04-feature.R mapping_files/00-list_images-D-red.csv
Rscript 04-feature.R mapping_files/00-list_images-D-green.csv
Rscript 04-feature.R mapping_files/00-list_images-D-blue.csv
Rscript 04a-merge_features.R mapping_files/00-list_images-D-green.csv
Rscript 05-random_forest.R mapping_files/00-list_images-D-green.csv mapping_files/00-list_image_mapping-D.csv
