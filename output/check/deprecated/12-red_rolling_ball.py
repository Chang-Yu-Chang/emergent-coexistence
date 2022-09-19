"""
This scripts takes the grey-scale image and apply rooling ball algorithm to subtract backgound

The details and examples are in https://scikit-image.org/docs/stable/auto_examples/segmentation/plot_rolling_ball.html

The rolling ball algorithm only works on grey-scale images
1. invert the image because the image has bright background
2. rolling ball
3. invert the result

To use this script, in terminal `python 02-rooling_ball.py LIST_IMAGES.csv`
"""

import os
import sys
import pandas as pd
import imageio
import skimage
from skimage import data, restoration, util, io, color

# The input has to be the csv generated from 00-list_images.R
list_images = pd.read_csv(str(sys.argv[1]))
#list_images = pd.read_csv('/Users/chang-yu/Desktop/Lab/emergent-coexistence/output/check/00-list_images-D.csv')

def rolling_ball_light(image):
    # invert the image
    image_inverted = util.invert(image)
    
    # rolling ball
    background_inverted = restoration.rolling_ball(image_inverted, radius = 80, num_threads = 10)
    
    # invert the result
    image_rolled_inverted = image_inverted - background_inverted
    image_rolled = util.invert(image_rolled_inverted)
    return image_rolled


for i in range(list_images.shape[0]):
    # File directory
    file_gray = list_images.iloc[i]['folder_red'] + list_images.iloc[i]['image_name'] + '.tiff'
    file_rolled = list_images.iloc[i]['folder_red_rolled'] + list_images.iloc[i]['image_name'] + '.tiff'
    
    # 
    image = io.imread(file_gray)
    image_rolled = rolling_ball_light(image)
    io.imsave(file_rolled, image_rolled)
    print("red channel\t" + str(i) + "/" + str(list_images.shape[0]) + "\t" + list_images.iloc[i]['image_name'])
    
    




