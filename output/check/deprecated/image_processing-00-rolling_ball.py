"""
This scripts takes the grey-scale image and apply rooling ball algorithm to subtract backgound

Only applied to the csv written in image_processing-01-example.csv

"""

import os
import sys
import pandas as pd
import imageio
import skimage
from skimage import data, restoration, util, io, color

list_images = pd.read_csv('/Users/chang-yu/Desktop/lab/emergent-coexistence/output/check/image_processing-01-example.csv')

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
    file_gray = list_images.iloc[i]['folder_green'] + list_images.iloc[i]['image_name'] + '.tiff'
    file_rolled = list_images.iloc[i]['folder_rolled'] + list_images.iloc[i]['image_name'] + '.tiff'
    
    # 
    image = io.imread(file_gray)
    image_rolled = rolling_ball_light(image)
    io.imsave(file_rolled, image_rolled)
    print("rolled\t" + str(i) + "/" + str(list_images.shape[0]) + "\t" + list_images.iloc[i]['image_name'])
    

















