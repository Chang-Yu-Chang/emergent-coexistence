"""
This scripts takes the grey-scale image and apply rooling ball algorithm to subtract backgound

The details and examples are in https://scikit-image.org/docs/stable/auto_examples/segmentation/plot_rolling_ball.html

Below are what this script does
1. read a grey-scale image
2. because the image has bright background, invert the image
3. rolling ball
4. invert the result
5. output the file as the same name

To use this script, in terminal `python 02-rooling_ball.py directory_to_the_image`
"""

import imageio
import matplotlib.pyplot as plt
import numpy as np
import pywt
import os
from skimage import data, restoration, util, io, color
import skimage


directory = '/Users/chang-yu/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/check/'
folder_source = directory + 'C2_green/'
list_img = os.listdir(folder_source) 
list_img.sort()

def rolling_ball_light(image):
    # 2. invert the image
    image_inverted = util.invert(image)
    
    # 3. rolling ball
    background_inverted = restoration.rolling_ball(image_inverted, radius = 80, num_threads = 10)
    
    # 4. invert the result
    filtered_image_inverted = image_inverted - background_inverted
    filtered_image = util.invert(filtered_image_inverted)
    #background = util.invert(background_inverted)
    return filtered_image
    

for file_name in list_img:
    # 1. read the input image
    full_name = os.path.join(folder_source, file_name)
    #file_name = 'C_T8_C11R1_1' + ".tiff"
    #full_name = os.path.join(folder_source, file_name)
    image = io.imread(full_name)#[:1000, :1000]
    print(file_name)
    print(image.shape)
    
    # 2-4. invert the image
    filtered_image = rolling_ball_light(image)

    # 5. output the image file
    folder_output = folder_source
    io.imsave(folder_output + file_name, filtered_image)


# Test indivual image
if False:
    file_name = 'C_T8_C11R1_1' + ".tiff"
    full_name = os.path.join(folder_source, file_name)
    image = io.imread(full_name)#[:1000, :1000]
    print(file_name)
    print(image.shape)
    
    # 2-4. invert the image
    filtered_image = rolling_ball_light(image)
    



