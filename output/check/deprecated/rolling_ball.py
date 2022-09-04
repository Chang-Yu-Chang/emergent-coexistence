"""
This scripts takes the grey-scale image and apply rooling ball algorithm to subtract backgound

The details and examples are in https://scikit-image.org/docs/stable/auto_examples/segmentation/plot_rolling_ball.html

The rolling ball algorithm only works on grey-scale images
1. invert the image because the image has bright background
2. rolling ball
3. invert the result

"""

import imageio
import os
import sys
import skimage
from skimage import data, restoration, util, io, color

file_gray = str(sys.argv[1])
file_rolled = str(sys.argv[2])

def rolling_ball_light(image):
    # 1. invert the image
    image_inverted = util.invert(image)
    
    # 2. rolling ball
    background_inverted = restoration.rolling_ball(image_inverted, radius = 80, num_threads = 10)
    
    # 3. invert the result
    image_rolled_inverted = image_inverted - background_inverted
    image_rolled = util.invert(image_rolled_inverted)
    #background = util.invert(background_inverted)
    return image_filtered
    
image = io.imread(file_gray)
image_rolled = rolling_ball_light(image)
io.imsave(file_rolled, image_rolled)

















