# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import cv2
from scipy.ndimage import gaussian_filter, binary_dilation, binary_closing, binary_erosion, binary_opening
from copy import deepcopy
import gc
# import math
#from skimage.draw import circle_perimeter
# import time
import os
# import re
from glob import glob

def downsize_image(image, factor_y=10, factor_x=10, img_type=np.uint8):
    """
    Here we downsize the image from a nxm to an n/f x m/f size with f being the downscale factor. If n % f != 0 and/or
    m % f != 0, the remainder pixels of the nxm image are not used in the downsizing to the n/f x m/f image.

    Parameters:
        image (ndarray):    The input image.
        factor (float):     The factor of downscaling.
        img_type (numpy):   The image type you want it in (not what it is).
                            
    Returns:
        new_img (ndarray):  The downscaled image.

    """
    if image.dtype == img_type:
        temp_img = np.zeros((int(image.shape[0] / factor_y), image.shape[1]), dtype=img_type)
        new_img = np.zeros((int(image.shape[0] / factor_y), int(image.shape[1] / factor_x)), dtype=img_type)

        for e in range(temp_img.shape[0]):
            temp_img[e, :] = np.mean(image[int(factor_y * e):int(factor_y * e + factor_y), :], axis=0)
        for j in range(new_img.shape[1]):
            new_img[:, j] = np.mean(temp_img[:, int(factor_x * j):int(factor_x * j + factor_x)], axis=1)

        return new_img
    
    elif image.dtype == np.uint16 and img_type == np.uint8:
        image = cv2.convertScaleAbs(image, alpha=(255/65535))

        temp_img = np.zeros((int(image.shape[0] / factor_y), image.shape[1]), dtype=img_type)
        new_img = np.zeros((int(image.shape[0] / factor_y), int(image.shape[1] / factor_x)), dtype=img_type)

        for e in range(temp_img.shape[0]):
            temp_img[e, :] = np.mean(image[int(factor_y * e):int(factor_y * e + factor_y), :], axis=0)
        for j in range(new_img.shape[1]):
            new_img[:, j] = np.mean(temp_img[:, int(factor_x * j):int(factor_x * j + factor_x)], axis=1)
    
        return new_img
    else:
        raise TypeError("downsize image() only works with 8 and 16 bit images, please convert it to one of the two first or contact Lars to ask to add your type to the function.")

ds_factor = 1 #10 # image should be perfectly divided in x and y by this factor, if not, errors will occur probably

filepath = r"\\store.erasmusmc.nl\department\gene\chien_data\Lab\Data_and_Analysis\Amke Bernaerts\chapter4_thesis_image_analysis\image registration\AB009"
txt = np.loadtxt(filepath+os.sep+"\ROI.txt", dtype=int)

x_pos = np.unique(txt[:,1])
y_pos = np.unique(txt[:,0])

idx_x = np.zeros(len(txt), dtype=int)
idx_y = np.zeros(len(txt), dtype=int)

N_cols = len(x_pos)
N_rows = len(y_pos)

channels = [460,532,640]

#For AB009: these are the ROIs that need to be removed to have the same boundaries removed as UFO2
rm_i = [5,6,15,16,25,26,35,36,45,46,55,56,65,66,75,76,85,86,95,96,105,106,115,116,125,126,135,136,145,146,155,156,165,166,175,176,185,186,195,196,205,206]
rm_i = np.append(rm_i, np.arange(41,61))
rm_i = np.append(rm_i, np.arange(101,111))
rm_i = np.append(rm_i, np.arange(151,171))
rm_i = rm_i-1

BIG_imgs = []

for c in channels:
    path_one_channel = filepath+"\\%03dnm_parts_ROIs" % c

    # WL_files = glob(filepath+"\\WhiteLight\\*.tif")
    WL_files = glob(path_one_channel+"\\*.tif")
    WL_files.sort()
    
    # size_image = 5120
    size_image = 1164
    BIG_img = np.zeros((int(N_rows*size_image/ds_factor), int(N_cols*size_image/ds_factor)), dtype=np.uint8) # we now assume all the images we use are 5120x5120 pixels
    
    
    try:
        os.mkdir(filepath+"\\%03dnm_parts_ROIs_no_boundaries" % c)
    except:
        print("folder already exists")
        
        
    ROI=0
    for i in range(len(txt)):
    # for i in range(128):
        img = cv2.imread(WL_files[i], cv2.IMREAD_UNCHANGED)
        # img = cv2.convertScaleAbs(img, alpha=(255/65535))
        if i not in rm_i:
            img_rot = cv2.rotate(img, cv2.ROTATE_90_CLOCKWISE) #to get parts images in same orientation as ufo2 images
            cv2.imwrite(filepath+"\\%03dnm_parts_ROIs_no_boundaries" % c +"\\parts_ROI%03d.tif" % ROI,img_rot)
            ROI += 1
            
        img = downsize_image(img, factor_y=ds_factor,factor_x=ds_factor, img_type=np.uint8)
        
        x = txt[i,1]
        y = txt[i,0]
        
        idx = np.where(x_pos==x)[0][0]
        idy = np.where(y_pos==y)[0][0]
        idx_x[i] = idx
        idx_y[i] = idy
        
        # plt.imshow(BIG_img)
        # plt.axis('off')
        # plt.show()
        
        BIG_img[int(idy*size_image/ds_factor):int((idy+1)*size_image/ds_factor), int(idx*size_image/ds_factor):int((idx+1)*size_image/ds_factor)] = img
    
        print("{:0.2f}% done".format((i+1)/len(txt)*100))
        gc.collect()
    
    plt.imshow(BIG_img)
    plt.axis('off')
    plt.show()
    BIG_imgs.append(BIG_img)
    cv2.imwrite(filepath+os.sep+"%03dnm_grid.tif"%c,BIG_img)
    


