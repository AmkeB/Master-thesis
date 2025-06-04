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



filepath = r"\\store.erasmusmc.nl\department\gene\chien_data\Lab\Data_and_Analysis\Amke Bernaerts\chapter4_thesis_image_analysis\image registration\AB009"

files = glob(filepath+"\\parts_images\\*.tif")
files.sort()


flip_even = True
flip_uneven = False


PARTS_startpoint = [905,160] #top left in parts image, x,y as defined in fiji >> for AB009
parts_pixelsize = 1.376 #um/pixel #for series3

channels = [460,532,640]

if len(files) != len(channels):
    print('error: number of channels given does not correspond to number of parts files')
    

# here the PARTS images are loaded for each channel and cropped into a snake-grid with same size ROIs as UFO2
for f in range(len(files)):
    parts_file = files[f]
    parts_img = cv2.imread(parts_file, cv2.IMREAD_UNCHANGED)
    
    parts_ROI_size = int((5120*0.31295)/parts_pixelsize)
    
    Nrows = int(np.floor(((np.shape(parts_img)[0])-PARTS_startpoint[1])/parts_ROI_size))
    Ncolumns = int(np.floor(((np.shape(parts_img)[1])-PARTS_startpoint[0])/parts_ROI_size))
    total_ROIs = Nrows*Ncolumns
    
    #making a snake-shaped grid on the PARTS image
    
    #first specify top right coordinates of all ROIs in the NrowsxNcolumns grid:
    top_left_x_coords = np.zeros([Nrows,Ncolumns])
    top_left_y_coords = np.zeros([Nrows,Ncolumns])
    for i in range(Nrows):
        for j in range(Ncolumns):
            top_left_x_coords[i,j] = PARTS_startpoint[0]+parts_ROI_size*j
            top_left_y_coords[i,j] = PARTS_startpoint[1]+parts_ROI_size*i
            
    snake_grid = np.reshape(np.arange(total_ROIs-1,-1,-1),(Ncolumns,Nrows)).T
    for i in range(Ncolumns):
        if flip_even == True:
            if (i%2)==0:
                column_flipped = np.flip(snake_grid[:,i])
                snake_grid[:,i] = column_flipped
        elif flip_uneven == True:
            if (i%2)!=0:
                column_flipped = np.flip(snake_grid[:,i])
                snake_grid[:,i] = column_flipped
            
   # create folder to store created ROIs for each channel
    try:
        os.mkdir(filepath+"\\%03dnm_parts_ROIs" % channels[f])
    except:
        print("folder already exists")
    
    for i in range(Nrows):
        for j in range(Ncolumns):
            ROI = snake_grid[i,j]
            x = int(top_left_x_coords[i,j])
            y = int(top_left_y_coords[i,j])
            ROI_parts = parts_img[y:(y+parts_ROI_size),x:(x+parts_ROI_size)]
            # ROI_parts_rot = cv2.rotate(ROI_parts, cv2.ROTATE_90_CLOCKWISE) #rotate 90 degrees to right to be same as UFO2
            
            # cv2.imwrite(filepath+"\PARTS_ROIs\parts_ROI%03d.tif" % ROI,ROI_parts)
            cv2.imwrite(filepath+"\\%03dnm_parts_ROIs"%channels[f]+"\\ROI%03d.tif" % ROI , ROI_parts)
            
    #creating a ROI.txt file for image stitching
    ROI_txt = np.zeros([total_ROIs,2])
    for i in range(Nrows):
        for j in range(Ncolumns):
            # print(j)
            ROI_number = snake_grid[i,j]
            ROI_txt[ROI_number,1] = top_left_x_coords[i,j]
            ROI_txt[ROI_number,0] = top_left_y_coords[i,j]
    # np.savetxt(filepath+"\ROI.txt", ROI_txt, fmt='%d', delimiter='\t')
    np.savetxt(filepath+"\ROI.txt", ROI_txt, fmt='%d', delimiter='\t')
        

    












