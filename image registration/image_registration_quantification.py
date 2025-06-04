# -*- coding: utf-8 -*-
"""
Created on Fri May 16 12:00:15 2025

@author: 106667
"""

import matplotlib.pyplot as plt
import cv2
from scipy.ndimage import gaussian_filter, binary_dilation, binary_closing, binary_erosion, binary_opening
from copy import deepcopy
# import gc
import math
# from scipy.signal import fftconvolve
# import logging
# import re
import numpy as np
# import time
# import os
from glob import glob

# import multiprocessing
# from AB010_Lars_image_analysis_functions import load_images_and_find_shift_for_PARTS_v3, auto_blur_contrast #, tophat, auto_blur_contrast, find_shift

import cv2
# import matplotlib.pyplot as plt

from skimage.metrics import normalized_mutual_information, hausdorff_distance
# from sklearn.metrics import jaccard_similarity_score

def auto_blur_contrast(image, minmax=False, thr_perc=0.4, blur=1):
    """
    This function applies an automatic blurring and contrasting to enhance the contrast the most of the input image. 
    Parameters can be changed to fine-tune for own application, or just the min and max intensity value can be used for contrast scaling.

    Parameters
    ----------
    image: (ndarray, dtype=uint)
        input image, can be grayscale or rgb.
    minmax: (bool), optional
        If True, only the minimum and maximum value of the image is used for the contrast enhancement. The default is False.
    thr_perc : float, optional
        the percentage on which the cutoff for contrast enchancement should be. 
        0.1*thr_perc for the lower percentage limit and thr_perc for the upper limit. 
        The default is 0.4 and the range 0 to 10 is recommended. 
    blur : int, optional
        The value for the sigma in the gaussian_filter blurring. The default is 1.

    Returns
    -------
    img : (ndarray, dtype=uint)
        The contrast enhanced and (slightly) blurred output image.
    """

    img = deepcopy(image)

    dt = img.dtype
    m = 255
    if dt != np.uint8:
        if dt == np.uint16:
            m = 65535
        elif dt == np.uint32:
            m = 4294967295
        elif dt == np.uint64:
            m = 18446744073709551616
        else:
            m = 65535
            dt = np.uint16  # if it is not uint yet, we're going to make it uint16 in the end
            #raise TypeError("Only uint images allowed")

    img = gaussian_filter(img, sigma=blur)

    if minmax:
        thr_min = np.min(img)
        thr_max = np.max(img)

    else:
        color = False
        if len(img.shape) == 3:
            img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
            color = True
        # elif img.dtype != np.uint8:
        #     img = cv2.convertScaleAbs(img, alpha=(255/m))

        histo, bins = np.histogram(img, bins=20000)
        histo = np.cumsum(histo)
        histo = histo / histo[-1] * 100

        min_where = np.where(histo < 0.1*thr_perc)
        max_where = np.where(histo > 100-thr_perc)

        if len(min_where[0]) != 0:
            plc_min = min_where[0][-1]
            thr_min = math.ceil((bins[plc_min] / 2 + bins[plc_min + 1] / 2))
        else:
            #print("no plc in histogram found to be <0.5%")
            thr_min = np.min(image)

        if len(max_where[0]) != 0:
            plc_max = max_where[0][0]
            thr_max = math.floor((bins[plc_max] / 2 + bins[plc_max + 1] / 2))
        else:
            #print("no plc in histogram found to be >99.5%")
            thr_max = np.max(image)

        if color:
            img = gaussian_filter(image, sigma=blur)

    #print("thr_min = {}".format(thr_min))
    #print("thr_max = {}".format(thr_max))
    img[img < thr_min] = thr_min
    img[img > thr_max] = thr_max
    # print(thr_max,thr_min)
    # if thr_max == 0 and thr_min == 0:
    #     thr_max = 1
    img = (img - thr_min) * np.float64(m/(thr_max-thr_min))
    img[img % 1 < 0.5] = np.floor(img[img % 1 < 0.5])
    img[img % 1 >= 0.5] = np.ceil(img[img % 1 >= 0.5])

    img = np.array(img, dtype=dt)

    return img

filepath_unshifted_ufo = r"\\store.erasmusmc.nl\department\gene\chien_data\Lab\Data_and_Analysis\Amke Bernaerts\chapter4_thesis_image_analysis\image registration\AB009\460nm_parts_ROIs_no_boundaries\ufo2_unshifted_imgs"
filepath_shifted_ufo_rot = r"\\store.erasmusmc.nl\department\gene\chien_data\Lab\Data_and_Analysis\Amke Bernaerts\chapter4_thesis_image_analysis\image registration\AB009\460nm_parts_ROIs_no_boundaries\ufo2_shifted_imgs_v3 rot"
filepath_shifted_ufo_norot = r"\\store.erasmusmc.nl\department\gene\chien_data\Lab\Data_and_Analysis\Amke Bernaerts\chapter4_thesis_image_analysis\image registration\AB009\460nm_parts_ROIs_no_boundaries\ufo2_shifted_imgs_v3 no rot"

filepath_parts = r"\\store.erasmusmc.nl\department\gene\chien_data\Lab\Data_and_Analysis\Amke Bernaerts\chapter4_thesis_image_analysis\image registration\AB009\460nm_parts_ROIs_no_boundaries\parts_ROIs_for_composite"
 
    
def DICE_COE(mask1, mask2):
    intersect = np.sum(mask1*mask2)
    fsum = np.sum(mask1)
    ssum = np.sum(mask2)
    dice = (2 * intersect ) / (fsum + ssum)
    # print('before',dice)
    # dice = np.mean(dice)
    return dice 

def calculate_dice_coeff(file_parts,file_ufo):

    img1 = cv2.imread(file_parts,cv2.IMREAD_UNCHANGED)
    registered_img2 = cv2.imread(file_ufo,cv2.IMREAD_UNCHANGED)
    
    #perform zeropadding on the regions where there is no image anymore:
    registered_img2_max = np.max(registered_img2)
    registered_img2[registered_img2 >= registered_img2_max] = 0
        
    
    
    # img1_blur = auto_blur_contrast(img1,thr_perc=0.4,blur=1)
    # img2_blur = auto_blur_contrast(registered_img2,thr_perc=0.4,blur=1)
    img1_blur = img1
    img2_blur = registered_img2
    ret1,img1_thr = cv2.threshold(img1_blur,0,1,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
    ret2,img2_thr = cv2.threshold(img2_blur,0,1,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
    # img1_thr[img1_thr==255]=1 
    # img2_thr[img2_thr==255]=1
    
    # # mi = normalized_mutual_information(img1_thr, img2_thr)
    # haus = hausdorff_distance(img1_thr,img2_thr)
    # # print(mi)
    # print(haus)
    
    # intersection = img1_thr.intersection(img2_thr)
    # union = img1_thr.union(img2_thr)
    # ji = len(intersection)/len(union)
    # print(ji)
    
    # pearson = np.corrcoef(img1_thr.flatten(),img2_thr.flatten())[0][1]
    # print(pearson)
 
    
    # img1_thr[img1_thr==255]=1 
    # img2_thr[img2_thr==255]=1
    dice = DICE_COE(img1_thr,img2_thr)
    
    # filepath1 = r"\\store.erasmusmc.nl\department\gene\chien_data\Lab\Data_and_Analysis\Amke Bernaerts\AB009_insituFUNseq_washing_optimization\Raw data\AB009\step2_imaging_and_segmentation_20250319at133800\460nm_renamed"
    # store = r"\\store.erasmusmc.nl\department\gene\chien_data\Lab\Data_and_Analysis\Amke Bernaerts\chapter4_thesis_image_analysis\image registration\AB009\460nm_parts_ROIs_no_boundaries\ufo2_unshifted_imgs"
    # WL_files1 = glob(filepath1+"\\*.tif")
    # WL_files1.sort()
    # parts_size=1164
    # for i in range(len(WL_files1)):
    #     file = WL_files1[i]
    #     img1_original = cv2.imread(file,cv2.IMREAD_UNCHANGED)
    #     img1_original = cv2.resize(img1_original,(parts_size,parts_size))
    #     cv2.imwrite(store+"\\ROI%03d.tif"%i,img1_original)
    return dice


ufo_files_unshifted = glob(filepath_unshifted_ufo+"\\*.tif")
ufo_files_unshifted.sort()
ufo_files_shifted_rot = glob(filepath_shifted_ufo_rot+"\\*.tif")
ufo_files_shifted_rot.sort()
ufo_files_shifted_norot = glob(filepath_shifted_ufo_norot+"\\*.tif")
ufo_files_shifted_norot.sort()
parts_files = glob(filepath_parts+"\\*.tif")
parts_files.sort()

dice_unshifted = []
dice_norot = []
dice_rot = []

for i in range(len(parts_files)):
    ufo_unshifted = ufo_files_unshifted[i]
    ufo_rot = ufo_files_shifted_rot[i]
    ufo_norot = ufo_files_shifted_norot[i]
    parts = parts_files[i]
    dice_unshifted.append(calculate_dice_coeff(parts,ufo_unshifted))
    dice_norot.append(calculate_dice_coeff(parts,ufo_norot))
    dice_rot.append(calculate_dice_coeff(parts,ufo_rot))
    
    
# np.savetxt(filepath_unshifted_ufo + "\\dice_shifted_norot_noblur.csv",dice_norot, delimiter=',')  

# Example Dice coefficient arrays (replace with your real data)
dice_unregistered = dice_unshifted[0:20]
dice_registered_no_rotation = dice_norot[0:20]
dice_registered_with_rotation = dice_rot[0:20]

num_images = len(dice_unregistered)
x = np.arange(num_images)  # Index of each image

plt.figure(figsize=(10, 6))

plt.plot(x, dice_unregistered, marker='o', label='Unregistered')
plt.plot(x, dice_registered_no_rotation, marker='s', label='Registered (No Rotation)')
plt.plot(x, dice_registered_with_rotation, marker='^', label='Registered (With Rotation)')

plt.xticks(x, [f'Image {i+1}' for i in x])
plt.xlabel('Image Index')
plt.ylabel('Dice Coefficient')
plt.title('Dice Coefficients Across Registration Conditions')
plt.ylim(0, 1)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()    
bar_width = 0.25
x = np.arange(num_images)

plt.figure(figsize=(10, 6))
plt.bar(x - bar_width, dice_unregistered, width=bar_width, label='Unregistered')
plt.bar(x, dice_registered_no_rotation, width=bar_width, label='Registered (No Rotation)')
plt.bar(x + bar_width, dice_registered_with_rotation, width=bar_width, label='Registered (With Rotation)')

plt.xticks(x, [f'Image {i+1}' for i in x])
plt.xlabel('Image Index')
plt.ylabel('Dice Coefficient')
plt.title('Dice Coefficients Comparison')
plt.ylim(0, 1)
plt.legend()
plt.grid(axis='y')
plt.tight_layout()
plt.show()



# file_parts = filepath_parts+"\\parts_ROI005.tif"
# file_ufo = filepath_shifted_ufo_rot+"\\ROI005.tif"
# img1 = cv2.imread(file_parts,cv2.IMREAD_UNCHANGED)
# registered_img2 = cv2.imread(file_ufo,cv2.IMREAD_UNCHANGED)

# #perform zeropadding on the regions where there is no image anymore:
# registered_img2_max = np.max(registered_img2)
# registered_img2[registered_img2 >= registered_img2_max] = 0
    


# # img1_blur = auto_blur_contrast(img1,thr_perc=0.4,blur=1)
# # img2_blur = auto_blur_contrast(registered_img2,thr_perc=0.4,blur=1)
# img1_blur = img1
# img2_blur = registered_img2
# ret1,img1_thr = cv2.threshold(img1_blur,0,1,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
# ret2,img2_thr = cv2.threshold(img2_blur,0,1,cv2.THRESH_BINARY+cv2.THRESH_OTSU)






