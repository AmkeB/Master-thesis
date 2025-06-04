# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 12:01:35 2025

@author: 106667
"""

import numpy as np
import matplotlib.pyplot as plt
import cv2
from scipy.ndimage import gaussian_filter, binary_dilation, binary_closing, binary_erosion, binary_opening
from copy import deepcopy
import gc
import math
from skimage.draw import circle_perimeter
import logging
from scipy.signal import fftconvolve
from Lars_image_analysis_functions import auto_blur_contrast

Nwell = 1
Nroi = 1
Nvar = 8
Npixels = 5120

LPb = 17.3 #value for 0.5
LPg = 20.7 #value for 1
LPr =  28.9 #value for 0.7
e465_460 = 0.38
e532_460 = 0.01
e465_532 = 0.01
e532_532 = 0.17
e640_640 = 1
e640_532 = 0.03
e532_640 = 0
e465_max = 83000
e532_max = 108000
e640_max = 102000 #not sure, kan niet vinden
e465_640 = 0
e640_460 = 0
n_465 = 0.93
n_532 = 0.87
n_640 = 0.11 #not sure, kan niet vinden

path = "\\\\store.erasmusmc.nl\\department\\gene\\chien_data\\Lab\\Labmembers\\Amke Bernaerts\\immune_ch_test"
# path = "\\\\store.erasmusmc.nl\\department\\gene\\chien_data\\Lab\\Labmembers\\Amke Bernaerts\\immune_ch_test\\ROI249"
# path = r"C:\Users\Amke\OneDrive\Documenten\MEP\thuiswerken\images"

file_blue = path+"\\ROI073_blue.tif"
file_green = path+"\\ROI074_green.tif"
file_red = path+"\\ROI072_red.tif"

# file_blue = path+"\\ROI249_blue image.tif"
# file_green = path+"\\ROI250_green image.tif"
# file_red = path+"\\ROI248_red image.tif"

img_blue = cv2.imread(file_blue,cv2.IMREAD_UNCHANGED)
img_green = cv2.imread(file_green,cv2.IMREAD_UNCHANGED)
img_red = cv2.imread(file_red,cv2.IMREAD_UNCHANGED)

# img_blue = deepcopy(image_blue)
# img_green = deepcopy(image_green)

img_blue_for_thresholding = auto_blur_contrast(img_blue,thr_perc=0.4,blur=1) #still change thr_perc to find best value!
img_green_for_thresholding = auto_blur_contrast(img_green,thr_perc=0.4,blur=1)
img_red_for_thresholding = auto_blur_contrast(img_red,thr_perc=0.4,blur=1)

img_blue_gaussian = gaussian_filter(img_blue, sigma=20)
img_green_gaussian = gaussian_filter(img_green, sigma=20)
img_red_gaussian = gaussian_filter(img_red, sigma=20)

BGb_mean = np.min(img_blue_for_thresholding)
BGg_mean = np.min(img_green_for_thresholding)
BGr_mean = np.min(img_red_for_thresholding)


ret_blue,th_blue = cv2.threshold(img_blue_for_thresholding,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
ret_green,th_green = cv2.threshold(img_green_for_thresholding,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
ret_red,th_red = cv2.threshold(img_red_for_thresholding,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)


#dilating background threshold:
bg_th_blue = np.copy(th_blue)
bg_th_blue[bg_th_blue==255]=1
bg_th_blue[bg_th_blue==0]=255
bg_th_blue[bg_th_blue==1]=0

kernel = np.ones((25,25))
dilated_bg_th_blue = cv2.erode(bg_th_blue,kernel)


Ib = np.zeros([Npixels,Npixels])
BGb = np.zeros([Npixels,Npixels])
Ig = np.zeros([Npixels,Npixels])
BGg = np.zeros([Npixels,Npixels])
Ir = np.zeros([Npixels,Npixels])
BGr = np.zeros([Npixels,Npixels])

N_signal_pixels_blue = 0
N_signal_pixels_green = 0
N_signal_pixels_red = 0

for i in range(Npixels):
    for j in range(Npixels):
        if th_blue[i,j] == 255:
            Ib[i,j] = img_blue[i,j]
            N_signal_pixels_blue+=1
        # elif th_blue[i,j] == 0:
        #     BGb[i,j] = img_blue[i,j]
            # enlarge -25 pixel for background??
        if dilated_bg_th_blue[i,j] == 255:
            BGb[i,j] = img_blue[i,j]

for i in range(Npixels):
    for j in range(Npixels):
        if th_green[i,j] == 255:
            Ig[i,j] = img_green[i,j]
            N_signal_pixels_green+=1
        elif th_green[i,j] == 0:
            BGg[i,j] = img_green[i,j]
            # enlarge -25 pixel for background??

for i in range(Npixels):
    for j in range(Npixels):
        if th_red[i,j] == 255:
            Ir[i,j] = img_red[i,j]
            N_signal_pixels_red+=1
        elif th_red[i,j] == 0:
            BGr[i,j] = img_red[i,j]
            # enlarge -25 pixel for background??

# BGb_mean = np.sum(BGb)/N_signal_pixels_blue
# BGg_mean = np.sum(BGg)/N_signal_pixels_green
# BGr_mean = np.sum(BGr)/N_signal_pixels_red
print('blue mean',BGb_mean)
print('green mean',BGg_mean)
print('red mean',BGr_mean)

G = (Ig - BGg_mean) / LPg
B = (Ib - BGb_mean) / LPb
R = (Ir - BGr_mean) / LPr

# N_532 = (G - e465_532 / e465_460 * B) / (e532_max * n_532 * (e532_532 - e465_532 / e465_460 * e532_460))
# N_465 = (B - e532_max * e532_460 * n_532 * N_532) / (e465_max * e465_460 * n_465)

# N_640 = (R - e532_max * e532_640 * n_532 * N_532) / (e640_max * e640_640 * n_640) ## AANPASSEN > 3 eqs 3 unknowns

V1b = e465_max*e465_460*n_465
V2b = e532_max*e532_460*n_532
V3b = e640_max*e640_460*n_640

V1g = e465_max*e465_532*n_465
V2g = e532_max*e532_532*n_532
V3g = e640_max*e640_532*n_640

V1r = e465_max*e465_640*n_465
V2r = e532_max*e532_640*n_532
V3r = e640_max*e640_640*n_640

N_460 = np.zeros([Npixels,Npixels])
N_532 = np.zeros([Npixels,Npixels])
N_640 = np.zeros([Npixels,Npixels])

# for i in range(Npixels):
#     for j in range(Npixels):
#         b = B[i,j]
#         g = G[i,j]
#         r = R[i,j]

#         a = np.array([[V1b,V2b,V3b],[V1g,V2g,V3g],[V1r,V2r,V3r]])
#         y = np.array([[b],[g],[r]])
#         N_all = np.linalg.solve(a,y)
        
#         N_460[i,j] = N_all[0][0]
#         N_532[i,j] = N_all[1][0]
#         N_640[i,j] = N_all[2][0]


# ## scale naar 0 tot 255 bijv en dan heb je je image
# # N_532 = N_all[1]

# print(np.max(N_532))
# green_end_img = (N_532/np.max(N_532)) * 255

# zelf oplossen:
K = (V3b*V1g-V3g*V1b)/(-V2b*V1g+V2g*V1b)
L = (G*V1b-B*V1g)/(-V2b*V1g+V2g*V1b)
rem = (V1r/V1b)*(B-V2b*L+V2r*L)

Nsyt = (R-rem)/(((V1r/V1b)+V2r)*K+V3r-(V1r*V3b)/V1b)
Ncd4 = (G*V1b-B*V1g-V3g*V1b*Nsyt+V3b*V1g*Nsyt)/(-V2b*V1g+V2g*V1b)
Npanck = (B-V2b*Ncd4-V3b*Nsyt)/V1b
print(np.max(Ncd4))
green_end_img_sol2 = (Ncd4/np.max(Ncd4)) * 255

end_img_zeros = np.copy(green_end_img_sol2)
end_img_zeros[end_img_zeros<0]=0
cv2.imwrite(path+"\\final_N_cd4_zeroes_lower_bg.tif",end_img_zeros)


# setting negative values to 0, because negative dye amount not possible:  >> niet doen, verlies informatie zo
# N_532[N_532<0] = 0
# N_465[N_465<0] = 0
# N_640[N_640<0] = 0

# I_blue_bt = e465_max * e465_532 * n_465 * N_465


# I_red_bt = e640_max * e640_532 * n_640 * N_640




# I_blue_bt_scaling_matrix = np.zeros([Npixels,Npixels])
# I_red_bt_scaling_matrix = np.zeros([Npixels,Npixels])
# for i in range(Npixels):
#     for j in range(Npixels):
#         if img_blue[i,j]!=0:
#             I_blue_bt_scaling_matrix[i,j] = I_blue_bt[i,j]/img_blue[i,j]
#         if img_red[i,j]!=0:
#             I_red_bt_scaling_matrix[i,j] = I_red_bt[i,j]/img_red[i,j]

# max_bt_value = np.max(I_blue_bt)
# I_blue_bt_scaling_matrix = I_blue_bt/max_bt_value

# I_red_bt = e640_max * e640_532 * n_640 * N_640
# max_bt_value_red = np.max(I_red_bt)
# I_red_bt_scaling_matrix = I_red_bt/max_bt_value_red



# print('final results')
# print(np.shape(N_532))
# print(np.max(N_532)) #gives 0.52
# print(N_532)
# print(np.shape(N_465))
# print(np.max(N_465)) #gives 0.12
# print(N_465)
# print('blue_bleedthrough',I_blue_bt)
# print(np.max(I_blue_bt)) #gives 94.82

#cv2.imwrite(filepath+os.sep+filename, img)

# cv2.imwrite(path+"\img_blue_for_thresholding.tif",img_blue_for_thresholding)
# cv2.imwrite(path+"\img_blue_thresholded.tif",th_blue)
# cv2.imwrite(path+"\img_green_for_thresholding.tif",img_green_for_thresholding)
# cv2.imwrite(path+"\img_green_thresholded.tif",th_green)
# cv2.imwrite(path+"\\blue_bt_for_scaling_2.tif",I_blue_bt_scaling_matrix)
# cv2.imwrite(path+"\\red_bt_for_scaling_2.tif",I_red_bt_scaling_matrix)
# cv2.imwrite(path+"\img_red_thresholded.tif",th_red)
# cv2.imwrite(path+"\\N_cd4.tif",green_end_img_sol2)













