# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 17:18:01 2025

@author: 106667
"""

from glob import glob
import numpy as np

path = r"\\store.erasmusmc.nl\department\gene\chien_data\Lab\Data_and_Analysis\Amke Bernaerts\AB009_insituFUNseq_washing_optimization\Analysis\PARTS_ROIs_no_boundaries\best_shifts_v3"

files = glob(path+'\ROI*.txt')
files.sort()

all_rots = np.zeros([len(files)])

for i in range(len(files)):
    text = np.loadtxt(files[i])
    all_rots[i] = text
    
np.savetxt(path+'\\all_rotations.txt',all_rots)



