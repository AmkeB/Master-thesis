# -*- coding: utf-8 -*-

import numpy as np
import time
import os
from glob import glob
import multiprocessing
from AB010_Lars_image_analysis_functions import load_images_and_find_shift_for_PARTS_v3 #, tophat, auto_blur_contrast, find_shift
import cv2
import matplotlib.pyplot as plt




if __name__ == '__main__':
    num_cores = multiprocessing.cpu_count()
    print("num_cores = {}".format(num_cores))

    filepath1 = r"\\store.erasmusmc.nl\department\gene\chien_data\Lab\Data_and_Analysis\Amke Bernaerts\AB009_insituFUNseq_washing_optimization\Raw data\AB009\step2_imaging_and_segmentation_20250319at133800\460nm_renamed"
    filepath2 = r"\\store.erasmusmc.nl\department\gene\chien_data\Lab\Data_and_Analysis\Amke Bernaerts\chapter4_thesis_image_analysis\image registration\AB009\460nm_parts_ROIs_no_boundaries"
    min_shift_necessary = 5
    # before running, check if there are any bubbles in the wells, those will make the program fail

    # in filepath1 is the original ROI tagging file (532nm z value) is present
    
    # ROIs_for_tagging = np.loadtxt(filepath1+"\\ROIs_for_tagging.txt", dtype=int)
    ROIs_for_tagging = np.loadtxt("\\\\store.erasmusmc.nl\\department\\gene\\chien_data\\Lab\\Data_and_Analysis\\Amke Bernaerts\\AB009_insituFUNseq_washing_optimization\\Raw data\\AB009\\step2_imaging_and_segmentation_20250319at133800\\ROIs_for_tagging.txt",dtype=int)


    # Creating folders    
    try:
        os.mkdir(filepath2+"\\shift plots v3")
        os.mkdir(filepath2+"\\corrected ROI files v3")
        os.mkdir(filepath2+"\\best_rotations_v3")
        os.mkdir(filepath2+"\\best_shifts_v3")
        os.mkdir(filepath2+"\\ufo2_shifted_imgs_v3")
    except:
        print("folder already exists")

    WL_files1 = glob(filepath1+"\\*.tif")
    WL_files1.sort()
    WL_files2 = glob(filepath2+"\\*.tif")
    WL_files2.sort()

    txt_done = glob(filepath2+"\\corrected ROI files v3\\tempROI*.txt")
    # shifts = np.zeros((len(WL_files1),2))

    begin = len(txt_done)
    temp_ending = len(WL_files2)
    final_ending = len(WL_files1)
    
    print("{} ROIs imaged, {} analyzed and {} going to be analyzed now.".format(temp_ending, begin, temp_ending-begin))

    m_dist = 20 # How much difference in pixels is possible between two alignments to be grouped together
    found_z_diff = 0#-250 #100 # manually found difference in z position (in unit of 0.1 um (10 = 1 um))
    
    T = time.perf_counter()
    to_do = np.arange(begin, temp_ending, 1)
    # to_do = [2,20,24,48,56]
    
    # here, load_images_and_find_shift_for_PARTS_v3 is run for each PARTS and UFO2 ROI using multiprocessing
    for i in range(int(len(to_do)/num_cores)+1): #np.arange(begin, temp_ending, num_cores):
        subprocess = []
        if num_cores*(i+1) < len(to_do):
            for j in to_do[i*num_cores:(i+1)*num_cores]:
                print(j)
                # t = time.perf_counter()
                ROIij = ROIs_for_tagging[j, :]
                proc = multiprocessing.Process(target=load_images_and_find_shift_for_PARTS_v3, kwargs={
                    "file_path": filepath2, "i": j, "file1": WL_files1[j], "file2": WL_files2[j], "coord": ROIij, "dz": found_z_diff, "rotations":[1,0.5,0,-0.5,-1], "include_rot":True, "endsize": final_ending, "max_dist": m_dist, "min_shift": min_shift_necessary,
                    "downsize_ufo2_img":True})
                # t2 = time.perf_counter()
                subprocess.append(proc)
                proc.start()
                # t3 = time.perf_counter()
                print("ID of process {}: {}".format(j, proc.pid))
                # print("it takes {} s to create multiproces and {} s to start it".format(t2-t, t3-t2))
        else:
            for j in to_do[i*num_cores:]:
                print(j)
                # t = time.perf_counter()
                ROIij = ROIs_for_tagging[j, :]
                proc = multiprocessing.Process(target=load_images_and_find_shift_for_PARTS_v3, kwargs={
                    "file_path": filepath2, "i": j, "file1": WL_files1[j], "file2": WL_files2[j], "coord": ROIij, "dz": found_z_diff, "rotations":[1,0.5,0,-0.5,-1], "include_rot":True, "endsize": final_ending, "max_dist": m_dist, "min_shift": min_shift_necessary,
                    "downsize_ufo2_img":True})
                # t2 = time.perf_counter()
                subprocess.append(proc)
                proc.start()
                # t3 = time.perf_counter()
                print("ID of process {}: {}".format(j, proc.pid))
                # print("it takes {} s to create multiproces and {} s to start it".format(t2-t, t3-t2))
        
        # print(len(subprocess))
        for proc in subprocess:
            proc.join()
    
    
    T2 = time.perf_counter()
    print("in total it took {} s".format(T2-T))
    
    if temp_ending == final_ending:
        ROI_final = np.zeros((final_ending,5), dtype=int)
        txt_done = glob(filepath2+"\\corrected ROI files v3\\tempROI*.txt")
        txt_done.sort()
        for ROI in txt_done:
            ROI_final += np.loadtxt(ROI, dtype=int)
            # os.remove(ROI)
            
        # np.savetxt(filepath2+"\\corrected ROI files\\ROIs_for_tagging_corrected.txt", ROI_final, fmt='%d', delimiter='\t')
            

