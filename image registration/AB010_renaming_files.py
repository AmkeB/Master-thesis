# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 09:38:56 2025

@author: 106667
"""

# renaming ROI files
import os
import re
from glob import glob

# Set the directory where the files are located
directory = "\\\\store.erasmusmc.nl\\department\\gene\\chien_data\\Lab\\Data_and_Analysis\\Amke Bernaerts\\AB009_insituFUNseq_washing_optimization\\Raw data\\AB009\\step2_imaging_and_segmentation_20250319at133800\\637nm_renamed"

# Get all tiff files in the directory
files = glob(directory+"\\*.tif")
files.sort()

# Iterate over the files and rename them
for i, file in enumerate(files):
    # Use regular expression to match and extract ROI number
    match = re.search(r'_ROI(\d+)', file)

    if match:
        # Get the original ROI number and adjust it
        original_roi_number = match.group(1)

        # Create a new ROI number (format it to have leading zeros)
        new_roi_number = f'{i:03d}'
        
        # Construct the new filename by replacing the old ROI number
        new_filename = file.replace(f'_ROI{original_roi_number}', f'_ROI{new_roi_number}')
        
        # Full file paths
        old_file_path = os.path.join(directory, file)
        new_file_path = os.path.join(directory, new_filename)
        
        # Rename the file
        os.rename(old_file_path, new_file_path)
        print(f'Renamed: {file} -> {new_filename}')