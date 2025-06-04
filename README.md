# Master-thesis
This Github contains the code used to perform image registration, cell segmentation and bleed-through subtractions as described in the master thesis.

Cell segmentation:
this code was developed by Li You and together with Amke Bernaerts adapted to segment tumor versus non-tumor cells.

Bleed-through subtractions:
this code was created by Amke Bernaerts and performs the calculations as shown in supplementary note 1

Image registration:
The code to create PARTS ROIs was created by Amke Bernaerts. The code to perform registration is based on code developed by Lars van Roemburg (AB010_Lars_image_analysis_functions.py) and I have added a new function to it that calculates the shift between PARTS and UFO2 ROIs using phase cross-correlation.

To obtain the results as presented in the thesis:
Data: 128 UFO2 ROIs (8 tissue sections) + series3 of PARTS image
-	running AB010_create_parts_ROIs.py with starting point [905,160] gives PARTS image for all channels in snake-shaped ROIs
-	running AB010_image_stitching_for_Parts.py gives only PARTS ROIs corresponding to the UFO2 ROIs without secureseal boundaries + rotates PARTS ROIs to be same orientation as UFO2
-	run parts_ufo2_image_registration_v3.py to get shifted ufo2 images that can be put into composite with parts_ROI_no_boundaries in fiji to assess shift
o	note: you can change ‘include_rot’ to True or False, depending on whether you want to find best rotation for each image or not
-	optional:
o	to shift another channel with the calculated shift: run combine shifts in one file v3.py to get shifts in all file and input this into shift_ufo2_nuclear_channel_v3.py
o	to get best_rotations in one file run: combine rotations in one file.py
