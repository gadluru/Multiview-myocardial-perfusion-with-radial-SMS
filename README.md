# Multiview-myocardial-perfusion-with-radial-SMS

This is the MATLAB package for reconstruction of radial simultaneous multi-slice (SMS) data. The sample datasets should be downloaded from
https://drive.google.com/drive/folders/1gKNN0XwiGs539vw7KlOcu3d3YRg5_Z6D?usp=sharing
and saved under /RawData/. We provide two datasets that are both dynamic contrast enhanced myocardial perfusion MRI, with one gated (*MID00121*) and one ungated (*MID00045*). There are two files for each dataset, one is the raw data with all the coils and another one with only 8 virtual coils (which is smaller in size). The 8 virtual coils were the first 8 principal components after performing principal component analysis (PCA) along the coil dimension of the raw data. Having either file in the /RawData/ folder is fine, however the reconstruction is performed with 8 virtual coils.

There are two *.m files named ‘Recon_gated_pixel_tracking.m’ and ‘Recon_ungated_pixel_tracking.m’ which are demonstrative codes for the gated and the ungated datasets. 

The codes were tested on MATLAB 2018b. The iteration part of the codes can be run on a MATLAB compatible GPU by setting ‘para.setting.ifGPU = 1’. If the GPU is not setting up correctly or does not have enough memory, this option may not working properly. 

The copyright of the codes is reserved by the cardiac MRI group of University of Utah. The codes are currently maintained by Ye Tian (ye.tian@utah.edu). 
