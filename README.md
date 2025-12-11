Standard tracking algorithms developed for the MDT chambers at CRF in LMU. The algorithms utilize small changes in the RecoUtility.cxx of MiniDAQ ROOT_PLOT framework to read the DAT files data produced by muon cosmics. What follows is a brief explanation of the repository scripts and setup. 

Setup -
Any python virtual environment setup as python -m venv New_virtual_environment then source New_virtual_environment/bin/activate. If any dependancies are missing then pip install uproot, pip install matplotlib, pip install pandas, pip install scipy and pip install Numba. Root executables can run with a root version or precompiled version, https://root.cern/install/. Cpp files need compilation described in file.

RecoUtility.cxx-
When stored in ATLAS_Online_Monitor/ROOT_plot/src/reco of https://github.com/romyers/ATLAS_Online_Monitor it is able to save the hit information in csv format for every 1 million hits. It then allows the employment of the other algorithms on the output hit files. 

Adc_hist.C -
Root executable that produces the adc time histogram for the fitting script in root file format. 

Drift_hist.C-
Root executable that produces the drift time histogram for the fitting script in root file format. 

Fit_t0_and_tail.C -
Root executable that fits t0 and ttail to calculate tmax for the drift time calibration. 

rt_rel_mon.py -
Creates monotonic autocalibrated rt relation which is saved as a root file. 

hit_radii.py - 
Assigns drift radii column in hits csv file using the rt relation.

muon_tracker_fixed.cpp - 
Finds perpendicular tracks with 6 hits using channel geometry for TDC pairs (mezzanine) using a seeding algorithm. 

pl_tr.py - 
Plots the first 50 or any unique track_id for debugging purposes. 

FILE FORMATS

hits.csv - Saved hit information from MiniDAQ DAT file, size 1 million hits. 

hits_with_radius.csv - Identical to hits.csv just containing drift_radius column. 

tracked_out.csv - Output of muon_tracker_fixed.cpp containing the 6 hits that form the track and some auxillary information about the track like id etc. 











