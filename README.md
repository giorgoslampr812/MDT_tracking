Standard tracking algorithms developed for the MDT chambers at CRF in LMU. The algorithms utilize small changes in the RecoUtility.cxx of MiniDAQ ROOT_PLOT framework to read the DAT files data produced by muon cosmics. What follows is a brief explanation of the repository scripts and setup. 

Setup 
Any python virtual environment setup as python -m venv New_virtual_environment

RecoUtility.cxx
When stored in ATLAS_Online_Monitor/ROOT_plot/src/reco of https://github.com/romyers/ATLAS_Online_Monitor it is able to save the hit information in csv format for every 1 million hits. It then allows the employment of the other algorithms on the output hit files. 

Adc_hist.C 
Root executable that produces the adc time histogram for the fitting script in root file format. 

Drift_hist.C
Root executable that produces the drift time histogram for the fitting script in root file format. 

Fit_t0_and_tail.C 
Root executable that fits t0 and ttail to calculate tmax for the drift time calibration. 

rt_rel_mon.py 
Creates monotonic autocalibrated rt relation. 

