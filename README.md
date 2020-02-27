# CSP MVAR ICA MATLAB

## Requirements
* EEGlab : https://sccn.ucsd.edu/eeglab/download.php
* add path to Luca's MVAR toolbox functions
* notBoxPlot : https://nl.mathworks.com/matlabcentral/fileexchange/26508-notboxplot

### This repository contains the required functions and subfolders for obtaining the results of information theory decomposition for epilepsy EEG/source activations signals: 

* ./codes/get_paper_results.m: which provides SCoT functionality and information theory measures calculation
* ./codes/funct_CSPVARICA_InfomaxICA_numComp.m : perform the CSPMVARICA approach, pretty much the same as SCoT does + automatic selection of components thanks to Riccardo
* ./codes/csp_numcomponents_rim_onlynum.m : CSP components number selection + plotting
* ./codes/plot_results_matlab_paper.m : boxplots of obtained results
* ./data/focal_data_filtered_for_SPS.mat : structure with all the data used in this study.
* ./data/x_baseline.mat, x_preictal.mat : sample trials, for fast testing, shape 19*625*300.
* ./results/mvar_source_results_MATLAB_paper_(COND1_COND2).mat : 3 structures with results for plot_results_matlab_paper.m
* ./functions_ : Luca's codes 
* ./Functions_MVAR : Luca's codes 
* ./eeglab
* ./notBoxPlot
