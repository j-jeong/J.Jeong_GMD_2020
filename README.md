# Code and data for the research paper: 'Using the International Tree-Ring Data Bank (ITRDB) records as century-long benchmarks for land-surface models'
### Author: Jina Jeong (j.jeong@vu.nl)

This repository contains codes, data, output, and configurations to reproduce the results and figures and results presented in the study. Detailed explanations for individual files are given below.

Note that this published repository doesn't contain the European biomass network dataset. The dataset is available after simple registration with your email address at: http://www.baci-h2020.eu/

After unzipping data in the directory, you can run benchmarks_verification.R.

#  This code is run with R 3.6.0. 

# Descriptions
## data
### ITRDB_obs.RData
List of 10 datasets selected from the ITRDB used for the first test, but for here, this file is used to run plot_Fig5_6_7_8.R.

## output
### BACI_sim.r5698.RData and BACI_sim.ccn.r5698.RData 
Lists of calculated TRW and number of individuals per size classes from 11-site simulations. This data is processed from output.nc with cal_simulated_trw.R
### BACI.DEO_r5698.nc
Example output nc file to run cal_simulated_trw.R

## R
### benchmarks_verification.R
Verification using modifiers to the simulated outputs resulting in the main result, Table 2. 
This code uses the European biomass network dataset. To run this, please download BACI datasets first.
### build_obs_Rdata.R
To build a list for convenience in further processes. To run this, please download BACI datasets first. This script should be run before using benchamarks_verification.R and plot_FigS2.R
### cal_simulated_trw.R
An example of calculating tree-ring width from ORCHIDEE output.
### functions_J.R
Set of small functions to process 4 tree-ring benchmarks.
### fun_cal_trw.R
Function to calculate tree-ring width from ORCHIDEE-output.
### plot_Fig5_6_7_8.R
Make explanatory figures for building benchmarks. This script reproduces Fig. 5 to 8.
### plot_FigS2.R
Make example benchmarking figures against the European biomass network dataset. This script reproduces Fig. S2 in the supplementary material.

## config
### xx.zip
set of each configuration for 11 sites used for running ORCHIDEE. Please note that site HD2 and TIC were located in the same pixel so spinup.zip has 10 configurations.
