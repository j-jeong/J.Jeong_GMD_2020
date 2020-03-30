# Code and data for the research paper: 'Using the International Tree-Ring Data Bank (ITRDB) records as century-long benchmarks for land-surface models'
### Author: Jina Jeong (j.jeong@vu.nl)

This repository contains codes, data, output, and configurations to reproduce figures and results for the paper. Detailed explanations for files are below.

Note that this published repository doesn't contain the BACI dataset. BACI dataset is available after simple registration with email at the site: http://www.baci-h2020.eu/

After unzipping data in the directory, you can run verify_virture.R.

#  This code is run with R 3.6.0. 

# Descriptions
## data
### 10site_info.csv
Meta data for selected 10 sites
### obs.RData
List of 10 rwl datasets selected from ITRDB

## output
### simulation.XX.RData
Lists of 10-sites simulations with each configuration. This data is processed from output.nc with trw.output.R
### brit019_basic.nc
Example output nc file to run trw.output.R

## R
### process_benchmark.R
Building 4 benchmarks with 2 metrics using ITRDB datasets, and visualizing the result (Fig. 9).
### verify_virture.R
Verification of one assumption laid under 4 benchmarks, the overestimation through big-tree selection (Fig. 4).
This code uses the BACI dataset. To run this, please download datasets.
### trw.output
An example of calculating tree-ring width from ORCHIDEE output. 
### functions_J.R
Small functions to process 4 tree-ring benchmarks.
### fun_cal_trw.R
Function to calculate tree-ring width from ORCHIDEE-output.

## figure
### comb_benchmarks.svg
Figure proceeded from verfy_virture.R
### verify_benchmarks.png
Figure proceeded from process_benchmark.R

## config
### xx.zip
set of each configuration for 10 sites used for running ORCHIDEE. Please note that the name of configurations are different from the paper, which is specified in the name of zip files
