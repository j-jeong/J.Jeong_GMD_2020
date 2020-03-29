# The code is written to provide an example of calculating tree-ring width from ORCHIDEE output. 
# This code is the old version corresponding to ORCHIDEE r5669
# Author: Jina Jeong (j.jeong@vu.nl)

rm(list=ls())
library(ncdf4)
library(dplR)
library(dplyr)

source('function/functions_J.R')
source('function/fun_sim_to_trw.R')
source('function/fun_cal_trw.R')


# site : target_site
# config : configuration
# fourd : if output is 4-dimentional or not. r5669 doesn't write 4-dim output.
# init : if consider initional diameter or not. 
site = 'brit019'
config = 'basic'

sim = cal_trw(site=site,config=config,fourd=F)

load('data/obs.RData')
rownames(sim) <- rownames(obs[[site]])
sim
