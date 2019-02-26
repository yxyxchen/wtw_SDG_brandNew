# library 
library('ggplot2')
library('dplyr')
library('tidyr')
load("wtwSettings.RData")
source('subFxs/simulationFxs.R') # 
source("subFxs/taskFxs.R")
source("subFxs/helpFxs.R")
source("subFxs/plotThemes.R")

# simulate
modelName = "full_model" 
nBlock = 1
nRep = 10
set.seed(123)
simulate(modelName, nBlock, nRep)



