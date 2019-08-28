# this script fits the RL model for each participant
# using Rstan
# here I change all my modelFitting function into the risk version
# while in stan, I have different expMofelfitting and modelFitting scripts for different things 
expModelFitting = function(modelName){
  # output directories
  dir.create("genData")
  dir.create("genData/expModelFitting")
  dir.create(sprintf("genData/expModelFitting/%s", modelName))
  
  # generate the output log file
  dir.create("outputs")
  writeLines("", sprintf("outputs/%s_log.txt", modelName))
  
  #load libraries
  library('plyr'); library(dplyr); library(ggplot2);library('tidyr');
  library('rstan');library("loo");library("coda") 
  library("doMC");library("foreach")
  source('subFxs/modelFittingFxs.R') # for fitting each single participant
  source('subFxs/loadFxs.R') # for load data
  source("subFxs/helpFxs.R") # for getParas
  load("expParas.RData")
  source("subFxs/analysisFxs.R")
  
  # compile the Rstan model 
  options(warn= 0, message =-1) 
  Sys.setenv(USE_CXX14=1)
  rstan_options(auto_write = TRUE) 
  model = stan_model(file = sprintf("stanModels/%s.stan", modelName))
  
  # determine parameters 
  paraNames = getParaNames(modelName)
  
  # load expData
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  ids = hdrData$id[hdrData$stress == "no_stress"]                 
  nSub = length(ids)                    
  
  # parallel compuation settings
  nCore = as.numeric(Sys.getenv("NSLOTS")) # needed for cluster
  if(is.na(nCore)) nCore = 1 # needed for cluster
  # nCore = parallel::detectCores() -1 
  # registerDoMC(nCore)
  
  # loop over participants 
  foreach(i = 1 : nSub) %dopar% {
      id = ids[[i]]
      thisTrialData = trialData[[id]]
      # truncate the last portion in each block 
      excludedTrials = which(thisTrialData$trialStartTime > (blockSec - max(tMaxs)))
      thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excludedTrials,]
      fileName = sprintf("genData/expModelFitting/%s/s%s", modelName, id)
      modelFitting(thisTrialData, fileName, paraNames, model, modelName)
  }
}
