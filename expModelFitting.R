# this script fits the RL model for each participant
# using Rstan
# here I change all my modelFitting function into the risk version
# while in stan, I have different expMofelfitting and modelFitting scripts for different things 
expModelFitting = function(modelName){
  # output directories
  dir.create("genData")
  dir.create("genData/expModelFitting")
  dir.create(sprintf("genData/expModelFitting/%s", modelName))
  
  #load libraries
  library('plyr'); library(dplyr); library(ggplot2);library('tidyr');library('rstan')
  library("loo")
  library("coda") 
  source('subFxs/modelFittingFxs.R') # for fitting each single participant
  source('subFxs/loadFxs.R') # for load data
  source("subFxs/helpFxs.R") # for getParas
  load("expParas.RData")
  source("subFxs/analysisFxs.R")
  
  #  set the environment for Rstan
  options(warn= 1, message =-1) # run without this for one participant to chec everything
  Sys.setenv(USE_CXX14=1) # needed in local computeres
  rstan_options(auto_write = TRUE) 
  
  # compile the stan model 
  model = stan_model(file = sprintf("stanModels/%s.stan", modelName))
  
  # load expData
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  idList = hdrData$id[hdrData$stress == "no_stress"]                 
  n = length(idList)                    
  
  # determine paras
  paraNames = getParaNames(modelName)
  if(paraNames == "wrong model name"){
    print(paraNames)
    break
  }
  
  # loop over participants 
  library("doMC")
  library("foreach")
  # nCore = as.numeric(Sys.getenv("NSLOTS")) # needed for cluster
  # if(is.na(nCore)) nCore = 1 # needed for cluster
  nCore = parallel::detectCores() -1 # only for the local computer
  registerDoMC(nCore)
  
  dir.create("outputs")
  writeLines("", sprintf("outputs/%s_log.txt", modelName))
  foreach(i = 1 : 2) %dopar% {
      thisID = idList[[i]]
      thisTrialData = trialData[[thisID]]
      cond = unique(thisTrialData$condition)
      cIdx = ifelse(cond == "HP", 1, 2)
      excludedTrials = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[cIdx]))
      thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excludedTrials,]
      fileName = sprintf("genData/expModelFitting/%s/s%s", modelName, thisID)
      modelFitting(thisTrialData, fileName, paraNames, model, modelName)
  }
}
