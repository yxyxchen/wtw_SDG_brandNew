# this script fits the RL model for each participant
# using Rstan
# here I change all my modelFitting function into the risk version
# while in stan, I have different simMofelfitting and modelFitting scripts for different things 
expModelFitting = function(modelName){
  # create outfiles
  nBlock = 3
  dir.create("genData")
  dir.create("genData/simModelFitting")
  dir.create(sprintf("genData/simModelFitting/%s", modelName))
  
  #load libraries
  library('plyr'); library(dplyr); library(ggplot2);library('tidyr');library('rstan')
  library("loo")
  library("coda") 
  source('subFxs/modelFittingFxs.R') # for fitting each single participant
  source('subFxs/loadFxs.R') # for load data
  source("subFxs/helpFxs.R") # for getParas
  load("wtwSettings.RData")
  source("subFxs/analysisFxs.R")
  
  #  set the environment for Rstan
  options(warn=-1, message =-1) # run without this for one participant to chec everything
  Sys.setenv(USE_CXX14=1) # needed in local computeres
  rstan_options(auto_write = TRUE) 
  
  # compile the stan model 
  model = stan_model(file = sprintf("stanModels/%s.stan", modelName))
  
  # load simData
  load("genData/simulation/simTrialData.RData")
  idList = factor(hdrData$ID[hdrData$stress == "no stress"], levels = hdrData$ID)
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

  foreach(i = 1 : n) %dopar% {
    thisID = idList[[i]]
    thisTrialData = simTrialData[[thisID]]
    cond = unique(thisTrialData$condition)
    fileName = sprintf("genData/simModelFitting/%s/s%s", modelName, thisID)
    modelFitting(thisTrialData, fileName, paraNames, model, modelName)
  }
}
