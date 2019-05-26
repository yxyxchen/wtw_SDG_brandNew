# this script fits the RL model for each participant
# using Rstan
# here I change all my modelFitting function into the risk version
# while in stan, I have different expMofelfitting and modelFitting scripts for different things 
expModelFitting = function(modelName, paras){
  # create outfiles
  dir.create("genData")
  dir.create("genData/expModelFittingSub")
  dir.create(sprintf("genData/expModelFittingSub/%s", modelName))
  #  load libraries and set environments
  options(warn=-1, message =-1) # default settings borrowed somewhere
  library('plyr'); library(dplyr); library(ggplot2);library('tidyr');library('rstan') #load libraries
  Sys.setenv(USE_CXX14=1) # making rstan working on this device 
  rstan_options(auto_write = TRUE) # default settings borrowed somewhere
  library("loo")
  # source scripts
  source('subFxs/modelFittingFxs.R') # for fitting single case 
  source('subFxs/loadFxs.R') # for load data
  source("subFxs/helpFxs.R")
  load("wtwSettings.RData")
  library("coda") # calculate psr in modelFittingFxs
  # compile the stan model 
  dir.create(sprintf("genData/expModelFitting/%s", modelName))
  model = stan_model(file = sprintf("stanModels/%s.stan", modelName))
  
  # determine wIni
  QHPApOptim = 5 / 6 * stepDuration / (1 - 0.9) 
  QLPApOptim = 0.93 * stepDuration / (1 - 0.9) 
  if(any(paras == "phiR")){
    wIni = (5/6 + 0.93) / 2 * stepDuration
  }else if(any(paras == "gamma")){
    wIni = (QHPApOptim + QLPApOptim) / 2
  }else{
    print("wrong model name!")
    break
  }
  
  # load expData
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  # list with a named element for each subject ID.
  allIDs = hdrData$ID                   # column of subject IDs

  # load empirical data
  load("genData/expDataAnalysis/blockData.RData")
  idList = unique(blockData$id[blockData$stress == "no stress"])
  n = length(idList)
  
  # loop over participants 
  library("doMC")
  library("foreach")
  nCore = parallel::detectCores() -1 # needed in local computeres
  registerDoMC(nCore)
  
  foreach(i = 1 : n) %dopar% {
    thisID = idList[[i]]
    thisTrialData = trialData[[thisID]]
    timeWaited = thisTrialData$timeWaited
    scheduledWait = thisTrialData$scheduledWait
    trialEarnings = thisTrialData$trialEarnings
    timeWaited[trialEarnings > 0] = scheduledWait[trialEarnings > 0]
    cond = unique(thisTrialData$condition)
    # add for the risk version
    wtwEarly = blockData$wtwEarly[blockData$id == thisID & blockData$blockNum == 1]
    fileName = sprintf("genData/expModelFittingSub/%s/s%d", modelName, thisID)
    modelFitting(cond, wIni, wtwEarly, timeWaited, trialEarnings, scheduledWait, fileName, paras, model)
  }
}
