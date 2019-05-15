# this script fits the RL model for each participant
# using Rstan
expModelFitting = function(modelName, paras){
  # create outfiles
  dir.create("genData")
  dir.create("genData/expModelFittingSb")
  dir.create(sprintf("genData/expModelFittingSb/%s", modelName))
  #  load libraries and set environments
  options(warn=-1, message =-1) # default settings borrowed somewhere
  library('plyr'); library(dplyr); library(ggplot2);library('tidyr');library('rstan') #load libraries
  Sys.setenv(USE_CXX14=1) # making rstan working on this device 
  rstan_options(auto_write = TRUE) # default settings borrowed somewhere
  options(mc.cores = parallel::detectCores())# enable multi-core precessors 
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
  if(modelName == "curiosityTrialRSp"){
    wIni = (5/6 + 0.93) / 2 * stepDuration
  }else if(modelName == "curiosityTrialSp"){
    wIni = (QHPApOptim + QLPApOptim) / 2
  }else{
    print("wrong model name!" )
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
  
  # loop over suvject
  for(i in 1 : n){
    thisID = idList[[i]]
    thisTrialData = trialData[[thisID]]
    thisTrialData = thisTrialData[-deleteTrials,]
    timeWaited = thisTrialData$timeWaited
    scheduledWait = thisTrialData$scheduledWait
    trialEarnings = thisTrialData$trialEarnings
    timeWaited[trialEarnings > 0] = scheduledWait[trialEarnings > 0]
    cond = unique(thisTrialData$condition)
    fileName = sprintf("genData/expModelFittingSub/%s/s%d", modelName, thisID)
    modelFitting(cond, wIni, timeWaited, trialEarnings, scheduledWait, fileName, paras, model)
  }
}
