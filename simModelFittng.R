simModelFitting = function(modelName, paras){
  dir.create("genData")
  dir.create("genData/simModelFitting")
  dir.create(sprintf("genData/simModelFitting/%s", modelName))
  #  load libraries and set environments
  options(warn=-1, message =-1) # default settings borrowed somewhere
  library('plyr'); library(dplyr); library(ggplot2);library('tidyr');library('rstan') #load libraries
  Sys.setenv(USE_CXX14=1) # making rstan working on this device 
  rstan_options(auto_write = TRUE) # default settings borrowed somewhere
  options(mc.cores = parallel::detectCores())# enable multi-core precessors 
  library("loo")
  source("subFxs/helpFxs.R")
  # source scripts
  source('subFxs/modelFittingFxs.R') # for fitting single case 
  source('subFxs/loadFxs.R') # for load data
  load("wtwSettings.RData")
  library("coda") # calculate psr in modelFittingFxs
  # compile the stan model 
  model = stan_model(file = sprintf("stanModels/%s.stan", modelName))
  fileName = sprintf("genData/simulation/%s/simParas.RData", modelName)
  load(fileName)
  # load
  for(c in 1 : 2){
  cond = conditions[c]
  fileName = sprintf("genData/simulation/%s/trial%sData.RData", modelName, cond)
  load(fileName)
  if(cond == "HP") trialData = trialHPData else trialData = trialLPData
    for(i in 1 : nComb){
      for(rIdx in 1 : nRep){
        thisTrialData = trialData[[simNo[i, rIdx]]]
        timeWaited = thisTrialData$timeWaited
        scheduledWait = thisTrialData$scheduledWait
        trialEarnings = thisTrialData$trialEarnings
        timeWaited[trialEarnings > 0] = scheduledWait[trialEarnings > 0]
        wIni = ifelse(cond == "HP", wInisTheory[1], wInisTheory[2]) # wIni is the theoratical initial values 
        fileName = sprintf("genData/simModelFitting/%s/%s_s%d_r%d", modelName, cond, i, rIdx)
        modelFittingSimple(cond, wIni, timeWaited, trialEarnings, scheduledWait, fileName, paras, model)
      }
    }
  }
}
