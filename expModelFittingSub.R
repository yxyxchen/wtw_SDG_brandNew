# this script fits the RL model for each participant
# using Rstan
# here I change all my modelFitting function into the risk version
# while in stan, I have different expMofelfitting and modelFitting scripts for different things 
expModelFitting = function(modelName){
  # create outfiles
  nBlock = 3
  dir.create("genData")
  dir.create("genData/expModelFittingSub")
  dir.create(sprintf("genData/expModelFittingSub/%s", modelName))
  
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
  dir.create(sprintf("genData/expModelFittingSub/%s", modelName))
  model = stan_model(file = sprintf("stanModels/%s.stan", modelName))
  
  # load expData
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  idList = hdrData$ID[hdrData$stress == "no stress"]                 
  n = length(idList)                    
  
  # determine paras
  paras = getParas(modelName)
  if(paras == "wrong model name"){
    print(paras)
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
    thisTrialData = trialData[[thisID]]
    cond = unique(thisTrialData$condition)
    cIdx = ifelse(cond == "HP", 1, 2)
    excludedTrials = lapply(1 : nBlock, function(i)
      which(thisTrialData$trialStartTime > (blockSecs - tMaxs[cIdx]) &
              (thisTrialData$blockNum == i)))
    includeStart = which(thisTrialData$trialNum == 1)
    includeEnd = sapply(1 : nBlock, function(i){
      if(length(excludedTrials[[i]] > 0)){
        min(excludedTrials[[i]])-1
      }else{
        max(which(thisTrialData$blockNum ==i))  
      }
    })
    tempt = lapply(1 : nBlock, function(i)
      truncateTrials(thisTrialData, includeStart[i], includeEnd[i]))
    thisTrialData = do.call("rbind", tempt)
    # thisTrialData = block2session(thisTrialData) actually not really necessary, but
    fileName = sprintf("genData/expModelFittingSub/%s/s%d", modelName, thisID)
    modelFitting(thisTrialData, fileName, paras, model)
  }
}
