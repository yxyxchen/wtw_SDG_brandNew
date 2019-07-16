# this script fits the RL model for each participant
# using Rstan
# here I change all my modelFitting function into the risk version
# while in stan, I have different expMofelfitting and modelFitting scripts for different things 
splitExpData = function(){
  source('subFxs/loadFxs.R') # for load data
  source("subFxs/helpFxs.R") # for getParas
  load("wtwSettings.RData")
  dir.create("genData/expModelFittingCV/split")
  # load expData
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  idList = hdrData$ID[hdrData$stress == "no stress"]                 
  n = length(idList)                    
  set.seed(123)
  for(i in 1 : n){
    thisID = idList[[i]]
    thisTrialData = trialData[[thisID]]
    cond = unique(thisTrialData$condition)
    cIdx = ifelse(cond == "HP", 1, 2)
    excludedTrials = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[cIdx]))
    thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excludedTrials,]
    # determine partitions 
    nPart = ceiling(nrow(thisTrialData) / nFold)
    partTable = sapply(1 : nPart, function(i) sample(1:nFold,replace = FALSE) + (i -1) * nFold)
    fileName = sprintf("genData/expModelFittingCV/split/s%d.RData",  thisID)
    save("partTable", file = fileName)
  }
}

expModelFitting = function(modelName){
  # model fitting parameters 
  nBlock = 3
  nFold = 10
  
  # create outfiles
  dir.create("genData")
  dir.create("genData/expModelFittingCV")
  dir.create(sprintf("genData/expModelFittingCV/%s", modelName))
  
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
  
  set.seed(123)
  foreach(i = 1 : n) %dopar%  {
    thisID = idList[[i]]
    thisTrialData = trialData[[thisID]]
    cond = unique(thisTrialData$condition)
    cIdx = ifelse(cond == "HP", 1, 2)
    excludedTrials = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[cIdx]))
    thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excludedTrials,]
    # determine partitions 
    load(sprintf("genData/expModelFittingCV/split/s%d.RData", thisID))
    
    # loop
    for(j in 1 : nFold){
      select = as.vector(partTable[-j,])
      thisTrialData = thisTrialData[(1 : nrow(thisTrialData)) %in% select,]
      fileName = sprintf("genData/expModelFittingCV/%s/s%d_f%d", modelName, thisID, j)
      modelFittingCV(thisTrialData, fileName, paras, model, modelName)
    }
  }
}
