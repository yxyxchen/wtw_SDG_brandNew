# this script fits the RL model for each participant
# using Rstan
expModelFitting = function(modelName, paras){
  # create outfiles
  dir.create("genData")
  dir.create("genData/expModelFitting")
  dir.create(sprintf("genData/expModelFitting/%s", modelName))
 #  load libraries and set environments
  options(warn=-1, message =-1) # run without this for one participant to chec everything
  library('plyr'); library(dplyr); library(ggplot2);library('tidyr');library('rstan') #load libraries
  Sys.setenv(USE_CXX14=1) # needed in local computeres
  rstan_options(auto_write = TRUE) # default settings borrowed somewhere
  # options(mc.cores = parallel::detectCores()) # not encouranged, better to specify
  # cores for each function
  library("loo")
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
    print("wrong model name!")
    break
  }

  # load expData
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  # list with a named element for each subject ID.
  allIDs = hdrData$ID                   # column of subject IDs
  n = length(allIDs)                    # n
  
  # load empirical data
  load("genData/expDataAnalysis/blockData.RData")
  idList = unique(blockData$id)

  # loop over participants 
  nCore = parallel::detectCores() -1 # only for the local computer
  library("doMC")
  library("foreach")
  registerDoMC(nCore)
  
  foreach(i = 1 : 10) %dopar% {
    thisID = idList[[i]]
    thisTrialData = trialData[[thisID]]
    thisTrialData = thisTrialData[thisTrialData$blockNum == 1,]
    # delete the last trial, since the decision is interuptted
    thisTrialData = thisTrialData[1 : (nrow(thisTrialData) - 1),]
    timeWaited = thisTrialData$timeWaited
    scheduledWait = thisTrialData$scheduledWait
    trialEarnings = thisTrialData$trialEarnings
    timeWaited[trialEarnings > 0] = scheduledWait[trialEarnings > 0]
    cond = unique(thisTrialData$condition)
    fileName = sprintf("genData/expModelFitting/%s/s%d", modelName, thisID)
    modelFitting(cond, wIni, timeWaited, trialEarnings, scheduledWait, fileName, paras, model)
  }

  # for(i in 1 : n){
  #   thisID = idList[[i]]
  #   thisTrialData = trialData[[thisID]]
  #   thisTrialData = thisTrialData[thisTrialData$blockNum == 1,]
  #   timeWaited = thisTrialData$timeWaited
  #   scheduledWait = thisTrialData$scheduledWait
  #   trialEarnings = thisTrialData$trialEarnings
  #   timeWaited[trialEarnings > 0] = scheduledWait[trialEarnings > 0]
  #   cond = unique(thisTrialData$condition)
  #   wIni = ifelse(cond == "HP", wInisTheory[1], wInisTheory[2]) # wIni is the theoratical initial values
  #   # load paras in the full_model
  #   nE = length(getParas("full_model")) 
  #   expParaMedian = vector(length = nE)
  #   fileName = sprintf("genData/expModelFitting/%s/s%d.txt", "full_model", thisID)
  #   junk = read.csv(fileName, header = F)
  #   expParaMedian  = apply(junk[,1:nE], MARGIN = 2, median)
  #   
  #   # fitting
  #   fileName = sprintf("genData/expModelFitting/%s/s%d", modelName, thisID)
  #   modelFittingNo(cond, wIni, timeWaited, trialEarnings, scheduledWait, fileName, paras, model, expParaMedian)
  # } 
}
