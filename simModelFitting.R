# simulate virtual participants (with different parameters. multiple times for each)
load("wtwSettings.RData") 
source("subFxs/helpFxs.R") # getPars
source("subFxs/taskFxs.R") # drawSamples
source("subFxs/repetitionFxs.R") # getRepFunction
source("subFxs/simulationFxs.R")

# load model names
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
ids = hdrData$ID[hdrData$stress == "no stress"]
nSub = length(ids)

modelName = "PRbs" 
repModelFun = getRepModelFun(modelName)
# initialize simTrialData
simTrialData = vector(mode = "list", length = nSub)
paraTable =  vector(mode = "list", length = nSub)
set.seed(123)
for(sIdx in 1 : nSub){
  thisID = ids[sIdx]
  thisTrialData = trialData[[thisID]]
  cond = unique(thisTrialData$condition)
  cIdx = ifelse(cond == "HP", 1, 2)
  excludedTrials = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[cIdx]))
  thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excludedTrials,]
  scheduledWait = thisTrialData$scheduledWait
  cond = unique(thisTrialData$cond)
  # read paras
  junk = read.csv(sprintf("genData/expModelFitting/%sdb/s%d_summary.txt", modelName, thisID), header = F)
  paras = junk[1 : (nrow(junk) - 1),1]
  
  simTrialData[[thisID]] = repModelFun(paras, cond, scheduledWait)
}
dirName = sprintf("genData/simulation/%s", modelName)
dir.create(dirName)
save( simTrialData, file = sprintf("%s/simTrialData.RData", dirName))

# simModelFitting
simModelFitting = function(modelName, paras, nTrial){
  # create outfiles
  dir.create("genData")
  dir.create("genData/simModelFitting")
  dir.create(sprintf("genData/simModelFitting/%s", modelName))
  dir.create(sprintf("genData/simModelFitting/%s/%dTrial", modelName, nTrial))
  
  #load libraries
  library('plyr'); library(dplyr); library(ggplot2);library('tidyr');library('rstan')
  library("loo")
  library("coda") 
  source('subFxs/modelFittingFxs.R') # for fitting each single participant
  source('subFxs/loadFxs.R') # for load data
  source("subFxs/helpFxs.R") # for getParas
  load("wtwSettings.RData")
  
  #  set the environment for Rstan
  options(warn=-1, message =-1) # run without this for one participant to chec everything
  Sys.setenv(USE_CXX14=1) # needed in local computeres
  rstan_options(auto_write = TRUE) 
  
  # compile the stan model 
  dir.create(sprintf("genData/expModelFittingSub/%s", modelName))
  model = stan_model(file = sprintf("stanModels/%s.stan", modelName))

  # load simulationData
  load(sprintf("genData/simulation/%s/%dTrial.RData", modelName, nTrial))
  load(sprintf("genData/simulation/%s/%dTrialPara.RData", modelName, nTrial))
  
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
  
  for(cIdx in 1 : 2){
    cond = conditions[cIdx]
    simTrialDataHere = simTrialData[[cIdx]]
    n = length(simTrialDataHere)
    foreach(i = 1 : n) %dopar% {
      thisTrialData = simTrialDataHere[[i]]
      sIdx = ceiling(i / nRep) # the index of virtual subject, each with a different para comb
      rIdx = i %% nRep # the index of repitition
      thisTrialData$condition = cond
      fileName = sprintf("genData/simModelFitting/%s/%dTrial/%s_s%d_r%d", modelName, nTrial, cond, sIdx, rIdx)
      modelFitting(thisTrialData, fileName, paras, model)
    }
  }


  # loop over participants 
  for(cIdx in 1:2){
    
    cond = conditions[cIdx]
    foreach(i = 1 : 1) %dopar% {
      lapply(1:nRep, function(rIdx){
        thisTrialData = simTrialDataHere[[simNo[rIdx, 1, i]]]
        timeWaited = thisTrialData$timeWaited
        scheduledWait = thisTrialData$scheduledWait
        trialEarnings = thisTrialData$trialEarnings
        timeWaited[trialEarnings > 0] = scheduledWait[trialEarnings > 0]
        
        modelFitting(cond, wIni, timeWaited, trialEarnings, scheduledWait, fileName, paras, model)
      })
    }
  }
}