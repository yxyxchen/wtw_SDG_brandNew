# simulate virtual participants (with different parameters. multiple times for each)
load("wtwSettings.RData") 
source("subFxs/helpFxs.R") # getPars
source("subFxs/taskFxs.R") # drawSamples
source("subFxs/repetitionFxs.R") # getRepFunction
source("subFxs/simulationFxs.R")
modelName = "curiosityTrialSp" 
nSeq = 1
nRep = 10
nTrial = 50
# initialize simTrialData
simTrialData = vector(mode = "list", length = 2)
for(cIdx in 1 : 2){
  cond = conditions[cIdx]
  paraTable = data.frame(phi = c(0.02, 0.05, 0.08), tau = c(5, 10, 15),
                         gamma = c(0.85, 0.90, 0.95))
  if(cond == "HP") paraTable$zeroPoint = c(12, 16, 20) else paraTable$zeroPoint = c(25, 30, 35) 
  scheduledWaitList = replicate(nSeq, replicate(nTrial, drawSample(cond)),  simplify = F)
  simTrialData[[cIdx]] = simulate(modelName, nRep, paraTable, scheduledWaitList, cond)
}
dirName = sprintf("genData/simulation/%s", modelName)
dir.create(dirName)
save( simTrialData, file = sprintf("%s/%dTrial.RData", dirName, nTrial))

# simModelFitting
simModelFitting = function(modelName, paras, nTrial){
  # create outfiles
  dir.create("genData")
  dir.create("genData/simModelFitting")
  dir.create(sprintf("genData/simModelFitting/%s", modelName))
  dir.create(sprintf("genData/simModelFitting/%s/%dTrial", modelName, nTrial))
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
  
  # load simulationData
  load(sprintf("genData/simulation/%s/%dTrial.RData", modelName, nTrial))
  load(sprintf("genData/simulation/%s/%dTrialPara.RData", modelName, nTrial))
  
  # pararell computation settings
  nCore = parallel::detectCores() -1 # only for the local computer
  library("doMC")
  library("foreach")
  registerDoMC(nCore)
  
  # loop over participants 
  for(cIdx in 1:2){
    simTrialDataHere = simTrialData[[cIdx]]
    cond = conditions[cIdx]
    foreach(i = 1 : 1) %dopar% {
      lapply(1:nRep, function(rIdx){
        thisTrialData = simTrialDataHere[[simNo[rIdx, 1, i]]]
        timeWaited = thisTrialData$timeWaited
        scheduledWait = thisTrialData$scheduledWait
        trialEarnings = thisTrialData$trialEarnings
        timeWaited[trialEarnings > 0] = scheduledWait[trialEarnings > 0]
        fileName = sprintf("genData/simModelFitting/%s/%dTrial/%s_s%d_r%d", modelName, nTrial, cond, i, rIdx)
        modelFitting(cond, wIni, timeWaited, trialEarnings, scheduledWait, fileName, paras, model)
      })
    }
  }
}