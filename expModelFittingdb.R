# this script fits the RL model for each participant
# using Rstan
# here I change all my modelFitting function into the risk version
# while in stan, I have different expMofelfitting and modelFitting scripts for different things 
expModelFitting = function(modelName){
  # create outfiles
  nBlock = 3
  dir.create("genData")
  dir.create("genData/expModelFitting")
  dir.create(sprintf("genData/expModelFitting/%s", modelName))
  
  #load libraries
  library('plyr'); library(dplyr); library(ggplot2);library('tidyr');library('rstan')
  library("stringr")
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
  nPara = length(paras)
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
  
  # idList = c(1, 13, 25, 26, 32, 44, 51, 65, 77, 87, 91, 94, 101, 102, 109) # Rlearn
  idList = c(1, 15, 25, 26, 31, 44, 51, 53, 69, 71, 79, 91, 94, 96, 109, 110) # RlearnL
  # idList = c(2, 16, 26, 44, 56, 62, 63,64, 77, 85, 100, 109, 110) $ uniPiror
  # idList = c(20, 56, 65, 106, 110)
  # idList = c(20, 45) # PRbsNC
  n = length(idList)
  foreach(i = 1 : n) %dopar% {
    thisID = idList[[i]]
    thisTrialData = trialData[[thisID]]
    cond = unique(thisTrialData$condition)
    cIdx = ifelse(cond == "HP", 1, 2)
    excludedTrials = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[cIdx]))
    thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excludedTrials,]
    fileName = sprintf("genData/expModelFitting/%s/s%d", modelName, thisID)
    # load upper and lower
    tempt = read.csv(sprintf("genData/expModelFitting/%s/s%d_summary.txt", substr(modelName, 1, nchar(modelName) -2), thisID),
                     header = F)
    low= tempt[1:nPara,4]
    up = tempt[1 : nPara,8]
    modelFittingdb(thisTrialData, fileName, paras, model, modelName, nPara, low, up)
  }
}
