# this script fits the RL model for each participant
# using Rstan
# here I change all my modelFitting function into the risk version
# while in stan, I have different expMofelfitting and modelFitting scripts for different things 
expModelFitting = function(modelName){
  #load libraries
  library('plyr'); library(dplyr); library(ggplot2);library('tidyr');
  library("stringr")
  library("loo")
  library("coda") 
  source('subFxs/modelFittingFxs.R') # for fitting each single participant
  source('subFxs/loadFxs.R') # for load data
  source("subFxs/helpFxs.R") # for getParas
  load("wtwSettings.RData")
  source("subFxs/analysisFxs.R")
  
  #  set the environment for Rstan
  library('rstan')
  options(warn=-1, message =-1) # run without this for one participant to chec everything
  Sys.setenv(USE_CXX14=1) # needed in local computeres
  rstan_options(auto_write = TRUE) 
  
  # loop over participants 
  library("doMC")
  library("foreach")
  # nCore = as.numeric(Sys.getenv("NSLOTS")) # needed for cluster
  # if(is.na(nCore)) nCore = 1 # needed for cluster
  nCore = parallel::detectCores() -1 # only for the local computer
  registerDoMC(nCore)
  
  # compile the stan model 
  model = stan_model(file = sprintf("stanModels/%s.stan", paste(modelName, "db", sep = "")))
  
  # load expData
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  ids = hdrData$ID[hdrData$stress == "no stress"]                 
  nSub = length(ids)                    
  
  # load nFits
  fitFile = sprintf("genData/expModelFitting/%s/fit.RData", modelName)
  if(file.exists(fitFile)){
    load(fitFile)
  }else{
    nFits = rep(1, nSub)
  }
  
  # determine paras
  paras = getParas(modelName)
  nPara = length(paras)
  if(paras == "wrong model name"){
    print(paras)
    break
  }
  
  # determine excID
  expPara = loadExpPara(paras,
                      sprintf("genData/expModelFitting/%s", modelName))
  useID = getUseID(expPara, paras)
  excID = ids[!ids %in% useID]
  
  # loop over excID
  n = length(excID)
  if(n > 0){
    for(i in 1 : n) {
      thisID = excID[[i]]
      thisTrialData = trialData[[thisID]]
      cond = unique(thisTrialData$condition)
      cIdx = ifelse(cond == "HP", 1, 2)
      excludedTrials = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[cIdx]))
      thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excludedTrials,]
      fileName = sprintf("genData/expModelFitting/%s/s%d", modelName, thisID)
      # enter the refit procedure
      converge = F
      nRefit = 0
      while(nRefit <= 2 & !converge){
        # load upper and lower
        tempt = read.csv(sprintf("genData/expModelFitting/%s/s%d_summary.txt", modelName, thisID),
                         header = F)
        low= tempt[1:nPara,4]
        up = tempt[1 : nPara,8]
        converge = modelFittingdb(thisTrialData, fileName, paras, model, modelName, nPara, low, up) 
      } # exit the refit procedure
      nFits[which(ids == thisID)] = nFits + nRefit
    }# loop over participants
    save(nFits, file = sprintf("genData/expModelFitting/%s/fit.RData", modelName))
  }
}
