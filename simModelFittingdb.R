# this script fits the RL model for each participant
# using Rstan
# here I change all my modelFitting function into the risk version
# while in stan, I have different simMofelfitting and modelFitting scripts for different things 
# current algorithm might, you know increase nFit yet didn't really refit
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
  
  # load simData
  load("genData/simulation/simTrialData.RData")
  idList = hdrData$ID
  n = length(idList)                
  
  originalFile = sprintf("genData/simModelFitting/%s", modelName)
  dbFile = sprintf("genData/simModelFitting/%sdb", modelName)
  if(!file.exists(dbFile)){
    dir.create(dbFile)
    allFiles = list.files(path = originalFile)
    nFile = length(allFiles)
    if(nFile == nSub * 3){
      lapply(1 : nFile, function(i) file.copy(sprintf("%s/%s", originalFile, allFiles[i]),
                                              sprintf("%s/%s", dbFile, allFiles[i])))
      print("creat the debug folder")
    }else{
      print("Wrong number of files in the original folder!")
      break
    }
  }
  
  # determine paras
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  if(paraNames == "wrong model name"){
    print(paraNames)
    break
  }
  
  # enter the refit process
  nLoop = 1
  while(nLoop < 15){
    # determine excID
    expPara = loadExpPara(paraNames,
                          sprintf("genData/simModelFitting/%sdb", modelName))
    useID = getUseID(expPara, paraNames)
    excID = ids[!ids %in% useID]
    
    # loop over excID
    n = length(excID)
    if(n > 0){
      text = sprintf("Start to refit %d participants", length(excID))
      print(text)
      foreach(i = 1 : n) %dopar% {
        thisID = excID[[i]]
        text = sprintf("refit s%d", thisID)
        print(text)
        # update nFits and converge
        fitFile = sprintf("genData/simModelFitting/%sdb/afit_s%d.RData", modelName, thisID)
        if(file.exists(fitFile)){
          load(fitFile); nFit = nFit  + 1; save(nFit, file = fitFile)
        }else{
          nFit = 2; save(nFit, file = fitFile)
        }
        
        # prepare
        thisTrialData = trialData[[thisID]]
        cond = unique(thisTrialData$condition)
        cIdx = ifelse(cond == "HP", 1, 2)
        excludedTrials = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[cIdx]))
        thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excludedTrials,]
        fileName = sprintf("genData/simModelFitting/%sdb/s%d", modelName, thisID)
        
        # load upper and lower
        tempt = read.csv(sprintf("genData/simModelFitting/%sdb/s%d_summary.txt", modelName, thisID),
                         header = F)
        low= tempt[1:nPara,4]
        up = tempt[1 : nPara,8]
        converge = modelFittingdb(thisTrialData, fileName, paraNames, model, modelName, nPara, low, up)
      }
      nLoop = nLoop + 1  
    }else{
      print("finished!")
      print(modelName)
      print(nSub)
      break
    }# loop over participants    
  }
    # evaluate useID again
    expPara = loadExpPara(paraNames,
                          sprintf("genData/simModelFitting/%sdb", modelName))
    useID = getUseID(expPara, paraNames)
    print("finished!")
    print(modelName)
    print(length(useID))
}
