expModelCmpCV = function(){
  # libraries and scripts
  library("ggplot2")
  library("dplyr")
  library("tidyr")
  source("subFxs/helpFxs.R")
  source("subFxs/loadFxs.R")
  source("subFxs/plotThemes.R")
  load("expParas.RData")
  
  # load empirical data
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  ids = hdrData$id[hdrData$stress == "no_stress"]                 
  nSub = length(ids) 
  
  # check fit
  modelNames = c("QL1", "QL2", "RL1", "RL2")
  nModel = length(modelNames)
  passCheck_ = matrix(NA, nrow = nSub, ncol = nModel)
  for(i in 1 : nModel){
    modelName = modelNames[i]
    paraNames = getParaNames(modelName)
    cvPara = loadExpPara(paraNames, sprintf("genData/expModelFitCV/%s", modelName))
    passCheck_[,i] = checkFit(paraNames, cvPara)
  }
  
  # use parameters estimated from the traning set to predict choices in the test set
  # and the total logliklihood sumed over all test sets serves as the loss function
  LLCV_ = matrix(NA, nrow = nSub, ncol = nModel)
  for(m in 1 : nModel){
    # name of the current model
    modelName = modelNames[m]
    
    # select the loss function based on the current model
    source(sprintf("subFxs/lossFuns/%s.R", modelName))
    lossFun = get(modelName)
    
    # loop over participants
    for(sIdx in 1 : nSub){
      # load kFold and empirical data 
      id = ids[sIdx]
      load(sprintf("genData/expModelFitCV/kFold/s%s.RData", id))
      thisTrialData = trialData[[id]]
      # truncate trials at the end of the block
      excluedTrials = which(thisTrialData$trialStartTime > (blockSec - max(tMaxs)))
      thisTrialData = thisTrialData[! 1 : nrow(thisTrialData) %in%  excluedTrials,]
      # initialize loglikelihood obtained from 10-fold cross validation
      LLCVs = vector(length = 10)
      # loop over k folds
      for(fIdx in 1 : 10){
        cvPara = read.csv(sprintf("genData/expModelFitCV/%s/s%s_f%d_summary.txt", modelName, id, fIdx), header = F)
        paras = cvPara[1 : (nrow(cvPara) - 1), 1]
        lossResults = lossFun(paras, thisTrialData$condition,
                              thisTrialData$trialEarnings, thisTrialData$timeWaited)
        LLCVs[fIdx] = sum(lossResults$LLTrial_[1 : length(lossResults$LLTrial_) %in% trialAllocation[fIdx,]])
      }
      # save loglikelihood from cross validation 
      LLCV_[sIdx, m] = LLCV
    }
  }
  
}
