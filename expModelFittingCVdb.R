# if there is an error, especially when the errors occurs after modelFitting
# nFits will not match summary.txt
# in expModelFittingCV.R, ids in useID, expPara, converges are 1,2,3,
# others, like all ids in the filename, are like original ids
expModelFitting = function(modelName){
  
  library('plyr'); library(dplyr); library(ggplot2);library('tidyr');library("stringr")
  source('subFxs/modelFittingFxs.R') # for fitting each single participant
  source('subFxs/loadFxs.R') # for load data
  source("subFxs/helpFxs.R") # for getParas
  load("wtwSettings.RData")
  
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
  
  # parameters
  nFold = 10
  
  # load expData
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  nSub = sum(hdrData$stress == "no stress")
  idsCV = 1 : (nSub * nFold) # id encoded in cvPara
  ids = hdrData$ID[hdrData$stress == "no stress"] # id encoded in trialData
  # initialize outputs
  
  # for a specific model 
  
  # detect the debug folder
  originalFile = sprintf("genData/expModelFittingCV/%s", modelName)
  dbFile = sprintf("genData/expModelFittingCV/%sdb", modelName)
  if(!file.exists(dbFile)){
    dir.create(dbFile)
    allFiles = list.files(path = originalFile)
    nFile = length(allFiles)
    if(nFile == length(idsCV)){
      lapply(1 : nFile, function(i) file.copy(sprintf("%s/%s", originalFile, allFiles[i]),
                                              sprintf("%s/%s", dbFile, allFiles[i])))
      print("creat the debug folder")
    }else{
      print("Wrong number of files in the original folder!")
      break
    }
  }
  
  # loop over models
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  
  # enter the refit procedure
  nLoop =1 
  while(nLoop < 15){
    # load cvPara
    cvPara = loadCVPara(paraNames, sprintf("genData/expModelFittingCV/%sdb", modelName),
                        "*_summary.txt")
    useID = getUseID(cvPara, paraNames)
    excID = idsCV[!idsCV %in% useID]
    # refit the mode
    if(length(excID) > 0){
      text = sprintf("Start to refit %d participants", length(excID))
      print(text)
      # compile the debug version of the model
      model = stan_model(file = sprintf("stanModels/%sdb.stan", modelName))
      foreach(i = 1 : length(excID)) %dopar% {
        # extract sIdx and fIdx from the id encoded in cvPara
        sIdx = ceiling(excID[i] / nFold)  # ceiling groups 1-10 together yet floor + 1 groups 0-9 together
        fIdx = excID[i] - (sIdx-1) * nFold
        text = sprintf("reFit s%d_f%d", ids[sIdx], fIdx)
        print(text)
        # update nFits and converge
        fitFile = sprintf("genData/expModelFittingCV/%sdb/afit_s%d_f%d.RData", modelName, ids[sIdx], fIdx)
        if(file.exists(fitFile)){
          load(fitFile)
          nFit = nFit  + 1
          save(nFit, file = fitFile)
        }else{
          nFit = 2
          save(nFit, file = fitFile)
        }
        # prepare data
        thisTrialData = trialData[[ids[sIdx]]]
        cond = unique(thisTrialData$condition)
        cIdx = ifelse(cond == "HP", 1, 2)
        excludedTrials = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[cIdx]))
        thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excludedTrials,]
        # select the training set
        load(sprintf("genData/expModelFittingCV/split/s%d.RData", ids[sIdx]))
        select = as.vector(partTable[-fIdx,])
        thisTrialData = thisTrialData[(1 : nrow(thisTrialData)) %in% select,]
        fileName = sprintf("genData/expModelFittingCV/%sdb/s%d_f%d", modelName,
                           ids[sIdx], fIdx)
        # refit
        # load upper and lower
        tempt = read.csv(sprintf("genData/expModelFittingCV/%sdb/s%d_f%d_summary.txt", modelName,
                                 ids[sIdx], fIdx),header = F)
        low= tempt[1:nPara,4]
        up = tempt[1 : nPara,8]
        converge = modelFittingCVdb(thisTrialData, fileName, paraNames, model, modelName, nPara, low, up)
      }
    }else{
      print("finished")
      print(modelName)
      print(nSub)
      break
    }# loop over participants
    nLoop = nLoop + 1
  } # end of the loops 
  # evaluate useID again
  cvPara = loadCVPara(paraNames, sprintf("genData/expModelFittingCV/%sdb", modelName),
                      "*_summary.txt")
  useID = getUseID(cvPara, paraNames)
  print("finished")
  print(modelName)
  print(length(useID))
}# end of the function




