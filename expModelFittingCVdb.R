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
    if(nFile == (nSub * nFold)){
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
    idsCV = cvPara$id
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
        id = excID[i]
        junk = str_locate(id, "s[0-9]+")
        sIdx = substr(id, (junk[1] + 1), junk[2]) # use to load fit.RData and trialData
        junk = str_locate(id, "f[0-9]+")
        fIdx =  as.double(substr(id, (junk[1] + 1), junk[2]))
        text = sprintf("reFit %s", id)
        print(text)
        # update nFits and converge
        fitFile = sprintf("genData/expModelFittingCV/%sdb/afit_%s.RData", modelName, id)
        if(file.exists(fitFile)){
          load(fitFile)
          nFit = nFit  + 1
          save(nFit, file = fitFile)
        }else{
          nFit = 2
          save(nFit, file = fitFile)
        }
        # prepare data
        thisTrialData = trialData[[sIdx]]
        cond = unique(thisTrialData$condition)
        cIdx = ifelse(cond == "HP", 1, 2)
        excludedTrials = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[cIdx]))
        thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excludedTrials,]
        # select the training set
        load(sprintf("genData/expModelFittingCV/split/s%s.RData", sIdx))
        select = c(1:5, as.vector(partTable[-fIdx,]))
        thisTrialData = thisTrialData[(1 : nrow(thisTrialData)) %in% select,]
        fileName = sprintf("genData/expModelFittingCV/%sdb/%s", modelName,
                           id)
        # refit
        # load upper and lower
        tempt = read.csv(sprintf("genData/expModelFittingCV/%sdb/%s_summary.txt", modelName, id),header = F)
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




