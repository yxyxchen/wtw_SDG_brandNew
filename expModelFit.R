expModelFit = function(modelName, isFirstFit){
  # generate output directories
  dir.create("genData")
  dir.create("genData/expModelFit")
  dir.create("stanWarnings")
  
  # load experiment parameters
  load("expParas.RData")
  
  # load sub-functions and packages
  library("dplyr"); library("tidyr")
  source("subFxs/loadFxs.R")
  source("subFxs/helpFxs.R")
  source('subFxs/modelFitGroup.R')
  
  # prepare inputs
  allData = loadAllData()
  hdrData = allData$hdrData
  trialData = allData$trialData[hdrData$stress == "no_stress"]
  outputDir = sprintf("genData/expModelFit/%s", modelName)
  config = list(
    nChain = 4,
    nIter = 100,
    adapt_delta = 0.99,
    max_treedepth = 11,
    warningFile = sprintf("stanWarnings/exp_%s.txt", modelName)
  )
  
  # if it is the first time to fit the model, fit all participants
  # otherwise, check model fitting results and refit those that fail any of the following criteria
  ## no divergent transitions 
  ## Rhat < 1.01 
  ## Effective Sample Size (ESS) > nChain * 100
  if(!isFirstFit){
    ids = names(trialData)
    paraNames = getParaNames(modelName)
    expPara = loadExpPara(paraNames, outputDir)

    # detect participants with divergent transitions
    ## read all warnings line by line
    warns = read.table(config$warningFile, sep = "\n", stringsAsFactors = F)
    ## select warnings reporting divergent transitions 
    dt_warns = warns[sapply(1 : nrow(warns), function(i) grepl('*divergent*',
                                                      warns[i,1])),1]
    ## extract ids associated with divergent transitions
    dt_ids = sub("s", "", str_extract(dt_warns, "s[0-9]*"))
    
    # detect participants with high Rhats 
    RhatCols = which(str_detect(colnames(expPara), "hat"))[1 : length(paraNames)] # columns recording Rhats
    high_Rhat_ids = ids[apply(expPara[,RhatCols] >= 1.01, MARGIN = 1, sum) > 0]
      
    # detect participants with low ESSs
    ESSCols = which(str_detect(colnames(expPara), "Effe"))[1 : length(paraNames)]# columns recording ESSs
    low_ESS_ids = ids[apply(expPara[,ESSCols] < (config$nChain * 100), MARGIN = 1, sum) > 0]
    
    # identify participants satisifying all three criteria:
    fail_check_ids = unique(c(dt_ids, high_Rhat_ids, low_ESS_ids))
    trialData = trialData[ids %in% fail_check_ids]
    
    # increase the num of Iterations 
    config = list(
      nChain = 4,
      nIter = 12000,
      adapt_delta = 0.99,
      max_treedepth = 11,
      warningFile = sprintf("stanWarnings/exp_refit_%s.txt", modelName)
    )
  }
  
  # fit the model 
  modelFitGroup(modelName, trialData, config, outputDir)
}

