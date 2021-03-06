expModelFit = function(encodeModel, decodeModel, isFirstFit){
  # generate output directories
  dir.create("genData")
  dir.create("genData/simModelFit")
  dir.create("stanWarnings")
  dir.create(sprintf("genData/simModelFit/%s", encodeModel))
  dir.create(sprintf("genData/simModelFit/%s/%s", encodeModel, decodeModel))
  
  # load experiment parameters
  load("expParas.RData")
  
  # load sub-functions and packages
  library("dplyr"); library("tidyr")
  source("subFxs/loadFxs.R")
  source("subFxs/helpFxs.R")
  source('subFxs/modelFitGroup.R')
  source("expSchematics.R")
  
  # load simulation data 
  load(sprintf("genData/simulation/%s.RData", encodeModel))
  
  # prepare inputs
  outputDir = sprintf("genData/simModelFit/%s/%s", encodeModel, decodeModel)
  config = list(
    nChain = 4,
    nIter = 100,
    adapt_delta = 0.99,
    max_treedepth = 11,
    warningFile = sprintf("stanWarnings/sim_%s_%s.txt",encodeModel, decodeModel)
  )
  
  # if it is the first time to fit the model, fit all participants
  # otherwise, check model fitting results and refit those that fail any of the following criteria
  ## no divergent transitions 
  ## Rhat < 1.01 
  ## Effective Sample Size (ESS) > nChain * 100
  if(!isFirstFit){
    ids = names(trialData)
    paraNames = getParaNames(decodeModel)
    simPara = loadExpPara(paraNames, outputDir)
    passCheck = checkFit(paraNames, simPara)
    trialData = trialData[passCheck]
    
    # increase the num of Iterations 
    config = list(
      nChain = 4,
      nIter = 12000,
      adapt_delta = 0.99,
      max_treedepth = 11,
      warningFile = sprintf("stanWarnings/sim_refit_%s_%s.txt", encodeModel, decodeModel)
    )
  }
  
  # fit the model for all participants
  modelFitGroup(decodeModel, trialData, config, outputDir, F)
}

