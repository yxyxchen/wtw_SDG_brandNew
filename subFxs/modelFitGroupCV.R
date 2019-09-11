# fit a model for multiple participants in Rstan 

# inputs:
# modelName : the model name
# trialData : a nSubx1 list, each element containing behavioral data for one participant
# config : a list containing the Rstab configuration 
# outputDir: the directory to save parameter estimations

# the config variable contains:
# nChain : number of chains 
# nIter : number of interations on each chain 
# adapt_delta: real number from 0 to 1, and increaseing it forces stan to explore the target distribution in a smaller step
# max_treedepth: maximal depth of the trees that stan evaluates during each iteration
# warningFile : file for saving warnings generated Rstan


modelFitGroupCV = function(modelName, trialData, config, outputDir){
  # create the output directory 
  dir.create(outputDir)
 
   # create the file for Rstan warnings and erros
  writeLines("", config[['warningFile']])
  
  #load libraries
  library('plyr'); 
  library('rstan');library("loo");library("coda") 
  library("doMC");library("foreach")
  source('subFxs/modelFitSingle.R') # for fitting each single participant
  load("expParas.RData")
  
  # compile the Rstan model 
  options(warn= 1) 
  Sys.setenv(USE_CXX14=1)
  rstan_options(auto_write = TRUE) 
  model = stan_model(file = sprintf("stanModels/%s.stan", modelName))
  
  # determine parameters 
  paraNames = getParaNames(modelName)
  
  # load expData
  ids = names(trialData)
  nSub = length(ids)                    
  
  # parallel compuation settings
  nCore = as.numeric(Sys.getenv("NSLOTS")) # needed for cluster
  if(is.na(nCore)) nCore = 1 # needed for cluster
  # nCore = parallel::detectCores() -1 
  # registerDoMC(nCore)
  foreach(i = 1 : nSub) %:% 
    foreach(j = 1 : 10) %dopar% {
        id = ids[[i]]
        thisTrialData = trialData[[id]]
        # truncate the last portion in each block 
        excludedTrials = which(thisTrialData$trialStartTime > (blockSec - max(tMaxs)))
        thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excludedTrials,]
        # Use 9 folds out of 10 as the training set. Noticeably, 
        # we always include the first 5 trials in the training set since 
        # they include essential info about prior beliefs. 
        load(sprintf("genData/expModelFitCV/kFold/s%s.RData", id))
        select = c(1:5, as.vector(trialAllocation[-j,]))
        thisTrialData = thisTrialData[(1 : nrow(thisTrialData)) %in% select,]
        # determine the output file
        outputFile = sprintf("%s/s%s_f%d", outputDir, id, j)
        # fit the model
        modelFitSingle(id, thisTrialData, modelName, paraNames, model, config, outputFile)
    }
}
