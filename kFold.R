# This function allocates empirial data into 10 folds. Specifically, 
# each fold gets one random trial within a sequence of 10 trials. Therefore,
# each fold contains a subset of trials evenly sampled from the learning process.

# We will conduct 10-fold cross validation afterwards based this allocation.

# Noticeably, we skip the first 5 trials. Since we want to always include the 
# first 5 trials in the training set.
 
kFold = function(){
  # num of folds
  nFold = 10
  
  # load libraries and sub-functions
  source('subFxs/loadFxs.R') 
  source("subFxs/helpFxs.R") 
  load("expParas.RData")
  dir.create("genData/expModelFitCV")
  dir.create("genData/expModelFitCV/kFold")
  
  # load empirical data
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  ids = hdrData$id[hdrData$stress == "no_stress"]
  nSub = length(ids)  
  
  # loop over participants
  set.seed(123)
  for(i in 1 : nSub){
    # load data for this participant
    id = ids[[i]]
    thisTrialData = trialData[[id]]
    # exclude trials at the end of the block
    excludedTrials = which(thisTrialData$trialStartTime > (blockSec - max(tMaxs)))
    thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excludedTrials,]
    # number of trials in one fold
    # since the total num of trials might not be divisble by 10
    # the exact size of each fold can be foldSize - 1
    foldSize = ceiling((nrow(thisTrialData) - 5) / nFold)
    # allocates trials into 10 folds
    # each row in the object trialAllocation contains trials for one fold
    trialAllocation = sapply(1 : foldSize, function(i) sample(1:nFold,replace = FALSE) + (i -1) * nFold + 5)
    fileName = sprintf("genData/expModelFitCV/kFold/s%s.RData",  id)
    save("trialAllocation", file = fileName)
  }
}