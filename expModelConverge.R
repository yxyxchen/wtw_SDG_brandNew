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

# loop over models
modelNames = c("PRbs", "PRbsNC", "Rlearn", "RlearnL", "reduce_gamma")

# select model
mIdx = 1
modelName = modelNames[mIdx]
paras = getParas(modelName)
nPara = length(paras)

# load nFits
fitFile = sprintf("genData/expModelFittingCV/%s/fit.RData", modelName)
if(file.exists(fitFile)){
  load(fitFile)
}else{
  nFits = matrix(rep(1, nrow * ncol), nrow = nFold, ncol = nSub)
}
  
# load cvPara
cvPara = loadCVPara(paras,
                  sprintf("genData/expModelFittingCV/%s", modelName),
                  "*.txt")
useID = getUseID(cvPara, paras)
excID = idsCV[!idsCV %in% useID]

# refit the mode
if(length(excID) > 0){
  # compile the debug version of the model
  model = stan_model(file = sprintf("stanModels/%s.stan", paste(modelName, "db", sep= "")))
  foreach(i = 1 : length(excID)) %dopar% {
    # extract sIdx and fIdx from the id encoded in cvPara
    sIdx = floor(excID[i] / nFold) + 1
    fIdx = excID[i] - (sIdx-1) * nFold
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
    fileName = sprintf("genData/expModelFittingCV/%s/s%d_f%d", modelName,
                       ids[sIdx], fIdx)
    # the refit loop
    nRefit = 0
    converge = F
    while(nRefit <= 2 & (!converge)){
      # load upper and lower
      tempt = read.csv(sprintf("genData/expModelFittingCV/%s/s%d_f%d_summary.txt", modelName,
                               ids[sIdx], fIdx),header = F)
      low= tempt[1:nPara,4]
      up = tempt[1 : nPara,8]
      converge = modelFittingdb(thisTrialData, fileName, paras, model, modelName, nPara, low, up)
      # update nRefit
      nRefit = nRefit + 1
    } # exit the refit procedure
    nFits[fIdx, sIdx] = nFits[fIdx, sIdx] + nRefit
  }# loop over participants
  save(nFits, sprintf("genData/expModelFittingCV/%s/fit.RData", modelName))
}

# loop over models



