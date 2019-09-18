# 
library('dplyr')
library('tidyr')
source("subFxs/helpFxs.R") # getParas
source("subFxs/loadFxs.R") # load scheduledWait from empirical data

load("expParas.RData")

# modelName 
modelName = "QL2"
source(sprintf("subFxs/gnrModels/%s.R", modelName))
gnrModel = get(modelName)

# load data
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
ids = hdrData$id[hdrData$stress == "no_stress"]                 
nSub = length(ids)    

# load expPara
paraNames = getParaNames(modelName)
nPara = length(paraNames)
dirName = sprintf("%s/%s","genData/expModelFit", modelName)
expPara = loadExpPara(paraNames, dirName)
passCheck = checkFit(paraNames, expPara)

# simulation
set.seed(123)
simTrialData = list()
for(sIdx in 1 : nSub){
  if(passCheck[sIdx]){
    id = ids[sIdx]
    paras = as.double(expPara[sIdx, 1 : nPara])
    # prepare input
    thisTrialData = trialData[[id]] # here we id instead of sIdx
    # excluded some trials
    excluedTrials = which(thisTrialData$trialStartTime > (blockSec - max(tMaxs)) &
                            thisTrialData$condition == conditions[2])
    thisTrialData = thisTrialData[!(1 : length(thisTrialData$trialEarnings)) %in% excluedTrials,]
    simTrialData[[id]] = gnrModel(paras, thisTrialData$condition, thisTrialData$scheduledWait)
  }
}
trialData = simTrialData
save(trialData, file = sprintf("genData/simulation/%s.RData", modelName))

