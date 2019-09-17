# 
library('dplyr')
library('tidyr')
source("subFxs/helpFxs.R") # getParas
source("subFxs/loadFxs.R") # load scheduledWait from empirical data

load("expParas.RData")

# modelName 
modelName = "RL2"
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
parentDir ="genData/expModelFitting"; dirName = sprintf("%s/%s",parentDir, modelName)
expPara = loadExpPara(paraNames, dirName)


# check fit



# simulation
set.seed(123)
simTrialData = list()
for(sIdx in 1 : nSub){
  id = ids[sIdx]
  paras = as.double(expPara[expPara$id == id, 1 : length(paraNames)])
  # prepare input
  thisTrialData = trialData[[id]] # here we id instead of sIdx
  # excluded some trials
  excluedTrials1 = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[1]) &
                           thisTrialData$condition == conditions[1])
  excluedTrials2 = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[2]) &
                           thisTrialData$condition == conditions[2])
  excluedTrials = c(excluedTrials1, excluedTrials2)
  thisTrialData = thisTrialData[!(1 : length(thisTrialData$trialEarnings)) %in% excluedTrials,]
  cond = unique(thisTrialData$condition)
  scheduledWait = thisTrialData$scheduledWait
  id = ids[sIdx]
  simTrialData[[id]] = repFun(paras, cond, scheduledWait)
}
save(simTrialData, hdrData, file = sprintf("genData/simulation/%s.RData", modelName))

