# 
library('ggplot2')
library('plyr')
library('dplyr')
library('tidyr')
load("wtwSettings.RData")
source('subFxs/repetitionFxs.R') # called by simulate 
source("subFxs/helpFxs.R") # getParas
source("subFxs/plotThemes.R")
source("subFxs/loadFxs.R") # load scheduledWait from empirical data
source("subFxs/analysisFxs.R") 

# modelName 
modelName = "RL2"
repFun = getRepFun(modelName)

# load expData
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
ids = hdrData$ID[hdrData$stress == "no stress"]                
nSub = length(ids)      

# load expPara
paraNames = getParaNames(modelName)
parentDir ="genData/expModelFitting"; dirName = sprintf("%s/%sdb",parentDir, modelName)
expPara = loadExpPara(paraNames, dirName)

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
  thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excluedTrials,]
  cond = unique(thisTrialData$condition)
  scheduledWait = thisTrialData$scheduledWait
  id = factor(ids[sIdx], levels = ids)
  simTrialData[[id]] = repFun(paras, cond, scheduledWait)
}
hdrData$ID = factor(hdrData$ID, levels = hdrData$ID)
save(simTrialData, hdrData, file = "genData/simulation/simTrialData.RData")

