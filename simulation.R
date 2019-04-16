# library 
library('ggplot2')
library('dplyr')
library('tidyr')
load("wtwSettings.RData")
source('subFxs/repetitionFxs.R') # 
source("subFxs/helpFxs.R")
source("subFxs/plotThemes.R")
source("subFxs/loadFxs.R") 
source("subFxs/analysisFxs.R")
###### simulate using the empirical schedualedWait
# load raw data 
allData = loadAllData()
hdrData = allData$hdrData  
allIDs = hdrData$ID     
expTrialData = allData$trialData   
idList = hdrData$ID
n = length(idList)

# simluation
paras = c(0.01, 10, 0.9) # remember to change the paras
set.seed(231)
modelName = "curiosityTrial"
repModelFun = getRepModelFun(modelName)
nRep = 1# number of repetitions
trialData = vector(length = n * nRep, mode ='list')
repNo = matrix(1 : (n * nRep), nrow = n, ncol = nRep)
for(sIdx in 1 : n){
  id = idList[[sIdx]] 
  cond = hdrData$cond[hdrData$ID == id]
  thisExpTrialData = expTrialData[[id]]
  thisExpTrialData = thisExpTrialData[thisExpTrialData$blockNum ==1, ]
  scheduledWait = thisExpTrialData$scheduledWait
  for(rIdx in 1 : nRep){
    tempt = repModelFun(paras, cond, scheduledWait)
    trialData[[repNo[sIdx, rIdx]]] = tempt
    # simDistMatrix[,rIdx] = abs(tempt$timeWaited - thisExpTrialData$timeWaited)
  }
  # simDist_[[sIdx]] = apply(simDistMatrix, 1, mean)
  # simDistSd_[[sIdx]] = apply(simDistMatrix, 1, sd)
}

for(i in 1:n){
  sIdx = repNo[i, 1]
  id = idList[sIdx]
  label = sprintf("%s,%s", hdrData$condition[hdrData$ID== i], id)
  thisTrialData = trialData[[sIdx]]
  trialPlots(thisTrialData, label)
  readline("continue")
  tMax = ifelse(label == "HP", 20, 40)
  kmGrid = seq(0, tMax, by=0.1) 
  thisTrialData = trialData[[sIdx]]
  kmscResults = kmsc(thisTrialData,tMax,label,T,kmGrid)
  readline("continue")
  # actionValueViewer(thisTrialData)
}

# simulate a long scheduledWait sequence to find out the asympotatic results
parasList = list(c(0.05, 10, 0.99), c(0.05, 10, 0.90), c(0.05, 10, 0.5))
nParas = length(parasList)
nSeq = 5
set.seed(123)
longScheduledWaitLPList = lapply(1 : nSeq, function(i) unlist(lapply(1:1000, function(x) drawSample("LP"))))
nRep = 10

# one sequences example
paras1 = c(0.02, 100, 0.8)
set.seed(123)
longScheduledWaitLP =  unlist(lapply(1:5000, function(x) drawSample("LP")))
tempt = repModelFun(paras1, "LP", longScheduledWaitLP)
label = sprintf("tau = %.2f", paras1[2])
trialPlots(tempt, label)
fileName = sprintf("zoomIn_%s.png", label)
trialPlots(truncateTrials(tempt, 4800, 5000), label)
ggsave(fileName, width = 6, height = 12)


# simulated scheduledWait
set.seed(123)
for(c in 1 : 2){
  cond  = conditions[c]
  if(cond == "HP"){
    scheduledWaitHPlist = lapply(1 : nRep, function(i) unlist(lapply(1:50, function(x) drawSample(cond))))
  }else{
    scheduledWaitLPlist = lapply(1 : nRep, function(i) unlist(lapply(1:50, function(x) drawSample(cond))))
  }
}
# define functions
simulate = function(modelName, nBlock, nRep, paraTable){
  dir.create("genData/simulation")
  dir.create(sprintf("genData/simulation/%s", modelName))
  # choose modelFun
  repModelFun = getRepModelFun(modelName)
  
  # determine paraComb
  paraComb = getParaComb(paraTable)
  nComb = nrow(paraTable) ^ ncol(paraTable)
  simNo = matrix(seq(1 : nComb * nRep), nrow = nComb, ncol = nRep)
  save("paraComb", "nComb", "nRep", "simNo", file = sprintf("genData/simulation/%s/simParas.RData", modelName))
  # initialize outputs
  trialData = vector(length = nComb * nRep, mode ='list')
  # loop over conditions
  for(condIdx in 1 : 2){
    cond = conditions[condIdx];
    if(cond == "HP") scheduledWaitList = scheduledWaitHPlist else scheduledWaitList = scheduledWaitLPlist
    # loop over repetions 
    for(h in 1 : nrow(paraComb)){
      para = as.double(paraComb[h,]);
      # calculate wIni
      for(j in 1 : nRep ){
        scheduledWait = scheduledWaitList[[j]]
        tempt =repModelFun(para, cond, scheduledWait)
        trialData[[simNo[h, j]]] = tempt
      }  
    }
    # save 
    if(cond == "HP"){
      trialHPData = trialData
      fileName = sprintf('genData/simulation/%s/trialHPData.RData', modelName)
      save(trialHPData,file = fileName)
    }else{
      trialLPData = trialData
      fileName =  sprintf('genData/simulation/%s/trialLPData.RData', modelName)
      save(trialLPData,file = fileName)
    }
  }
}


# simulate
modelName = "full_model" 
nBlock = 1
nRep = 10
paraTable = data.frame(phi = c(0.02, 0.05, 0.08), tau = c(5, 10, 15),
                       gamma = c(0.85, 0.90, 0.95), QwaitIni = c(2, 3, 4))
simulate(modelName, nBlock, nRep, paraTable)


modelName = "R_learning2" 
nBlock = 1
nRep = 10
paraTable = data.frame(phi1 = c(0.02, 0.05, 0.08), phi2 = c(0.02, 0.05, 0.08),tau = c(5, 10, 15))
simulate(modelName, nBlock, nRep, paraTable)



