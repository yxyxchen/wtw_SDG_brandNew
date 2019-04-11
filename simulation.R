# library 
library('ggplot2')
library('dplyr')
library('tidyr')
load("wtwSettings.RData")
source('subFxs/repetitionFxs.R') # 
source("subFxs/taskFxs.R")
source("subFxs/helpFxs.R")
source("subFxs/plotThemes.R")

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

# simulate a long scheduledWait sequence to find out the asympotatic results
parasList = list(c(0.05, 10, 0.99), c(0.05, 10, 0.90), c(0.05, 10, 0.5))
nParas = length(parasList)
nSeq = 5
set.seed(123)
longScheduledWaitLPList = lapply(1 : nSeq, function(i) unlist(lapply(1:1000, function(x) drawSample("LP"))))
nRep = 10

# determine model name here 
modelName = "curiosityTrial"
repModelFun = getRepModelFun(modelName)
for(i in 1 : nParas){
  paras = parasList[[i]]
  for(j in 1 : nSeq){
    longScheduledWaitLP = longScheduledWaitLPList[[j]]
      for(k in 1 : nRep){
        tempt = repModelFun(paras, "LP", longScheduledWaitLP)
      }
  }
}
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



