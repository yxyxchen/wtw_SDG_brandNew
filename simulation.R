# library 
library('ggplot2')
library('dplyr')
library('tidyr')
load("wtwSettings.RData")
source('subFxs/repetitionFxs.R') # 
source("subFxs/taskFxs.R")
source("subFxs/helpFxs.R")
source("subFxs/plotThemes.R")

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
    # loop over repetions 
    for(h in 1 : nrow(paraComb)){
      para = paraComb[h,];
      # calculate wIni
      for(j in 1 : nRep ){
        set.seed(j)
        scheduledWait = unlist(lapply(1:50, function(x) drawSample(cond)))
        
        para = c(0.1, 0.2, 1)
        tempt =repModelFun(para, cond, scheduledWait)
        trialPlots(tempt, "")
        
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


modelName = "R_learning" 
nBlock = 1
nRep = 10
paraTable = data.frame(phi1 = c(0.02, 0.05, 0.08), phi2 = c(0.02, 0.05, 0.08), tau = c(5, 10, 15))
simulate(modelName, nBlock, nRep, paraTable)



