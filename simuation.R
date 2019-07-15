# it is a playground script where I can:
# visually check whether my model cam reproduce the empirical data
# probe the effect of different parameters 
# understand how my model works
# library 
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
source("subFxs/taskFxs.R") # drawSample
source("subFxs/simulationFxs.R")
###### simulate using the empirical schedualedWait
# load raw data 
allData = loadAllData()
hdrData = allData$hdrData  
allIDs = hdrData$ID     
expTrialData = allData$trialData   
idList = hdrData$ID[hdrData$stress = "no stress"]
n = length(idList)
load("genData/expDataAnalysis/sessionData.RData")
##### get a sense of the model ######
# simluation for one para and one scheduledWait from simulation or ...
# error prone.. try 3 and everything changes 
sIdx = 2
paras =  c(0.02, 0.01, 10, 30, 0.02, 0.01)
modelName = "Rlearn"
repModelFun = getRepModelFun(modelName)
id = idList[[sIdx]]
cond = hdrData$cond[hdrData$ID == id]
thisExpTrialData = expTrialData[[id]]
# thisExpTrialData = thisExpTrialData[thisExpTrialData$blockNum ==1, ]
scheduledWait = thisExpTrialData$scheduledWait
set.seed(123)
tempt = repModelFun(paras, cond, scheduledWait)
trialPlots(tempt, cond)
fileName = sprintf("simulation_s%d.png", sIdx)
ggsave(fileName, width = 5, height = 4)

# for outliers 
id = 39
sIdx = which(useID  == id)

sIdx = 39
id = expPara$id[sIdx]
cond = unique(summaryData$condition[summaryData$id == id])
label = sprintf("Sub %d, %s", id, cond)
trialPlots(block2session(expTrialData[[id]]), label)
trialPlots(repTrialData[[repNo[1,sIdx]]],label)

## test the effect of tau
# when tau is 0
tMax = tMaxs[2]
nStep = tMax / stepDuration
waitSteps = 1 : nStep
pWaits = rep(0.5, nStep)
wdResults = getWD(pWaits, stepDuration)

# assume my model quit at every one stepDuration. Otherwise, it will be larger
paras = c(0.03, 3.5, 0.95, 35)
modelName = "curiosityTrialSp"
repModelFun = getRepModelFun(modelName)
cond = "LP"
tMax = tMaxs[2]
lenSeq = 100
set.seed(123)
scheduledWait = unlist(lapply(1:lenSeq, function(x) drawSample(cond)))
tempt = repModelFun(paras, cond, scheduledWait)
trialPlots(tempt, cond)

# when tau is positive inf
muWD = (which(tempt$Qwaits[,lenSeq] < tempt$Qquits[lenSeq])[1] - 1) * stepDuration
varWD = 0

# when tau = 10
tauList = 1:20
nTau = length(tauList)
stdWds = vector(length = length(tauList))
muWds = vector(length = length(tauList))
nStep = tMax / stepDuration
for(i in 1 : nTau){
  tau = tauList[i]
  pWaits = sapply(1 : nStep, function(i) softMax(c(tempt$Qwaits[i,lenSeq],
                                                   tempt$Qquits[lenSeq]), tau))
  wdResults = getWD(pWaits, stepDuration)
  stdWds[i] = wdResults$std
  muWds[i] = wdResults$mu
}
plot(tauList, stdWds, xlab = "Tau", ylab = "Std wtw")
plot(tauList, muWds, xlab = "Tau", ylab = "AUC")
plot(muWds, stdWds, xlab = "AUC", ylab = "Std wtw")


## the effect of paras in Qwaits and Qquits 
paras = c(0.03, 3.4, 0.95, 35)
modelName = "curiosityTrialSp"
repModelFun = getRepModelFun(modelName)
cond = "LP"
tMax = tMaxs[2]
lenSeq = 1000
set.seed(123)
scheduledWait = unlist(lapply(1:lenSeq, function(x) drawSample(cond)))
tempt = repModelFun(paras, cond, scheduledWait)
trialPlots(tempt, cond)