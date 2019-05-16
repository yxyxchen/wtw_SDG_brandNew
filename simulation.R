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
idList = hdrData$ID
n = length(idList)

##### get a sense of the model ######
# simluation for one para and one scheduledWait from simulation or ...
# error prone.. try 3 and everything changes 
sIdx = 1
paras = c(0.03, 3.4, 0.95, 35)
modelName = "curiosityTrialSp"
repModelFun = getRepModelFun(modelName)
id = idList[[sIdx]]
cond = hdrData$cond[hdrData$ID == id]
thisExpTrialData = expTrialData[[id]]
thisExpTrialData = thisExpTrialData[thisExpTrialData$blockNum ==1, ]
scheduledWait = thisExpTrialData$scheduledWait
set.seed(123)
# cond = "LP"
# scheduledWait = unlist(lapply(1:800, function(x) drawSample(cond)))
tempt = repModelFun(paras, cond, scheduledWait)
trialPlots(tempt, cond)
fileName = sprintf("simulation_s%d.png", sIdx)
ggsave(fileName, width = 5, height = 4)

## test the effect of tau
# when tau is 0
tMax = tMaxs[2]
nStep = tMax / stepDuration
waitSteps = 1 : nStep
pWaits = rep(0.5, nStep)
wdResults = getWD(pWaits, stepDuration)

# assume my model quit at every one stepDuration. Otherwise, it will be larger
# like in exp
paras = c(0.03, 3.4, 0.95, 35)
modelName = "curiosityTrialSp"
repModelFun = getRepModelFun(modelName)
cond = "LP"
lenSeq = 200
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
for(i in 1 : nTau){
  tau = tauList[i]
  pWaits = sapply(1 : nStep, function(i) softMax(c(tempt$Qwaits[i,lenSeq],
                                                   tempt$Qquits[lenSeq]), tau))
  wdResults = getWD(pWaits, stepDuration)
  stdWds[i] = wdResults$std
  muWds[i] = wdResults$mu
}
plot(tauList, stdWds)
plot(tauList, muWds)

# will there be a lot of difference in Qwait as time goes by? I am not sure, let check!


# first of all, I would like to check how AUC chanes aross blocks 
load("genData/expDataAnalysis/blockData.RData")
load("genData/expDataAnalysis/kmOnGridSess.RData")
load("genData/expDataAnalysis/sessionData.RData")
ggplot(blockData, aes(blockNum, AUC, group = id)) + geom_line() +
  facet_grid(~condition)
group_by(blockData,condition, blockNum) %>% summarise(muAUC = mean(AUC))


summaryData = sessionData
muKms = vector(mode = "list", length = 2)
stdKms = vector(mode = "list", length = 2)
tGrid = seq(0, blockSecs, by = 0.1)
for(cIdx in 1 : 2){
  cond = conditions[cIdx]
  junk = kmOnGrid_[which(summaryData$condition == cond)]
  kmOnGrid = sapply(1 : length(junk), function(i) junk[[i]])
  muKms[[cIdx]] = apply(kmOnGrid, 1, mean)
  stdKms[[cIdx]] = apply(kmOnGrid, 1, sd)
}
plotData = data.frame(time = rep(tGrid,2), mu = unlist(muKms), std = unlist(stdKms),
                      condition = rep(conditions, each = length(tGrid)))

