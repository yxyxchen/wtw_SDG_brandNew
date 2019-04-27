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
source("subFxs/taskFxs.R")
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
paras = c(0.02, 15, 0.96)
modelName = "curiosityTrial"
repModelFun = getRepModelFun(modelName)
sIdx = 2
id = idList[[sIdx]]
cond = hdrData$cond[hdrData$ID == id]
thisExpTrialData = expTrialData[[id]]
thisExpTrialData = thisExpTrialData[thisExpTrialData$blockNum ==1, ]
scheduledWait = thisExpTrialData$scheduledWait
set.seed(123)
# cond = "LP"
# scheduledWait = unlist(lapply(1:1000, function(x) drawSample(cond)))
tempt = repModelFun(paras, cond, scheduledWait)
trialPlots(tempt, cond)

##### get a sense of the model from a lot scheduledWait######
# simluation for one para and a lot scheduledWait
paras = c(0.02, 10, 0.001) 
modelName = "curiosityTrialR"
repModelFun = getRepModelFun(modelName)
nRep = 1 # number of repetitions
trialData = vector(length = n * nRep, mode ='list')
repNo = matrix(1 : (n * nRep), nrow = n, ncol = nRep)
set.seed(231)
for(sIdx in 1 : n){
  id = idList[[sIdx]] 
  cond = hdrData$cond[hdrData$ID == id]
  thisExpTrialData = expTrialData[[id]]
  thisExpTrialData = thisExpTrialData[thisExpTrialData$blockNum ==1, ]
  scheduledWait = thisExpTrialData$scheduledWait
  for(rIdx in 1 : nRep){
    tempt = repModelFun(paras, cond, scheduledWait)
    trialData[[repNo[sIdx, rIdx]]] = tempt
  }
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


# one sequences example
paras = c(0.01, 15, 1)
set.seed(123)
modelName = "curiosityTrial"
repModelFun = getRepModelFun(modelName)
longScheduledWaitLP =  unlist(lapply(1:5000, function(x) drawSample("LP")))
tempt = repModelFun(paras, "LP", longScheduledWaitLP)
label = sprintf("tau = %.2f", paras[2])
trialPlots(tempt, label)
fileName = sprintf("zoomIn_%s.png", label)
trialPlots(truncateTrials(tempt, 4800, 5000), label)
ggsave(fileName, width = 6, height = 12)
actionValueViewer(truncateTrials(tempt, 4000, 4100))
plot(tempt$Qwaits[,4000])
which.min(tempt$Qwaits[,2000])

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



