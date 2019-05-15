 # be careful to always to use id in code, instead of expTrialData
library("ggplot2") 
library("dplyr")
library("tidyr")
source("subFxs/plotThemes.R")
load("wtwSettings.RData") 
source("subFxs/helpFxs.R") # getPars
source("subFxs/loadFxs.R") # load  expPara
source("subFxs/taskFxs.R") # drawSamples
source("subFxs/repetitionFxs.R") # getRepFunction
source("subFxs/analysisFxs.R") # kmsc, trialPlot

# input 
dataType = "sess"

# load blockData or seesData, since we need empirical AUC and nTrial. 
if(dataType == "block"){
  load("genData/expDataAnalysis/blockData.RData")
  load("genData/expDataAnalysis/kmOnGridBlock.RData")
  summaryData = blockData[blockData$blockNum == 1, ]
}else{
  load("genData/expDataAnalysis/sessionData.RData")
  load("genData/expDataAnalysis/kmOnGridSess.RData")
  summaryData = sessionData
}

# load trialData since we need scheduledWait 
allData = loadAllData()
hdrData = allData$hdrData           
expTrialData = allData$trialData       
allIDs = hdrData$ID 

# load expPara
modelName = "curiosityTrialSp"
paras = getParas(modelName)
parentDir = ifelse(dataType == "block", "genData/expModelFitting", "genData/expModelFittingSub")
dirName = sprintf("%s/%s",parentDir, modelName)
expPara = loadExpPara(paras, dirName)
useID = getUseID(expPara, paras)

# simulate nRep times for each participants, with different parameter samples
repModelFun = getRepModelFun(modelName)
nComb = 10
nSub = length(useID)
trialData = vector(length = nSub * nComb, mode ='list')
repNo = matrix(1 : (nSub * nComb), nrow = nComb, ncol = nSub)
set.seed(231)
for(sIdx in 1 : nSub){
  id = useID[[sIdx]]
  # load para samples
  paraSamples = read.table(sprintf("%s/%s/s%d.txt", parentDir, modelName, id),sep = ",", row.names = NULL)
  # load behavioral inputs
  thisExpTrialData = expTrialData[[id]]
  cond = unique(thisExpTrialData$condition)
  if(dataType == "block") thisExpTrialData = thisExpTrialData[thisExpTrialData$blockNum ==1, ]
  scheduledWait = thisExpTrialData$scheduledWait
  # simulate
  for(cbIdx in 1 : nComb){
    paraSample = as.double(paraSamples[sample(1 : nrow(paraSamples), 1), 1 : length(paras)])
    tempt = repModelFun(paraSample, cond, scheduledWait)
    trialData[[repNo[cbIdx, sIdx]]] = tempt
  }
}

# calculate AUC and timeWaited from the reproduced data
plotKMSC = F
# initialize 
AUCRep_ = matrix(0, nComb, nSub)
timeWaitedRep_ = vector(mode = "list", length = nSub)
kmOnGridRep_ = vector(mode = "list", length = nSub)
for(sIdx in 1 : nSub){
  # prepare inputs
  id = useID[[sIdx]]
  cond = summaryData$condition[summaryData$id == id]
  tMax = ifelse( cond == "HP", tMaxs[1], tMaxs[2])
  nTrial = summaryData$nTrial[summaryData$id == id]
  kmGrid = seq(0, tMax, by=0.1) 
  label = sprintf("sub%d", id)
  # initialize output variables 
  kmOnGridRep = matrix(0, length(kmGrid), nComb)
  timeWaitedRep = matrix(0, nTrial, nComb)
  # calculate outputs for each participant 
  for(cbIdx in 1 : nComb){
    thisTrialData = trialData[[repNo[cbIdx, sIdx]]]
    # auc 
    kmscResults = kmsc(thisTrialData, tMax, label ,plotKMSC, kmGrid)
    AUCRep_[cbIdx, sIdx] = kmscResults[['auc']]
    kmOnGridRep[,cbIdx] = kmscResults$kmOnGrid
    # timeWaited
    timeWaitedRep[, cbIdx] = thisTrialData$timeWaited
  }
  # save outputs for each participant 
  kmOnGridRep_[[sIdx]] = kmOnGridRep
  timeWaitedRep_[[sIdx]] = timeWaitedRep
}

# compare emipircal and reproduced trialPlot, quit different, try another way?
trialPlots(trialData[[repNo[1,2]]], "s")
trialPlots(block2session(expTrialData[[2]]), "s")

# compare emipirical and reproduced AUC
muAUCRep = apply(AUCRep_, MARGIN = 2, mean)
stdAUCRep = apply(AUCRep_, MARGIN = 2, sd)
minAUCRep = muAUCRep - stdAUCRep
maxAUCRep = muAUCRep + stdAUCRep
plotData = data.frame(muAUCRep, minAUCRep, maxAUCRep,
                      AUC = summaryData$AUC[summaryData$id %in% useID],
                      condition = summaryData$condition[summaryData$id %in% useID])
ggplot(plotData,
       aes(AUC, muAUCRep)) +  geom_errorbar(aes(ymin = minAUCRep, ymax = maxAUCRep), color = "grey") +
  geom_point() + facet_grid(~condition) + 
  geom_abline(slope = 1, intercept = 0) + saveTheme + xlim(c(-2, 45)) + ylim(c(-2, 45)) +
  ylab("Predicted / s") + xlab("Observed / s") + ggtitle("Predicted vs Observed AUC")
fileName = sprintf("figures/expModelRepitation/AUC_AUCRep_%s.pdf", modelName) 
ggsave(filename = fileName,  width = 6, height = 4)


# survival curve prediction
for(sIdx in 1 : n){
  thisID = idList[sIdx]
  if(thisID %in% useID){
    para = as.double(expPara[sIdx, 1 : length(paras)])
    label = sprintf('Subject %s, %s, %s, LL = %.1f',thisID, hdrData$condition[sIdx], hdrData$stress[sIdx], expPara$LL_all[sIdx])
    label = paste(label, paste(round(para, 3), collapse = "", seq = " "))
    # prepara data 
    cond = hdrData$condition[hdrData$ID == thisID]
    kmOnGrid = kmOnGrid_[[which(blockData$id == thisID & blockData$blockNum == 1)]]
    tMax = ifelse(cond == "HP", tMaxs[1], tMaxs[2])
    kmGrid = seq(0, tMax, by = 0.1)
    kmOnGridRep = kmOnGridRep_[[which(idList== thisID)]]
    junk = data.frame(time = kmGrid, exp = kmOnGrid, rep= kmOnGridRep)
    plotData = gather(junk, source, survival_rate, -time)
    p = ggplot(plotData, aes(time, survival_rate, color = source)) + geom_line(size = 2) + ggtitle(label) + displayTheme
    print(p)
    readline("continue")
    fileName = sprintf("%s_s%d.png", modelName, thisID)
    #ggsave(fileName, width = 4, height =4)
  }
}

######### simluation for sequences
expPara$condition = blockData$condition[blockData$blockNum == 1]
tempt = summarise(group_by(expPara, condition), phi = mean(phi), tau = mean(tau),
                  gamma = mean(gamma))
medianParaHP = as.double(tempt[1,])[-1]
medianParaLP = as.double(tempt[2,])[-1]

medianParaHP = c(0.05, 15)
medianParaLP = c(0.05, 30)
repModelFun = getRepModelFun("heuristicRL")
set.seed(231)
nRep = 10
trialData = vector(length = n * nRep, mode ='list')
repNo = matrix(1 : (n * nRep), nrow = n, ncol = nRep)
for(sIdx in 1 : n){
  id = idList[[sIdx]]
  cond = unique(blockData$condition[blockData$id == id])
  if(cond == "HP") para = medianParaHP else para =  medianParaLP
  thisExpTrialData = expTrialData[[id]]
  thisExpTrialData = thisExpTrialData[thisExpTrialData$blockNum ==1, ]
  scheduledWait = thisExpTrialData$scheduledWait
  for(rIdx in 1 : nRep){
   if(cond == "HP"){
     para = medianParaHP
   }else{
     para = medianParaLP
   }
    tempt = repModelFun(para, cond, scheduledWait)
    trialData[[repNo[sIdx, rIdx]]] = tempt
    # simDistMatrix[,rIdx] = abs(tempt$timeWaited - thisExpTrialData$timeWaited)
  }
  # simDist_[[sIdx]] = apply(simDistMatrix, 1, mean)
  # simDistSd_[[sIdx]] = apply(simDistMatrix, 1, sd)
}

# calculate AUC and timeWaited
plotKMSC = F
# initialize 
totalEarningsRep_ = matrix(0, n, nRep)
AUCRep_ = matrix(0, n, nRep)
timeWaitedRep_ = vector(mode = "list", length = n)
timeWaitedRepSd_ = vector(mode = "list", length = n)
kmOnGridRep_ = vector(mode = "list", length = n)
kmOnGridRepSd_ = vector(mode = "list", length = n)
for(sIdx in 1 : n){
  id = idList[[sIdx]]
  thisExpTrialData = expTrialData[[id]]
  thisExpTrialData = thisExpTrialData[thisExpTrialData$blockNum ==1, ] 
  tMax = ifelse( unique(thisExpTrialData$condition) == "HP", tMaxs[1], tMaxs[2])
  kmGrid = seq(0, tMax, by=0.1) 
  label = "asda"
  nTrial = nrow(thisExpTrialData)
  # initialize
  timeWaitedMatrix = matrix(0, nTrial, nRep)
  kmOnGridMatrix = matrix(0, length(kmGrid), nRep)
  for(rIdx in 1 : nRep){
    thisTrialData = trialData[[repNo[sIdx, rIdx]]]
    junk = thisTrialData$timeWaited
    timeWaitedMatrix[,rIdx] = junk
    
    totalEarningsRep_[sIdx, rIdx] =  sum(thisTrialData$trialEarnings)
    kmscResults = kmsc(thisTrialData,tMax,label,plotKMSC,kmGrid)
    AUCRep_[sIdx, rIdx] = kmscResults[['auc']]
    kmOnGridMatrix[,rIdx] = kmscResults$kmOnGrid
  }
  timeWaitedRep_[[sIdx]] = apply(timeWaitedMatrix, MARGIN = 1, mean)
  timeWaitedRepSd_[[sIdx]] = apply(timeWaitedMatrix, MARGIN = 1, sd)
  kmOnGridRep_[[sIdx]] = apply(kmOnGridMatrix, MARGIN = 1, mean)
  kmOnGridRepSd_[[sIdx]] = apply(kmOnGridMatrix, MARGIN = 1, sd)
}

# AUC prediction
plotData = blockData[blockData$blockNum == 1,]
plotData$AUCRep = apply(AUCRep_, MARGIN = 1, FUN = mean)
plotData$AUCRepSd = apply(AUCRep_, MARGIN = 1, FUN = sd)
plotData$AUCRepMin = plotData$AUCRep - plotData$AUCRepSd 
plotData$AUCRepMax = plotData$AUCRep + plotData$AUCRepSd 
ggplot(plotData[plotData$id %in% useID, ],
       aes(AUC, AUCRep)) +  geom_errorbar(aes(ymin = AUCRepMin, ymax = AUCRepMax), color = "grey")  + geom_point() + facet_grid(~condition) + 
  geom_abline(slope = 1, intercept = 0) + saveTheme + xlim(c(0, 40)) + ylim(c(0, 40))


which(plotData$AUCRep > 25)







