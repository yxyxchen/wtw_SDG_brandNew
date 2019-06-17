# be careful to always to use id in code, instead of expTrialData
# determine if truncated
isTrun = T
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
nBlock = 3

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
modelName = "uniPrior"
dirName = sprintf("figures/expModelRepitation//%s", modelName)
dir.create(dirName)
paras = getParas(modelName)
parentDir = ifelse(dataType == "block", "genData/expModelFitting", "genData/expModelFittingSub")
dirName = sprintf("%s/%s",parentDir, modelName)
tempt = loadExpPara(paras, dirName)
useID = getUseID(tempt, paras)
expPara = merge(x=tempt,y=summaryData, by="id",all.x=TRUE)


# simulate nRep ztimes for each participants, with different parameter samples
repModelFun = getRepModelFun(modelName)
nComb = 10
nSub = length(useID)
repTrialData = vector(length = nSub * nComb, mode ='list')
repNo = matrix(1 : (nSub * nComb), nrow = nComb, ncol = nSub)
set.seed(231)
for(sIdx in 1 : nSub){
  id = useID[[sIdx]]
  # load para samples
  paraSamples = read.table(sprintf("%s/%s/s%d.txt", parentDir, modelName, id),sep = ",", row.names = NULL)
  # load behavioral inputs
  thisExpTrialData = expTrialData[[id]] # here we useID
  cond = unique(thisExpTrialData$condition)
  if(dataType == "block") thisExpTrialData = thisExpTrialData[thisExpTrialData$blockNum ==1, ]
  if(isTrun){
    excludedTrials = lapply(1 : nBlock, function(i)
      which(thisExpTrialData$trialStartTime > (blockSecs - tMaxs[i]) &
              (thisExpTrialData$blockNum == i)))
    includeStart = which(thisExpTrialData$trialNum == 1)
    includeEnd = sapply(1 : nBlock, function(i){
      if(length(excludedTrials[[i]] > 0)){
        min(excludedTrials[[i]])-1
      }else{
        max(which(thisExpTrialData$blockNum ==i))  
      }
    })
    tempt = lapply(1 : nBlock, function(i)
      truncateTrials(thisExpTrialData, includeStart[i], includeEnd[i]))
    thisExpTrialData = do.call("rbind", tempt)
  }
  scheduledWait = thisExpTrialData$scheduledWait
  # simulate
  for(cbIdx in 1 : nComb){
    paraSample = as.double(paraSamples[sample(1 : nrow(paraSamples), 1), 1 : length(paras)])
    tempt = repModelFun(paraSample, cond, scheduledWait)
    repTrialData[[repNo[cbIdx, sIdx]]] = tempt
  }
}

# calculate AUC and timeWaited from the reproduced data
plotKMSC = F
# initialize 
AUCRep_ = matrix(0, nComb, nSub)
stdWdRep_ = matrix(0, nComb, nSub)
# timeWaitedRep_ = vector(mode = "list", length = nSub)
kmOnGridRep_ = vector(mode = "list", length = nSub)
for(sIdx in 1 : nSub){
  # prepare inputs
  id = useID[[sIdx]]
  cond = unique(summaryData$condition[summaryData$id == id])
  tMax = ifelse( cond == "HP", tMaxs[1], tMaxs[2])
  nTrial = summaryData$nTrial[summaryData$id == id]
  label = sprintf("sub%d", id)
  # initialize output variables 
  kmOnGridRep = matrix(0, length(kmGrid), nComb)
  timeWaitedRep = matrix(0, nTrial, nComb)
  # calculate outputs for each participant 
  for(cbIdx in 1 : nComb){
    thisRepTrialData = repTrialData[[repNo[cbIdx, sIdx]]]
    # auc 
    kmscResults = kmsc(thisRepTrialData, min(tMaxs), label ,plotKMSC, kmGrid)
    AUCRep_[cbIdx, sIdx] = kmscResults[['auc']]
    stdWdRep_[cbIdx, sIdx] = kmscResults$stdWd
    kmOnGridRep[,cbIdx] = kmscResults$kmOnGrid
    # # timeWaited
    # timeWaitedRep[, cbIdx] = thisRepTrialData$timeWaited
  }
  # save outputs for each participant 
  kmOnGridRep_[[sIdx]] = kmOnGridRep
  # timeWaitedRep_[[sIdx]] = timeWaitedRep
}

# compare emipirical and reproduced AUC
muAUCRep = apply(AUCRep_, MARGIN = 2, mean);stdAUCRep = apply(AUCRep_, MARGIN = 2, sd)
minAUCRep = muAUCRep - stdAUCRep;maxAUCRep = muAUCRep + stdAUCRep
muStdWdRep = apply(stdWdRep_, MARGIN = 2, mean);stdStdWdRep = apply(stdWdRep_, MARGIN = 2, sd)
minStdWdRep = muStdWdRep - stdStdWdRep;maxStdWdRep = muStdWdRep + stdStdWdRep
data.frame(muAUCRep, minAUCRep, maxAUCRep,muStdWdRep, minStdWdRep, maxStdWdRep,
                      AUC = summaryData$AUC[summaryData$id %in% useID], stdWD = summaryData$stdWd[summaryData$id %in% useID],
                      condition = summaryData$condition[summaryData$id %in% useID]) %>%
  ggplot(aes(AUC, muAUCRep)) +  geom_errorbar(aes(ymin = minAUCRep, ymax = maxAUCRep), color = "grey") +
  geom_point(size = 2) + facet_grid(~condition) + 
  geom_abline(slope = 1, intercept = 0) + saveTheme + xlim(c(-2, 22)) + ylim(c(-2, 22)) +
  ylab("Model-generated (s)") + xlab("Observed (s)") + ggtitle("Average WTW") +
  myThemeBig + theme(plot.title = element_text(face = "bold", hjust = 0.5))
fileName = sprintf("figures/expModelRepitation/%s/AUC_AUCRep.png", modelName)
ggsave(filename = fileName,  width = 6, height = 4)

#
data.frame(muAUCRep, minAUCRep, maxAUCRep,muStdWdRep, minStdWdRep, maxStdWdRep,
           AUC = summaryData$AUC[summaryData$id %in% useID], stdWd = summaryData$stdWd[summaryData$id %in% useID],
           condition = summaryData$condition[summaryData$id %in% useID]) %>%
  ggplot(aes(stdWd, muStdWdRep)) + geom_point() + geom_errorbar(aes(ymin = minStdWdRep, ymax = maxStdWdRep), color = "grey") +
  geom_point(size = 2) + facet_grid(~condition) + 
  geom_abline(slope = 1, intercept = 0) + saveTheme  +
  ylab(expression(bold(paste("Model-generated (s"^2,")")))) +
  xlab(expression(bold(paste("Observed (s"^2,")")))) +ggtitle("std WTW") +
  myThemeBig + theme(plot.title = element_text(face = "bold", hjust = 0.5))
fileName = sprintf("figures/expModelRepitation/%s/std_stdRep.png", modelName)
ggsave(filename = fileName,  width = 6, height = 4)

# compare emipircal and reproduced trialPlot, for one participant 
id = 7
sIdx = which(useID  == id)
cond = unique(summaryData$condition[summaryData$id == id])
label = sprintf("Sub %d, %s", id, cond)
if(isTrun){
  junk = lastTrunc(expTrialData[[id]])
}
junk = block2session(junk)
trialPlots(junk, "Observed Data") 
ggsave(sprintf("figures/expModelRepitation/%s/actual_data_%d.png", modelName, id),
       width = 5, height = 4)

tempt = repTrialData[[repNo[1,sIdx]]]
tempt$timeWaited =  matrix(unlist(lapply(1:nComb, function(i) repTrialData[[repNo[i,sIdx]]]$timeWaited)), ncol = nComb) %>%
  apply(MARGIN  = 1, FUN = mean) 
tempt = within(tempt, sapply(1 : length(timeWaited), function(i) ifelse(timeWaited[i] >= scheduledWait[i], tokenValue, 0)))
tempt$blockNum = junk$blockNum
# tempt = repTrialData[[repNo[1,sIdx]]]
trialPlots(tempt,"Model-generated Data")
ggsave(sprintf("figures/expModelRepitation/%s/sim_data__%d.png", modelName, id),
       width = 5, height = 4)

# survival curve prediction
idList = hdrData$ID[hdrData$stress == "no stress"]
for(sIdx in 1 : nSub){
  thisID = idList[sIdx]
  if(thisID %in% useID){
    para = as.double(expPara[sIdx, 1 : length(paras)])
    label = sprintf('Subject %s, %s, %s, LL = %.1f',thisID, hdrData$condition[sIdx], hdrData$stress[sIdx], expPara$LL_all[sIdx])
    label = paste(label, paste(round(para, 3), collapse = "", seq = " "))
    # prepara data 
    cond = hdrData$condition[hdrData$ID == thisID]
    kmOnGrid = kmOnGrid_[[which(hdrData$ID == thisID)]]
    tMax = ifelse(cond == "HP", tMaxs[1], tMaxs[2])
    kmOnGridRep = apply(kmOnGridRep_[[which(idList== thisID)]], MARGIN = 1, mean)
    junk = data.frame(time = kmGrid, exp = kmOnGrid, rep= kmOnGridRep)
    plotData = gather(junk, source, survival_rate, -time)
    p = ggplot(plotData, aes(time, survival_rate, color = source)) + geom_line(size = 2) + ggtitle(label) + displayTheme
    print(p)
    readline("continue")
    fileName = sprintf("%s_s%d.png", modelName, thisID)
    #ggsave(fileName, width = 4, height =4)
  }
}
