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
modelName = "baseline"
dirName = sprintf("figures/expModelRepitation//%s", modelName)
dir.create(dirName)
paras = getParas(modelName)
parentDir = ifelse(dataType == "block", "genData/expModelFitting", "genData/expModelFittingSub")
dirName = sprintf("%s/%s",parentDir, modelName)
tempt = loadExpPara(paras, dirName)
useID = getUseID(tempt, paras)
expPara = merge(x=tempt,y=summaryData, by="id",all.x=TRUE)


# simulate nRep times for each participants, with different parameter samples
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
timeWaitedRep_ = vector(mode = "list", length = nSub)
kmOnGridRep_ = vector(mode = "list", length = nSub)
for(sIdx in 1 : nSub){
  # prepare inputs
  id = useID[[sIdx]]
  cond = unique(summaryData$condition[summaryData$id == id])
  tMax = ifelse( cond == "HP", tMaxs[1], tMaxs[2])
  nTrial = summaryData$nTrial[summaryData$id == id]
  kmGrid = seq(0, tMax, by=0.1) 
  label = sprintf("sub%d", id)
  # initialize output variables 
  kmOnGridRep = matrix(0, length(kmGrid), nComb)
  timeWaitedRep = matrix(0, nTrial, nComb)
  # calculate outputs for each participant 
  for(cbIdx in 1 : nComb){
    thisRepTrialData = repTrialData[[repNo[cbIdx, sIdx]]]
    # auc 
    kmscResults = kmsc(thisRepTrialData, tMax, label ,plotKMSC, kmGrid)
    AUCRep_[cbIdx, sIdx] = kmscResults[['auc']]
    kmOnGridRep[,cbIdx] = kmscResults$kmOnGrid
    # timeWaited
    timeWaitedRep[, cbIdx] = thisRepTrialData$timeWaited
  }
  # save outputs for each participant 
  kmOnGridRep_[[sIdx]] = kmOnGridRep
  timeWaitedRep_[[sIdx]] = timeWaitedRep
}

# compare emipircal and reproduced trialPlot, for one participant 
id = 1
sIdx = which(useID  == id)
cond = unique(summaryData$condition[summaryData$id == id])
label = sprintf("Sub %d, %s", id, cond)
junk = block2session(expTrialData[[id]])
trialPlots(junk, "Observed Data") 
ggsave(sprintf("figures/expModelRepitation/%s/actual_data_%d.png", modelName, id),
       width = 5, height = 4)

tempt = repTrialData[[repNo[1,sIdx]]]
tempt$timeWaited =  matrix(unlist(lapply(1:nComb, function(i) repTrialData[[repNo[i,sIdx]]]$timeWaited)), ncol = nComb) %>%
  apply(MARGIN  = 1, FUN = mean) 
tempt = within(tempt, sapply(1 : length(timeWaited), function(i) ifelse(timeWaited[i] >= scheduledWait[i], tokenValue, 0)))
tempt$blockNum = junk$blockNum

# tempt = repTrialData[[repNo[1,sIdx]]]
trialPlots(tempt,"Model-predicted Data")
ggsave(sprintf("figures/expModelRepitation/%s/sim_data__%d.png", modelName, id),
       width = 5, height = 4)


# compare emipirical and reproduced AUC
muAUCRep = apply(AUCRep_, MARGIN = 2, mean)
stdAUCRep = apply(AUCRep_, MARGIN = 2, sd)
minAUCRep = muAUCRep - stdAUCRep
maxAUCRep = muAUCRep + stdAUCRep
data.frame(muAUCRep, minAUCRep, maxAUCRep,
                      AUC = summaryData$AUC[summaryData$id %in% useID], 
                      condition = summaryData$condition[summaryData$id %in% useID]) %>%
  ggplot(aes(AUC, muAUCRep)) +  geom_errorbar(aes(ymin = minAUCRep, ymax = maxAUCRep), color = "grey") +
  geom_point(size = 2) + facet_grid(~condition) + 
  geom_abline(slope = 1, intercept = 0) + saveTheme + xlim(c(-2, 45)) + ylim(c(-2, 45)) +
  ylab("Predicted (s)") + xlab("Observed (s)") + ggtitle("AUC") +
  myThemeBig + theme(plot.title = element_text(face = "bold", hjust = 0.5))
fileName = sprintf("figures/expModelRepitation/%s/AUC_AUCRep.png", modelName)
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
