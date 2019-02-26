# be careful to always to use id in code, instead of expTrialData
library("ggplot2")
library(stringr)
source("subFxs/plotThemes.R")
load("genData/expDataAnalysis/blockData.RData")
source("subFxs/taskFxs.R") # used in repetition
source("subFxs/repetitionFxs.R")
source("subFxs/analysisFxs.R") # for analysis
load("wtwSettings.RData") # used in repetition
source("subFxs/loadFxs.R") # 
load("wtwSettings.RData")
load("genData/expDataAnalysis/kmOnGrid.RData")
source("subFxs/helpFxs.R")
# inputs
modelName = "full_model"
pars = getParas(modelName)
# load expPara
expPara = loadExpPara(modelName, pars)
#tempt= loadExpParaExtra(modelName, pars)
#expParaMode = tempt$expParaMode
#expParaMedian = tempt$expParaMedian
idList = unique(blockData$id)
n = length(idList)
RhatCols = which(str_detect(colnames(expPara), "hat"))[1 : length(pars)]
EffeCols = which(str_detect(colnames(expPara), "Effe"))[1 : length(pars)]
useID = idList[apply(expPara[,RhatCols] < 1.1, MARGIN = 1, sum) == length(pars) & 
                         apply(expPara[,EffeCols] >100, MARGIN = 1, sum) == length(pars)]

# load raw data 
allData = loadAllData()
hdrData = allData$hdrData  
allIDs = hdrData$ID     
expTrialData = allData$trialData       
n = length(idList)

# tempt
load("scheduledWait.RData")
# simluation 
set.seed(231)
repModelFun = getRepModelFun(modelName)
nRep = 10# number of repetitions
trialData = vector(length = n * nRep, mode ='list')
repNo = matrix(1 : (n * nRep), nrow = n, ncol = nRep)
# simDist_ =  vector(mode = "list", length = length(useID))
# simDistSd_ =  vector(mode = "list", length = length(useID))
for(sIdx in 1 : n){
  id = idList[[sIdx]]
  #para = as.double(expPara[sIdx, 1 : length(pars)])
  #para = as.double(expParaMedian[sIdx, 1 : length(pars)])
  paraList = read.table(sprintf("genData/expModelFitting/%s/s%d.txt", modelName, id),sep = ",", row.names = NULL)
  #cond = unique(blockData$condition[blockData$id == id])
  cond = "HP"
  #thisExpTrialData = expTrialData[[id]]
  #thisExpTrialData = thisExpTrialData[thisExpTrialData$blockNum ==1, ]
  #scheduledWait = thisExpTrialData$scheduledWait
  scheduledWait = scheduledWait.HP
  for(rIdx in 1 : nRep){
    para = as.double(paraList[sample(1 : nrow(paraList), 1), 1 : length(pars)])
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
  #thisExpTrialData = expTrialData[[id]]
  #thisExpTrialData = thisExpTrialData[thisExpTrialData$blockNum ==1, ] 
  #tMax = ifelse( unique(thisExpTrialData$condition) == "HP", tMaxs[1], tMaxs[2])
  tMax = 16
  kmGrid = seq(0, tMax, by=0.1) 
  label = "asda"
  #nTrial = nrow(thisExpTrialData)
  nTrial = 300
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
ggplot(plotData,aes(AUCRep)) + geom_histogram(bins = 6, fill = conditionColors[2]) + saveTheme + xlab("AUC Rep by Individual Paras") +
  saveTheme
ggsave("AUCRep_LP.pdf", width = 6, height = 4)
AUCSummary = plotData
ggplot(plotData[plotData$id %in% useID, ],
       aes(AUC, AUCRep)) +  geom_errorbar(aes(ymin = AUCRepMin, ymax = AUCRepMax), color = "grey")  + geom_point() + facet_grid(~condition) + 
  geom_abline(slope = 1, intercept = 0) + saveTheme + xlim(c(-2, 45)) + ylim(c(-2, 45))
fileName = sprintf("figures/expModelRepitation/AUC_AUCRep_%s.pdf", modelName) 
ggsave(filename = fileName,  width = 6, height = 4)

# plot the boxplot of AUCRep
plotData = data.frame(AUC = rep(blockData$AUC[blockData$blockNum == 1], nRep),
                      AUCRep = as.vector(AUCRep_), id = rep(expPara$id, nRep),
                      condition = rep(blockData$condition[blockData$blockNum == 1], nRep))

AUCRank = AUCSummary %>% group_by(condition) %>%  mutate(AUCRank = order(AUC, decreasing=TRUE))
plotData$AUCRank = rep(as.factor(AUCRank$AUCRank), each = nRep)
ggplot(plotData[plotData$id %in% useID, ], aes(AUC, AUCRep)) + geom_point() + facet_grid(~condition)

# survival curve prediction
for(sIdx in 20 : n){
  thisID = idList[sIdx]
  if(thisID %in% useID){
    para = as.double(expPara[sIdx, 1 : length(pars)])
    label = sprintf('Subject %s, %s, %s, LL = %.1f',thisID, hdrData$condition[sIdx], hdrData$stress[sIdx], expPara$LL_all[sIdx])
    label = paste(label, paste(round(para, 3), collapse = "", seq = " "))
    # prepara data 
    cond = subData$condition[[which(subData$id == thisID)]]
    kmOnGrid = kmOnGrid_[[which(subData$id == thisID)]]
    tMax = ifelse(cond == "HP", tMaxs[1], tMaxs[2])
    kmGrid = seq(0, tMax, by = 0.1)
    kmOnGrid = kmOnGrid_[[which(subData$id == thisID)]]
    kmOnGridRep = kmOnGridRep_[[which(idList== thisID)]]
    junk = data.frame(time = kmGrid, exp = kmOnGrid, rep= kmOnGridRep)
    plotData = gather(junk, source, survival_rate, -time)
    p = ggplot(plotData, aes(time, survival_rate, color = source)) + geom_line() + ggtitle(label) + displayTheme
    print(p)
    readline("continue")
  }
}

zoomInID = unique(blockData$id[blockData$AUC> 11 & blockData$AUC < 25 & blockData$condition == "LP"])
# trialData prediction
for(sIdx in 20 : n){
  thisID = idList[sIdx]
  #if(thisID %in% zoomInID){
    thisExpTrialData = expTrialData[[thisID]]
    thisExpTrialData = thisExpTrialData[thisExpTrialData$blockNum == 1,]
    nTrial = nrow(thisExpTrialData)
    if(thisID %in% useID){
      para = as.double(expPara[expPara$id == thisID, 1 : length(pars)])
      
      # prepara data 
      timeWaited = thisExpTrialData$timeWaited
      trialEarnings = thisExpTrialData$trialEarnings
      scheduledWait = thisExpTrialData$scheduledWait
      timeWaited[trialEarnings >0] = scheduledWait[trialEarnings >0]
      nAction = sum(round(ifelse(trialEarnings >0, ceiling(timeWaited / stepDuration), floor(timeWaited / stepDuration) + 1)))
      label = sprintf('Subject %s, %s, %s, -LL = %.2f',thisID, unique(blockData$condition[blockData$id == thisID]),
                      unique(blockData$stress[blockData$id == thisID]), -expPara$LL_all[expPara$id == thisID] / nAction)
      #label = paste(label, paste(round(para, 3), collapse = "", seq = " "))
      label = paste(label, sprintf("AUC = %.2f, AUCR = %.2f", AUCSummary$AUC[expPara$id == thisID], AUCSummary$AUCRep[expPara$id == thisID]))
      # prepare plotData
      plotData = data.frame(trialNum = rep(1 : nTrial, 2), timeWaited = c(timeWaited,
                                                                          timeWaitedRep_[[which(expPara$id == thisID)]]),
                            quitIdx = rep(trialEarnings == 0, 2), source = rep(c("exp", "rep"), each = nTrial))
      p = ggplot(plotData, aes(trialNum, timeWaited)) + geom_line(aes(color = source, alpha = source))  +
        geom_point(data = plotData[plotData$quitIdx == 1 & plotData$source == "exp", ], aes(trialNum, timeWaited)) + ggtitle(label)+
        displayTheme + scale_alpha_manual(values = c(0.8, 0.5))
      print(p)
      readline("continue")
    }
    
  #}
}

# simluation for sequences
expPara$condition = blockData$condition[blockData$blockNum == 1]
tempt = summarise(group_by(expPara, condition), phi = mean(phi), tau = mean(tau),
                  gamma = mean(gamma), QwaitIni = mean(QwaitIni))
medianParaHP = as.double(tempt[1,])[-1]
medianParaLP = as.double(tempt[2,])[-1]

set.seed(231)
paraHPList = 
paraLPPossible = data.frame(phi = )
paraLPist = getCombination(paraLPPossible)
nRep = 30# number of repetitions
trialData = vector(length = n * nRep, mode ='list')
repNo = matrix(1 : (n * nRep), nrow = n, ncol = nRep)
# simDist_ =  vector(mode = "list", length = length(useID))
# simDistSd_ =  vector(mode = "list", length = length(useID))
for(sIdx in 1 : n){
  id = idList[[sIdx]]
  cond = unique(blockData$condition[blockData$id == id])
  if(cond == "HP") para = medianParaHP else para =  medianParaLP
  thisExpTrialData = expTrialData[[id]]
  thisExpTrialData = thisExpTrialData[thisExpTrialData$blockNum ==1, ]
  scheduledWait = thisExpTrialData$scheduledWait
  for(rIdx in 1 : nRep){
   # if(cond == "HP"){
   #   para = paraHPList[rIdx, ]
   # }else{
   #   para = paraLPList[rIdx, ]
   # }
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

