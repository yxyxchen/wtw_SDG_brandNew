 # be careful to always to use id in code, instead of expTrialData
library("ggplot2")
library("stringr")
library("dplyr")
library("tidyr")
source("subFxs/plotThemes.R")
source("subFxs/helpFxs.R")
source("subFxs/loadFxs.R") 

load("genData/expDataAnalysis/blockData.RData")
source("subFxs/taskFxs.R") # used in repetition
source("subFxs/repetitionFxs.R")
source("subFxs/analysisFxs.R") # for analysis
load("wtwSettings.RData") # used in repetition
load("wtwSettings.RData")
load("genData/expDataAnalysis/kmOnGridBlock.RData")


# load raw data 
allData = loadAllData()
hdrData = allData$hdrData  
allIDs = hdrData$ID     
expTrialData = allData$trialData   
idList = hdrData$ID
n = length(idList)

# inputs
modelName = "curiosityTrialR"
# paras = getParas(modelName)
paras = c("phi", "tau", "phiR")
# load expPara
expPara = loadExpPara(modelName, paras)
#tempt= loadExpParaExtra(modelName, pars)
#expParaMode = tempt$expParaMode
#expParaMedian = tempt$expParaMedian
idList = unique(blockData$id)
n = length(idList)
RhatCols = which(str_detect(colnames(expPara), "hat"))[1 : length(paras)]
EffeCols = which(str_detect(colnames(expPara), "Effe"))[1 : length(paras)]
useID = idList[apply(expPara[,RhatCols] < 1.1, MARGIN = 1, sum) == length(paras) & 
                         apply(expPara[,EffeCols] >100, MARGIN = 1, sum) == length(paras)]



# simluation using sample para
set.seed(231)
repModelFun = getRepModelFun(modelName)
nRep = 10# number of repetitions
trialData = vector(length = n * nRep, mode ='list')
repNo = matrix(1 : (n * nRep), nrow = n, ncol = nRep)
for(sIdx in 1 : n){
  id = idList[[sIdx]]
  #para = as.double(expPara[sIdx, 1 : length(paras)])
  paraList = read.table(sprintf("genData/expModelFitting/%s/s%d.txt", modelName, id),sep = ",", row.names = NULL)
  cond = unique(blockData$condition[blockData$id == id])
  thisExpTrialData = expTrialData[[id]]
  thisExpTrialData = thisExpTrialData[thisExpTrialData$blockNum ==1, ]
  scheduledWait = thisExpTrialData$scheduledWait
  for(rIdx in 1 : nRep){
    para = as.double(paraList[sample(1 : nrow(paraList), 1), 1 : length(paras)])
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
AUCSummary = plotData
ggplot(plotData[plotData$id %in% useID, ],
       aes(AUC, AUCRep)) +  geom_errorbar(aes(ymin = AUCRepMin, ymax = AUCRepMax), color = "grey")  + geom_point() + facet_grid(~condition) + 
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
    p = ggplot(plotData, aes(time, survival_rate, color = source)) + geom_line() + ggtitle(label) + displayTheme
    print(p)
    readline("continue")
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







