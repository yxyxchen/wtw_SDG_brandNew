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

# load summaryData
nBlock = 3
nComb = 1
load("genData/expDataAnalysis/sessionData.RData")
load("genData/expDataAnalysis/kmOnGridSess.RData")
summaryData = sessionData
# load trialData since we need scheduledWait 
allData = loadAllData()
hdrData = allData$hdrData 
expTrialData = allData$trialData       
allIDs = hdrData$ID 


# re-simulate data
modelName = "PRbs"
dir.create(sprintf("figures/expModelRepitation/%s",modelName))
rep.PRbs = modelRepitation(modelName, summaryData, expTrialData, nComb)

modelName = "PRbsNC"
dir.create(sprintf("figures/expModelRepitation/%s",modelName))
rep.PRbsNC = modelRepitation(modelName, summaryData, expTrialData, nComb)

modelName = "Rlearn"
dir.create(sprintf("figures/expModelRepitation/%s",modelName))
rep.Rlearn = modelRepitation(modelName, summaryData, expTrialData, nComb)

modelName = "RlearnL"
dir.create(sprintf("figures/expModelRepitation/%s",modelName))
rep.RlearnL = modelRepitation(modelName, summaryData, expTrialData, nComb)
plotKMSC = F

modelName = "Rlearndb"
dir.create(sprintf("figures/expModelRepitation/%s",modelName))
rep.Rlearndb = modelRepitation(modelName, summaryData, expTrialData, nComb)

modelName = "reduce_gamma"
dir.create(sprintf("figures/expModelRepitation/%s",modelName))
rep.gamma = modelRepitation(modelName, summaryData, expTrialData, nComb)

modelName = "RlearnL"
expPara = loadExpPara(getParas(modelName), sprintf("genData/expModelFitting/%s", modelName))
useID = getUseID(expPara, getParas(modelName))
length(useID)
ids = hdrData$ID[hdrData$stress == "no stress"]
ids[!ids %in% useID]
expPara = loadExpParaExtra(getParas(modelName), sprintf("genData/expModelFitting/%s", modelName))
hist(expPara$phi97.5[!expPara$id %in% useID])

# to understand why Rlearn didn't converge 
expPara_PR = rep.PRbs$expPara
expPara_Rlearn = rep.RlearnL$expPara
useID_PR = getUseID(expPara_PR, getParas("PRbs"))
useID_Rlearn = getUseID(expPara_Rlearn, getParas("RlearnL"))

expPara_db = rep.Rlearndb$expPara
useID_db = getUseID(expPara_db, paras = getParas("Rlearndb"))

# unconverge and gamma
unConverge = expPara_PR$gamma[!expPara_PR$id %in% useID_Rlearn]
converge = expPara_PR$gamma[expPara_PR$id %in% useID_Rlearn]
a = hist(unConverge)
breaks = a$breaks
medians = a$mids
a = hist(unConverge, breaks = breaks)
b = hist(converge, breaks = breaks)
unConvergeRatio = a$counts / (a$counts + b$counts)
plot(medians, unConvergeRatio, xlab = "Gamma",
     ylab = "Unconverged ratio") 
hist(expPara_PR$gamma, breaks = breaks)

# unconverge and nExclude
unConverge = expPara_PR$nExclude[!expPara_PR$id %in% useID_Rlearn]
converge = expPara_PR$nExclude[expPara_PR$id %in% useID_Rlearn]
breaks = c(0, 5, 10, 15, 20, 25, 30, 35)
medians = c(2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5)
a = hist(unConverge, breaks = breaks)
b = hist(converge, breaks = breaks)
unConvergeRatio = a$counts / (a$counts + b$counts)
plot(medians, unConvergeRatio, xlab = "nExclude",
     ylab = "Unconverged ratio") 
hist(expPara_PR$nExclude)

ids = expPara_PR$id
ids[!ids %in% useID_Rlearn]

# initialize 
thisRep = rep.gamma
expPara = thisRep$expPara
repTrialData = thisRep$repTrialData
modelName = "reduce_gamma"
paras = getParas(modelName)

useID = getUseID(expPara, paras)
repNo = thisRep$repNo
nSub =(length(useID))
AUCRep_ = matrix(NA, nrow = nComb , ncol = nSub)
stdWdRep_ = matrix(NA, nrow = nComb, ncol = nSub)
kmOnGridRep_ = vector(mode = "list", length = nSub)
plotKMSC = F
for(sIdx in 1 : nSub){
  # prepare inputs
  id = useID[[sIdx]]
  cond = unique(summaryData$condition[summaryData$id == id])
  tMax = ifelse( cond == "HP", tMaxs[1], tMaxs[2])
  nTrial = summaryData$nTrial[summaryData$id == id]
  label = sprintf("sub%d", id)
  kmOnGridMatrix = matrix(NA, nrow = length(kmGrid), ncol = nComb)
  for(cIdx in 1 : nComb){
    thisRepTrialData = repTrialData[[repNo[cIdx, which(thisRep$useID == id)]]]
    kmscResults = kmsc(thisRepTrialData, min(tMaxs), label ,plotKMSC, kmGrid)
    AUCRep_[cIdx,sIdx] = kmscResults[['auc']]
    stdWdRep_[cIdx, sIdx] = kmscResults$stdWd
    kmOnGridMatrix[,cIdx] = kmscResults$kmOnGrid
  }
  kmOnGridRep_[[sIdx]] = kmOnGridMatrix
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

# I don't know
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
id = 1
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

with(thisRep,{
id = 1
sIdx = which(useID  == id)
tempt = repTrialData[[repNo[1,sIdx]]]
tempt$timeWaited =  matrix(unlist(lapply(1:nComb, function(i) repTrialData[[repNo[i,sIdx]]]$timeWaited)), ncol = nComb) %>%
  apply(MARGIN  = 1, FUN = mean) 
tempt = within(tempt, sapply(1 : length(timeWaited), function(i) ifelse(timeWaited[i] >= scheduledWait[i], tokenValue, 0)))
tempt$blockNum = junk$blockNum
trialPlots(tempt,"Model-generated Data")
ggsave(sprintf("figures/expModelRepitation/%s/sim_data__%d.png", modelName, id),
       width = 5, height = 4)
})


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
    kmOnGridRep = apply(kmOnGridRep_[[which(useID== thisID)]], MARGIN = 1, mean)
    junk = data.frame(time = kmGrid, exp = kmOnGrid, rep= kmOnGridRep)
    plotData = gather(junk, source, survival_rate, -time)
    p = ggplot(plotData, aes(time, survival_rate, color = source)) + geom_line(size = 2) + ggtitle(label) + displayTheme
    p = p + ylim(c(-0.1,1.1))
    print(p)
    readline("continue")
    fileName = sprintf("%s_s%d.png", modelName, thisID)
    #ggsave(fileName, width = 4, height =4)
  }
}
