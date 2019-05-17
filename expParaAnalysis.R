library("ggplot2")
library("dplyr")
library("Hmisc")
library("coin")
source("subFxs/plotThemes.R")
source("subFxs/loadFxs.R") # load blockData and expPara
source("subFxs/helpFxs.R") # getParas
source("subFxs/analysisFxs.R") # plotCorrelation and getCorrelation
load("wtwSettings.RData")

# input 
dataType = "sess"
modelName = "curiosityTrialSp"

# create output directories
dir.create("figures/expParaAnalysis")
saveDir = sprintf("figures/expParaAnalysis/%s", modelName)
dir.create(saveDir)

# load blockdata data
if(dataType == "block"){
  load("genData/expDataAnalysis/blockData.RData")
  load("genData/expDataAnalysis/kmOnGridBlock.RData")
  summaryData = blockData[blockData$blockNum == 1, ] # so summaryData only have something relevant
}else{
  load("genData/expDataAnalysis/sessionData.RData")
  load("genData/expDataAnalysis/kmOnGridSess.RData")
  summaryData = sessionData
}

# maybe I shoud truncate summaryData to make them the same as balabala 

# load expPara
paras = getParas(modelName)
nPara = length(paras)
parentDir = ifelse(dataType == "block", "genData/expModelFitting", "genData/expModelFittingSub")
dirName = sprintf("%s/%s",parentDir, modelName)
expPara = loadExpPara(paras, dirName)
useID = getUseID(expPara, paras)

# trait analysis # put it away
personality = read.csv("data/SDGdataset.csv")
traits = c("Delay.of.Gratification", "Barratt.Impulsiveness",
           "Intolerance.of.Uncertainty", "Trait.Anxiety..STAIT.")
traitNames = c("DelayGra", "Impulsive",
               "Uncertain", "Anxiety")
nTrait = length(traits)
traitParaCorr = vector(mode = "list", length = length(paras) * length(traits))
corrNo = matrix(1:(nPara * nTrait), nrow = nPara, ncol = nTrait)
for(pIdx in 1 : length(paras)){
  para = paras[pIdx]
  paraColor = paraColors[pIdx]
  for(trIdx in 1 : length(traits)){
    trait = traits[trIdx]
    traitName = traitNames[trIdx]
    # combine personlaity, expPara and condtion. useID only
    input = data.frame(personality[summaryData$id %in% useID,trait], expPara[expPara$id %in% useID,para],
                       summaryData$condition[summaryData$id %in% useID])
    traitParaCorr[[corrNo[pIdx, trIdx]]]= getCorrelation(input)
  }
}
dimNames = list(paras, traitNames)
rhoTable = lapply(1:2, function(j) matrix(sapply(1: (nTrait * nPara), function(i) traitParaCorr[[i]]$rhos[j]),
                                          nrow = nPara, dimnames = dimNames))
pTable = lapply(1:2, function(j) matrix(sapply(1: (nTrait * nPara), function(i) traitParaCorr[[i]]$ps[j]),
                                        nrow = nPara, dimnames = dimNames))
# plot trait analysis
library("corrplot")
for(i in 1 : 2){
  cond = conditions[i]
  fileName = sprintf("traitPara%s.png", cond)
  png(fileName)
  corrplot(rhoTable[[i]], 
           p.mat = pTable[[i]], 
           is.corr = T, 
           method = "color",
           insig = "label_sig",
           tl.col = "black", tl.srt = 15, tl.cex = 1.5) 
  mtext("Parameter", side = 2, line = 1.5, cex = 2)
  mtext("Trait", side = 3, line = 2, cex = 2, at = 2.5)
  dev.off()
}

# plot hist 
expPara$condition = summaryData$condition[summaryData$id %in% expPara$id]
for(i in 1 : length(paras)){
  para = paras[i]
  paraColor = paraColors[[i]]
  p = ggplot(expPara[expPara$id %in% useID,], aes_string(para)) + geom_histogram(bins = 6, fill = paraColor) + facet_grid(~condition)+
    saveTheme + ylab("Count") + xlab(capitalize(para))
  parentDir = ifelse(dataType == "sess", "figures/expParaAnalysisSub", "figures/expParaAnalysis")
  dir.create(parentDir)
  fileName = sprintf("%s/hist_%s.pdf", parentDir, para)
  ggsave(fileName, width = 6, height = 3)
}

# parse behavioral data
load("genData/expDataAnalysis/blockData.RData")
# I will construct a dataset of interest here
data = data.frame(tau = expPara$tau,
                  cvQuit = blockData$cvQuitTime[(blockData$blockNum == 3) & (blockData$id %in% expPara$id)],
                  muQuit =  blockData$muQuitTime[(blockData$blockNum == 3) & (blockData$id %in% expPara$id)],
                  stdQuit =  blockData$stdQuitTime[(blockData$blockNum == 3) & (blockData$id %in% expPara$id)],
                  nQuit =  blockData$nQuit[(blockData$blockNum == 3) & (blockData$id %in% expPara$id)],
                  stdWd = blockData$stdWd[(blockData$blockNum == 3) & (blockData$id %in% expPara$id)],
                  cvWd = blockData$cvWd[(blockData$blockNum == 3) & (blockData$id %in% expPara$id)],
                  auc = blockData$AUC[(blockData$blockNum == 3) & (blockData$id %in% expPara$id)],
                  cond = expPara$condition)
data = data[expPara$id %in% useID,]
ggplot(data, aes(nQuit)) + geom_dotplot(binwidth = 1) + facet_grid(~cond)
ggplot(data, aes(auc, stdWd)) + geom_point() + facet_grid(~cond)

ggplot(data, aes(tau, cvWd)) + geom_point() + facet_grid(~cond)


# plot correlation
data =  data.frame(tau = expPara$tau,
                   cvWd = blockData$cvWd[(blockData$blockNum == 3) & (blockData$id %in% expPara$id)],
                   cond = expPara$condition)
data = data[expPara$id,]
getCorrelation(data)
# I think I shoud test them, a while range of famaily ?
fit = lm(data = data[data$cond == "LP" & (data$tau) < 14 & (data$cvWd < 1)], cvWd ~ tau + I(1 / auc))
summary(fit)

# is it still true, 
cor.test(log(data$tau[data$cond == "LP"]), log(data$cvWd[data$cond== "LP"] + 2))

# since stdQuitTime is NA is the participant didn't quit
# I think we should look at these point right?
data = data.frame(x = expPara$tau, y = blockData$cvQuitTime[(blockData$blockNum == 3) & (blockData$id %in% expPara$id)],
                  cond = expPara$condition)
data = data[expPara$id %in% useID & !is.na(data$y),]
p = plotCorrelation(data, "grey", T)
p + xlab("Tau rank") + ylab("quitTime std rank") + saveTheme
ggsave("varQuitTime.png", width = 8, height = 4)
plot(blockData$muQuitTime[blockData$condition == "HP"],  blockData$stdQuitTime[blockData$condition == "HP"])

ggplot(data, aes(log(x), log(y))) + geom_point() + facet_grid(~cond)

# 
for(pIdx in 1 : length(paras)){
  para = paras[pIdx]
  # plotParaAUC(expPara, para, blockData, useID)
  input = data.frame(x = expPara[[para]], y = blockData$AUC, cond = blockData$condition)
  input = input[expPara$id %in% useID,]
  p = plotCorrelation(input, paraColors[i], 1)
  xName = sprintf("%sRank", para)
  yName = "AUC-Rank"
  p + xlab(xName) + ylab(yName) + saveTheme
  ggsave(sprintf("%s/AUC_%s.pdf", saveDir, para), width =8 ,height = 4)
}



# 
nPara = length(paras)
for(i in 1 : nPara){
  paraX = paras[i]
  for(j in (i+1) : nPara){
    paraY = paras[j]
    plotParaPara(expPara, paraX, paraY, useID)
    fileName = sprintf("figures/expParaAnalysis/cor_%s_%s.pdf", paraX, paraY)
    ggsave(fileName, width = 6, height = 4)
  }
}



# divide into more


# cvQuitTime
data = data.frame(expPara$tau,
                  y = blockData$cvQuitTime[blockData$blockNum == 1],
                  expPara$condition)
data = data[expPara$id %in% useID & !is.na(data$y),]
p = plotCorrelation(data, "grey", T)
p + xlab("Tau rank") + ylab("cvQuitTime rank") + saveTheme
ggsave("cvQuitTime.png", width = 8, height = 4)

# wtwEarly
data = data.frame(expPara$zeroPoint,
                  y = blockData$wtwEarly[blockData$blockNum == 1],
                  expPara$condition)
data = data[expPara$id %in% useID & !is.na(data$y),]
p = plotCorrelation(data, "grey", T)
p + xlab("zeroPoint rank") + ylab("wtwEarly rank") + saveTheme
ggsave("wtwEarly.png", width = 8, height = 4)


