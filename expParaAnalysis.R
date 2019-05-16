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
parentDir = ifelse(dataType == "block", "genData/expModelFitting", "genData/expModelFittingSub")
dirName = sprintf("%s/%s",parentDir, modelName)
expPara = loadExpPara(paras, dirName)
useID = getUseID(expPara, paras)

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

expPara$condition = summaryData$condition[summaryData$id %in% expPara$id]
for(i in 1 : length(paras)){
  para = paras[i]
  paraColor = paraColors[[i]]
  p = ggplot(expPara[expPara$id %in% useID,], aes_string(para)) + geom_histogram(bins = 6, fill = paraColor) + facet_grid(~condition)+
    saveTheme + ylab("Count") + xlab(capitalize(para))
  if(para == "gamma") p = p + scale_x_continuous(breaks = c(0.7, 0.8, 0.9, 1.0))
  if(para == "QwaitIni") p = p + scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10))
  fileName = sprintf("figures/expParaAnalysis/hist_%s.pdf", para)
  ggsave(fileName, width = 6, height = 3)
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

# trait analysis
personality = read.csv("data/SDGdataset.csv")
plotData = data.frame(expPara, personality)
plotData$stress = hdrData$stress
traits = c("Delay.of.Gratification", "Barratt.Impulsiveness",
           "Intolerance.of.Uncertainty", "Trait.Anxiety..STAIT.")
traitNames = c("DelayGra", "Impulsive",
           "Uncertain", "Anxiety")
plotData = plotData[plotData$id %in% useID,]
rhoTableHP = matrix(NA, nrow = length(paras), ncol = length(traits))
rhoTableLP = matrix(NA, nrow = length(paras), ncol = length(traits))
pTableHP = matrix(NA, nrow = length(paras), ncol = length(traits))
pTableLP = matrix(NA, nrow = length(paras), ncol = length(traits))
for(pIdx in 1 : length(paras)){
  para = paras[pIdx]
  paraColor = paraColors[pIdx]
  for(trIdx in 1 : length(traits)){
    trait = traits[trIdx]
    traitName = traitNames[trIdx]
    input = data.frame(personality[summaryData$id %in% useID,trait], expPara[expPara$id %in% useID,para],
                       summaryData$condition[summaryData$id %in% useID])
    junk = getCorrelation(input)
    # xName = sprintf("%sRank", para)
    # yName = sprintf("%sRank", traitName)
    # p + xlab(xName) + ylab(yName) + saveTheme
    # ggsave(sprintf("figures/expParaAnalysis/%s/%s_%s.png", modelName, traitName, para), width =8 ,height = 4)
    rhoTableHP[pIdx, trIdx] = junk$rhos[1]
    rhoTableLP[pIdx, trIdx] = junk$rhos[2]
    pTableHP[pIdx, trIdx] = junk$p[1]
    pTableLP[pIdx, trIdx] = junk$p[2]   
  }
}
# a lot of ties
# plot correlations 
rownames(pTableHP) = paras
colnames(pTableHP) = traitNames
rownames(pTableLP) = paras
colnames(pTableLP) = traitNames
rownames(rhoTableHP) = paras
colnames(rhoTableHP) = traitNames
rownames(rhoTableLP) = paras
colnames(rhoTableLP) = traitNames
library("corrplot")
png('spHP.png')
corrplot(rhoTableHP, 
         p.mat = pTableHP, 
         is.corr = F, 
         method = "color",
         insig = "label_sig") 
mtext("parameter", side = 2, line = 1.5, cex = 1.5)
mtext("personality", side = 3, line = 2, cex = 1.5, at = 2.5)
dev.off()

# divide into more
data = data.frame(expPara$tau,
                  y = blockData$varQuitTime[blockData$blockNum == 1],
                  expPara$condition)
data = data[expPara$id %in% useID & !is.na(data$y),]
p = plotCorrelation(data, "grey", T)
p + xlab("Tau rank") + ylab("varQuitTime rank") + saveTheme
ggsave("varQuitTime.png", width = 8, height = 4)

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


