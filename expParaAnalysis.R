library("ggplot2")
source("subFxs/plotThemes.R")
source("subFxs/loadFxs.R")
source("subFxs/helpFxs.R")
library("dplyr")
library("Hmisc")
library("coin")
source("subFxs/analysisFxs.R")
load("wtwSettings.RData")

# load blockdata data
load("genData/expDataAnalysis/blockData.RData")
blockData = blockData[blockData$blockNum == 1,]
idList = unique(blockData$id) 
n = length(idList)


dir.create("figures/expParaAnalysis")
saveDir = sprintf("figures/expParaAnalysis/%s", modelName)
dir.create(saveDir)
# full_model
modelName = "curiosityTrialRSp"
paras = getParas(modelName)
expPara = loadExpPara(modelName, paras)
useID = getUseID(blockData, expPara, paras)
blockData = blockData[blockData$blockNum == 1, ]
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

expPara$condition = blockData$condition
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

# # plot diagonise statistics
# meanData = as.double(dplyr::summarise(expPara[expPara$id %in% useID,], phi = mean(phiRhat), tau = mean(tauRhat),
#                   gamma = mean(gammaRhat), QwaitIni = mean(QwaitIniRhat)))
# stdData = as.double(dplyr::summarise(expPara[expPara$id %in% useID,], phi = sd(phiRhat), tau = sd(tauRhat),
#                           gamma = sd(gammaRhat), QwaitIni = sd(QwaitIniRhat)))
# 
# plotData = data.frame(para = factor(paras, levels = paras), mean = meanData, std = stdData,
#                       ymin = meanData - stdData, ymax = meanData + stdData
#                       )
# ggplot(plotData, aes(para, mean))+ geom_bar(stat = "identity") +
#   geom_errorbar(aes(ymin = ymin,ymax = ymax), width = 0.2) + saveTheme +
#   xlab("Parameter")+ ylab("Rhat")
# ggsave("figures/expParaAnalysis/Rhat.pdf", width = 6, height = 4)



# meanData = as.double(dplyr::summarise(expPara[expPara$id %in% useID,], phi = mean(phiEffe) / 5000,
#                                       tau = mean(tauEffe) / 5000,
#                                       gamma = mean(gammaEffe) / 5000, QwaitIni = mean(QwaitIniEffe)/ 5000))
# stdData = as.double(dplyr::summarise(expPara[expPara$id %in% useID,], phi = sd(phiEffe) / (5000^2),
#                                      tau = sd(tauEffe) / (5000^2),
#                                      gamma = sd(gammaEffe) / (5000^2), QwaitIni = sd(QwaitIniEffe) / (5000^2)))
# 
# plotData = data.frame(para = factor(paras, levels = paras), mean = meanData, std = stdData,
#                       ymin = meanData - stdData, ymax = meanData + stdData
# )
# ggplot(plotData, aes(para, mean))+ geom_bar(stat = "identity") +
#   geom_errorbar(aes(ymin = ymin,ymax = ymax), width = 0.2) + saveTheme +
#   xlab("Parameter")+ ylab("Effective Sample Proportion / %")
# ggsave("figures/expParaAnalysis/Effe.pdf", width = 6, height = 4)

## parameter correlation
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

    input = data.frame(personality[,trait], expPara[[para]], blockData$condition)
    junk = getCorrelation(input, paraColors[i], T)
    # xName = sprintf("%sRank", para)
    # yName = sprintf("%sRank", traitName)
    # p + xlab(xName) + ylab(yName) + saveTheme
    # ggsave(sprintf("figures/expParaAnalysis/%s/%s_%s.png", modelName, traitName, para), width =8 ,height = 4)
    junk = getCorrelation(input, paraColors[i], T)
    rhoTableHP[pIdx, trIdx] = junk$rhos[1]
    rhoTableLP[pIdx, trIdx] = junk$rhos[2]
    pTableHP[pIdx, trIdx] = junk$p[1]
    pTableLP[pIdx, trIdx] = junk$p[2]   
  }
}

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
                  y = blockData$varQuitTime[blockData$blockNum == 3],
                  expPara$condition)
data = data[expPara$id %in% useID & !is.na(data$y),]
p = plotCorrelation(data, "green", T)
p + xlab("Tau rank") + ylab("varQuitTime rank") + saveTheme
ggsave("varQuitTime.png", width = 8, height = 4)

data = data.frame(expPara$tau,
                  y = personality$Delay.of.Gratification,
                  expPara$condition)
data = data[expPara$id %in% useID,]
p = plotCorrelation(data, "green", "spearman", T)
p + xlab("Tau rank") + ylab("Impulsiveness") + saveTheme
ggsave("delay.png", width = 8, height = 4)