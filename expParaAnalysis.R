library("ggplot2")
source("subFxs/plotThemes.R")
source("subFxs/loadFxs.R")
source("subFxs/helpFxs.R")
library("dplyr")
library("Hmisc")
load("wtwSettings.RData")

# load blockdata data
load("genData/expDataAnalysis/blockData.RData")
blockData = blockData[blockData$blockNum == 1,]
idList = unique(blockData$id) 
n = length(idList)


plotParaAUC = function(expPara, paraName, blockData, useID){
  paraColor = paraColors[[paraName]]
  # subset data
  expPara = expPara[(expPara$id %in% useID),]
  blockData = blockData[(blockData$id %in% useID),]
  
  # construct AUCRank
  expPara$AUC = blockData$AUC
  tempt = blockData %>% group_by(condition) %>% mutate(AUCRank = rank(AUC))
  expPara$AUCRank = tempt$AUCRank
  # construct paraRank
  rankName = sprintf("%sRank", paraName)
  expParaHP =expPara[expPara$condition == "HP", ]
  expParaHP[[rankName]] = rank(expParaHP[,paraName])
  expParaLP =expPara[expPara$condition == "LP", ]
  expParaLP[[rankName]] = rank(expParaLP[,paraName])
  tempt = rbind(expParaLP, expParaHP)
  junk = lapply(1 : nrow(expPara), function(i) tempt[tempt$id == expPara$id[i],rankName])
  expPara[[rankName]] = junk
  # calculate correlations
  corHP = cor.test(expPara[expPara$condition == "HP", paraName], expPara[expPara$condition == "HP",]$AUC, method = "spearman")
  corLP = cor.test(expPara[expPara$condition == "LP", paraName], expPara[expPara$condition == "LP",]$AUC, method = "spearman")
  rhoHP = round(corHP$estimate, 3)
  rhoLP= round(corLP$estimate, 3)
  pHP = round(corHP$p.value, 3)
  pLP = round(corLP$p.value, 3)
  textColors = c(ifelse(pHP < 0.05, "red", "blue"), ifelse(pLP < 0.05, "red", "blue"))
  textData = data.frame(label = c(paste(rhoHP, "(p =", pHP, ")"), paste(rhoLP, "(p =", pLP, ")")),
                        condition = c("HP", "LP"), color = textColors)
  # prepare plotData
  plotData = expPara
  plotData[,rankName] = unlist(plotData[,rankName])
  
  p = ggplot(plotData, aes_string(rankName, "AUCRank")) + geom_point(size = 4, color = paraColor, fill = paraColor) +
    facet_grid(~condition)  +
    geom_text(data = textData,aes(x = -Inf,y = -Inf, label = label),
    hjust   = -0.2,vjust = -1,color = "blue",size = 5, fontface = 2, color = textColors) + saveTheme  +
    ylim(c(-8, 68)) + ylab("AUC rank") + xlab(capitalize(rankName))
  print(p)
} 


plotParaPara = function(expPara, paraX, paraY, useID){
  expPara = expPara[expPara$id %in% useID,]
  # calculate correlations
  corHP = cor.test(expPara[expPara$condition == "HP", paraX], expPara[expPara$condition == "HP", paraY], method = "spearman")
  corLP = cor.test(expPara[expPara$condition == "LP", paraX], expPara[expPara$condition == "LP", paraY], method = "spearman")
  rhoHP = round(corHP$estimate, 3)
  rhoLP= round(corLP$estimate, 3)
  pHP = round(corHP$p.value, 3)
  pLP = round(corLP$p.value, 3)
  textColors = c(ifelse(pHP < 0.05, "red", "blue"), ifelse(pLP < 0.05, "red", "blue"))
  textData = data.frame(label = c(paste(rhoHP, "(p =", pHP, ")"), paste(rhoLP, "(p =", pLP, ")")),
                        condition = c("HP", "LP"))
  plotData = expPara
  p = ggplot(plotData, aes_string(paraX, paraY)) + geom_point() +
    facet_grid(~condition)  +
    geom_text(data = textData,aes(x = -Inf,y = -Inf, label = label),
              hjust   = -0.2,vjust = -1,color = "blue",size = 5, fontface = 2) +
    saveTheme + ylab(capitalize(paraY)) + xlab(capitalize(paraX))
  print(p)
}

dir.create("figures/expParaAnalysis")
saveDir = sprintf("figures/expParaAnalysis/%s", modelName)
dir.create(saveDir)
# full_model
modelName = "curiosityTrial"
paras = getParas(modelName)
expPara = loadExpPara(modelName, paras)
useID = getUseID(blockData, expPara, paras)
for(pIdx in 1 : length(paras)){
  para = paras[pIdx]
  plotParaAUC(expPara, para, blockData, useID)
  ggsave(sprintf("%s/AUC_%s.pdf", saveDir, para), width =8 ,height = 4)
}


for(i in 1 : length(paras)){
  para = paras[i]
  paraColor = paraColors[[para]]
  p = ggplot(expPara, aes_string(para)) + geom_histogram(bins = 6, fill = paraColor) + facet_grid(~condition)+
    saveTheme + ylab("Count") + xlab(capitalize(para))
  if(para == "gamma") p = p + scale_x_continuous(breaks = c(0.7, 0.8, 0.9, 1.0))
  if(para == "QwaitIni") p = p + scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10))
  fileName = sprintf("figures/expParaAnalysis/hist_%s.pdf", para)
  ggsave(fileName, width = 6, height = 3)
}

# plot diagonise statistics
meanData = as.double(dplyr::summarise(expPara[expPara$id %in% useID,], phi = mean(phiRhat), tau = mean(tauRhat),
                  gamma = mean(gammaRhat), QwaitIni = mean(QwaitIniRhat)))
stdData = as.double(dplyr::summarise(expPara[expPara$id %in% useID,], phi = sd(phiRhat), tau = sd(tauRhat),
                          gamma = sd(gammaRhat), QwaitIni = sd(QwaitIniRhat)))

plotData = data.frame(para = factor(paras, levels = paras), mean = meanData, std = stdData,
                      ymin = meanData - stdData, ymax = meanData + stdData
                      )
ggplot(plotData, aes(para, mean))+ geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = ymin,ymax = ymax), width = 0.2) + saveTheme +
  xlab("Parameter")+ ylab("Rhat")
ggsave("figures/expParaAnalysis/Rhat.pdf", width = 6, height = 4)



meanData = as.double(dplyr::summarise(expPara[expPara$id %in% useID,], phi = mean(phiEffe) / 5000,
                                      tau = mean(tauEffe) / 5000,
                                      gamma = mean(gammaEffe) / 5000, QwaitIni = mean(QwaitIniEffe)/ 5000))
stdData = as.double(dplyr::summarise(expPara[expPara$id %in% useID,], phi = sd(phiEffe) / (5000^2),
                                     tau = sd(tauEffe) / (5000^2),
                                     gamma = sd(gammaEffe) / (5000^2), QwaitIni = sd(QwaitIniEffe) / (5000^2)))

plotData = data.frame(para = factor(paras, levels = paras), mean = meanData, std = stdData,
                      ymin = meanData - stdData, ymax = meanData + stdData
)
ggplot(plotData, aes(para, mean))+ geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = ymin,ymax = ymax), width = 0.2) + saveTheme +
  xlab("Parameter")+ ylab("Effective Sample Proportion / %")
ggsave("figures/expParaAnalysis/Effe.pdf", width = 6, height = 4)

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
plotData = plotData[plotData$id %in% useID,]
para = "tau"
for(trIdx in 1 : length(traits)){
  trait = traits[trIdx]
  p = ggplot(plotData,
             aes_string(para, trait)) +
    geom_point() + facet_grid(~condition)
  print(p)
  print(cor.test(plotData[plotData$condition == "HP",trait],
                 plotData[plotData$condition == "HP",para], method = "spearman"))
  print(cor.test(plotData[plotData$condition == "LP",trait],
                 plotData[plotData$condition == "LP",para], method = "spearman"))
  readline("continue")
}


data = data.frame(expPara$tau,
                  y = blockData$varQuitTime[blockData$blockNum == 3],
                  expPara$condition)
data = data[expPara$id %in% useID & !is.na(data$y),]
p = plotCorrelation(data, "green", "spearman", T)
p + xlab("Tau rank") + ylab("varQuitTime rank") + saveTheme
ggsave("varQuitTime.png", width = 8, height = 4)

data = data.frame(expPara$tau,
                  y = personality$Delay.of.Gratification,
                  expPara$condition)
data = data[expPara$id %in% useID,]
p = plotCorrelation(data, "green", "spearman", T)
p + xlab("Tau rank") + ylab("Impulsiveness") + saveTheme
ggsave("delay.png", width = 8, height = 4)