library("ggplot2")
library("dplyr")
library("tidyr")
library("Hmisc")
library("coin")
source("subFxs/plotThemes.R")
source("subFxs/loadFxs.R") # load blockData and expPara
source("subFxs/helpFxs.R") # getParas
source("subFxs/analysisFxs.R") # plotCorrelation and getCorrelation
load("wtwSettings.RData")

# input 
dataType = "sess"
modelName = "PR"

# create output directories
if(dataType == "sess"){
  dir.create("figures/expParaAnalysisSub")
  saveDir = sprintf("figures/expParaAnalysisSub/%s", modelName)
  dir.create(saveDir)
}else{
  dir.create("figures/expParaAnalysis")
  saveDir = sprintf("figures/expParaAnalysis/%s", modelName)
  dir.create(saveDir)
}

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

# load expPara
paras = getParas(modelName)
nPara = length(paras)
parentDir = ifelse(dataType == "block", "genData/expModelFitting", "genData/expModelFittingSub")
dirName = sprintf("%s/%s",parentDir, modelName)
tempt = loadExpPara(paras, dirName)
useID = getUseID(tempt, paras)
expPara = merge(x=tempt[,c(paras, "id")],y=summaryData, by="id",all.x=TRUE)
expPara$AUC1 = blockData$AUC[blockData$id %in% expPara$id & blockData$blockNum == 1]
expPara$AUC2 = blockData$AUC[blockData$id %in% expPara$id & blockData$blockNum == 2]
expPara$AUC1_2  = (expPara$AUC1 + expPara$AUC2) / 2
# plot them separately is ugly
traitParaCorr = vector(mode = "list", length = nPara)
for(pIdx in 1 : length(paras)){
  para = paras[pIdx]
  # plotParaAUC(expPara, para, blockData, useID)
  input = data.frame(x = expPara[[para]], y = expPara$AUC,
                     cond = expPara$condition)
  input = input[expPara$id %in% useID,]
  traitParaCorr[[pIdx]] = getCorrelation(input)
}
rhoTable = sapply(1:2, function(j) sapply(1: nPara, function(i) traitParaCorr[[i]]$rhos[j]))
rownames(rhoTable) = c("LR", "LP", "Tau", "Gamma", "P")
colnames(rhoTable) = c("AUC_HP", "AUC_LP")
pTable = sapply(1:2, function(j) sapply(1: nPara, function(i) traitParaCorr[[i]]$ps[j]))

library("corrplot")
col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061")))

parentDir = ifelse(dataType == "sess", "figures/expParaAnalysisSub", "figures/expParaAnalysis")

fileName = sprintf("%s/%s/AUCPara_%s.png", parentDir, modelName, dataType)
png(fileName)
corrplot(t(rhoTable), 
         p.mat = t(pTable), 
         is.corr = T, 
         method = "color",
         insig = "p-value",
         tl.col = "black", tl.srt = 0, tl.cex = 1.5,
         col = col2(50)) 
dev.off()


# plotParaAUC(expPara, para, blockData, useID)
input = data.frame(expPara[paras], AUC = expPara$AUC,
                   cond = expPara$condition, id = expPara$id) %>% filter(id %in% useID)
junk = input %>% group_by(cond) %>% dplyr::summarise(muAUC = mean(AUC), stdAUC = sd(AUC))
input = input %>% mutate(AUC = AUC - junk$muAUC[1] * (cond == "HP") - junk$muAUC[2] * (cond == "LP"))
inputHP = input[input$cond == "HP",] %>% mutate(phi = deMean(phi), phiP = deMean(phiP), tau = deMean(tau), gamma = deMean(gamma), zeroPoint = deMean(zeroPoint)) 
inputLP = input[input$cond == "LP",] %>% mutate(phi = deMean(phi), phiP = deMean(phiP), tau = deMean(tau), gamma = deMean(gamma), zeroPoint = deMean(zeroPoint))
plotData = rbind(inputHP, inputLP)

names = c(paras, "AUC")
rankData =  lapply(1 : length(names), function(i) rank(plotData[,names[i]]))
plotData[,names] = rankData

ps = vector(length = nPara)
rhos = vector(length = nPara)
textColors = vector(length = nPara)
for(i in 1 : nPara){
  tempt = cor.test(plotData[,paras[i]], plotData$AUC, method = "kendall")
  textColors[i] = ifelse(tempt$p.value < 0.05, "black", "grey")
  rhos[i] = round(tempt$estimate, 2)
  if(tempt$p.value < 0.001){
    ps[i] = "p < 0.001"
  }else{
    ps[i] = sprintf("p = %.3f", round(tempt$p.value, 3)) 
  }
}
textData = data.frame(label = paste0(rhos, "(", ps, ")"), textColors, para=factor(paras, levels = paras, labels = c("LR", "LP", "tau", "gamma", "P")))
plotData %>% gather(-c("cond", "id", "AUC"), key = "para", value = value) %>%
  mutate(para = factor(para, levels = paras, labels = c("LR", "LP", "tau", "gamma", "P"))) %>%
  ggplot(aes(value, AUC, color = para)) + geom_point(size = 3) +
    facet_wrap(~para, nrow = 1, labeller = label_parsed) +
  scale_color_manual(values = unlist(paraColors[1 : length(paras)])) + myTheme + ylab("Para rank") +
  ylab("WTW average") + xlab("Parameter") + scale_x_continuous(breaks = c(0, 60), limits = c(0, 60)) + 
  scale_y_continuous(breaks = c(0, 60), limits = c(0, 60)) + theme(legend.position = "None") + geom_text(data = textData,aes(x = -Inf,y = -Inf, label = label),
            hjust = -0.2 ,vjust = -12,color = "#252525", size = 5, fontface = 2) + geom_smooth(method = lm, se = F)
ggsave(sprintf("%s/%s/single_AUC.png", "figures/expParaAnalysisSub", modelName), width = 10 ,height = 3) 


# plot for HP
for(c in 1:2){
  thisCond = conditions[c]
  input = data.frame(expPara[paras], AUC = expPara$AUC,
                     cond = expPara$condition, id = expPara$id) %>% filter(id %in% useID & cond == thisCond)
  names = c(paras, "AUC")
  rankData =  lapply(1 : length(names), function(i) rank(input[,names[i]]))
  plotData = input
  plotData[,names] = rankData
  
  ps = vector(length = nPara)
  rhos = vector(length = nPara)
  textColors = vector(length = nPara)
  for(i in 1 : nPara){
    tempt = cor.test(input[,paras[i]], input$AUC, method = "kendall")
    textColors[i] = ifelse(tempt$p.value < 0.05, "black", "grey")
    rhos[i] = round(tempt$estimate, 2)
    if(tempt$p.value < 0.001){
      ps[i] = "p < 0.001"
    }else{
      ps[i] = sprintf("p = %.3f", round(tempt$p.value, 3)) 
    }
  }
  textData = data.frame(label = paste0(rhos, "(", ps, ")"), textColors, para=factor(paras, levels = paras, labels = c("LR", "LP", "tau", "gamma", "P")))
  plotData %>% gather(-c("cond", "id", "AUC"), key = "para", value = value) %>%
    mutate(para = factor(para, levels = paras, labels = c("LR", "LP", "tau", "gamma", "P"))) %>%
    ggplot(aes(value, AUC, color = para)) + geom_point(size = 3) +
    facet_wrap(~para, nrow = 1, labeller = label_parsed) +
    scale_color_manual(values = unlist(paraColors[1 : length(paras)])) + myTheme + ylab("Para rank") +
    ylab("AUC rank") + xlab("Parameter rank") + scale_x_continuous(breaks = c(0, 30), limits = c(0, 30)) + 
    scale_y_continuous(breaks = c(0, 30), limits = c(-5, 35)) + theme(legend.position = "None") + geom_text(data = textData,aes(x = -Inf,y = -Inf, label = label),
                                                                                                            hjust = -0.2 ,vjust = -12,color = "#252525", size = 5, fontface = 2) 
  ggsave(sprintf("%s/%s/single_AUC_%s.png", "figures/expParaAnalysisSub", modelName, thisCond), width = 10 ,height = 3) 
}


# add blockData
load("genData/expDataAnalysis/blockData.RData")
expPara$stdQuitTime3 = blockData$stdQuitTime[blockData$blockNum == 3 &
                                               blockData$stress == "no stress"]
expPara$stdWd3 = blockData$stdWd[blockData$blockNum == 3 &
                                   blockData$stress == "no stress"]
expPara$AUCEarly = sessionData$AUCEarly[sessionData$id %in% expPara$id]
# wtwEarly to prior
expPara %>% filter(id %in% useID) %>%ggplot(aes(zeroPoint, AUCEarly)) +
  geom_point() + facet_wrap(~condition, scales = "free")  + myTheme
  

input = data.frame(expPara$zeroPoint, expPara$AUCEarly, expPara$condition)
input = input[expPara$id %in% useID, ]
plotCorrelation(input, dotColor = "grey", F)

para = "tau"
ggplot(expPara, aes(tau, stdQuitTime)) + geom_point() +
  facet_grid(~ condition, scales = "free") 

para = "tau"
ggplot(expPara, aes(tau, stdQuitTime3)) + geom_point() +
  facet_grid(~ condition, scales = "free") 


library("ppcor")
tempt = expPara[expPara$condition == "HP"&
                  expPara$id %in% useID,] 
pcor.test(tempt$tau, tempt$stdWd, tempt$AUC,
          method = "kendall")

tempt = expPara[expPara$condition == "LP"&
                  expPara$id %in% useID,] 
pcor.test(tempt$tau, tempt$stdWd, tempt$AUC,
          method = "kendall")

data = expPara[,c("tau", "stdQuitTime", "condition")]
plotCorrelation(data, "grey", T)
data = expPara[,c("tau", "cvQuitTime", "condition")]
plotCorrelation(data, "grey", T)
data = expPara[,c("tau", "stdWd3", "condition")]
plotCorrelation(data, "grey", T)
data = expPara[,c("tau", "cvWd3", "condition")]
plotCorrelation(data, "grey", T)
data = expPara[,c("phi", "stdWd", "condition")]
plotCorrelation(data, "grey", T)
data = expPara[,c("tau", "cvWd", "condition")]
plotCorrelation(data, "grey", T)

# look linear
# maybe not about you know, the final block, the overall.
fit = lm(data = expPara[expPara$condition == "LP",], cvWd ~ 
           tau + I(1 / AUC))
summary(fit)

# parse behavioral data
load("genData/expDataAnalysis/blockData.RData")
# I will construct a dataset of interest here
summaryData$adaptation = summaryData$wtwEarly -  sapply(1 : nrow(summaryData), function(i) winAUC_[[i]][3])
summaryData$id[summaryData$adaptation > 20 & summaryData$id %in% useID]

data = merge(x = expPara[,1:nPara],
             y = summaryData, by="id",all.x=TRUE) 
data = data[expPara$id %in% useID,]

ggplot(data, aes(phi, adaptation)) + geom_point() + facet_grid(~condition)

# time to stabalize 
sIdx =  15
trialPlots(block2session(trialData[[sIdx]]))
plot(timeStdWd_[[sIdx]])
ggplot(data, aes(auc, stdWd)) + geom_point() + facet_grid(~cond)

ggplot(data, aes(tau, cvWd)) + geom_point() + facet_grid(~cond)


# plot correlation
data =  data.frame(tau = expPara$tau,
                   cvWd = blockData$cvWd[(blockData$blockNum == 3) & (blockData$id %in% expPara$id)],
                   cond = expPara$condition,
                   auc = blockData$AUC[(blockData$blockNum == 3) & (blockData$id %in% expPara$id)])
data = data[expPara$id,]
# I think I shoud test them, a while range of famaily ?
fit = lm(data = data[data$cond == "LP",], cvWd ~ tau + I(1 / auc))
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


# correlation between AUC



# correlation between para
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


