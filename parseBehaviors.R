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
  load("genData/expDataAnalysis/blockData.RData")
  summaryData = data.frame(cbind(sessionData, blockData[blockData$blockNum == 3 & blockData$id %in% expPara$id,
                                                 c("cvWd", "stdWd", "AUC", "cvQuitTime")]))
  colnames(summaryData) = c(names(sessionData),c("cvWd3", "stdWd3", "AUC3", "cvQuitTime3"))
}

# load expPara
paras = getParas(modelName)
nPara = length(paras)
parentDir = ifelse(dataType == "block", "genData/expModelFitting", "genData/expModelFittingSub")
dirName = sprintf("%s/%s",parentDir, modelName)
tempt = loadExpPara(paras, dirName)
useID = getUseID(tempt, paras)
# load summaryData
expPara = merge(x=tempt[,c(paras, "id")],y=summaryData, by="id",all.x=TRUE)


# link noise to tau 
para = "tau"
ggplot(expPara, aes(tau, cvQuitTime)) + geom_point() +
  facet_grid(~ condition, scales = "free") 

para = "tau"
ggplot(expPara, aes(tau, stdQuitTime)) + geom_point() +
  facet_grid(~ condition, scales = "free") 

para = "zeroPoint"
ggplot(expPara, aes(zeroPoint, stdWd)) + geom_point() +
  facet_grid(~ condition, scales = "free") 

para = "tau"
ggplot(expPara, aes(tau, cvWd3)) + geom_point() +
  facet_grid(~ condition, scales = "free") 

para = "tau"
ggplot(expPara, aes(gamma, stdWd)) + geom_point() +
  facet_grid(~ condition, scales = "free") 

para = "tau"
ggplot(expPara, aes(tau, cvWd)) + geom_point() +
  facet_grid(~ condition, scales = "free") 

library("ppcor")
tempt = expPara[expPara$condition == "HP"&
                expPara$id %in% useID,]
pcor.test(tempt$zeroPoint, tempt$stdWd, tempt$AUC,
          method = "kendall")

tempt = expPara[expPara$condition == "LP"&
                  expPara$id %in% useID,]
pcor.test(tempt$gamma, tempt$stdWd, tempt$AUC,
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


