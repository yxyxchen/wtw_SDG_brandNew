# 
library("ggplot2")
library("Hmisc")
source("subFxs/analysisFxs.R")
source("subFxs/helpFxs.R")
source("subFxs/loadFxs.R")
load("wtwSettings.RData")

modelName = "R_learning"
dir.create("figures/simDataAnalysis")
dirName = sprintf("genData/simulation/%s", modelName)
load(sprintf("%s/trialHPData.RData", dirName))
load(sprintf("%s/trialLPData.RData", dirName))
load(sprintf("%s/simParas.RData", dirName))
source("subFxs/plotThemes.R")
dir.create("figures/simDataAnalysis")
dirName = sprintf("figures/simDataAnalysis/%s", modelName)
dir.create(dirName)
for(c in 1 : 2){
  cond = conditions[c]
  if(cond == "HP") trialData = trialHPData else trialData = trialLPData
  tMax = ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  # calculate AUC and timeWaited
  plotKMSC = F
  label = ""
  kmGrid = seq(0, tMax, by=0.1) 
  # initialize 
  totalEarnings_ = matrix(0, nComb, nRep)
  AUC_ = matrix(0, nComb, nRep)
  for(sIdx in 1 : nComb){
    for(rIdx in 1 : nRep){
      thisTrialData = trialData[[simNo[sIdx, rIdx]]]
      kmscResults = kmsc(thisTrialData,tMax,label,plotKMSC,kmGrid)
      AUC_[sIdx, rIdx] = kmscResults[['auc']]
    }
  }
  AUC = apply(AUC_, MARGIN = 1, FUN = mean)
  if(cond == "HP") AUCHP = AUC else AUCLP = AUC
}

paras = getParas("R_learning")
nPara = length(paras)
for(c in 1 : 2){
  cond = conditions[c]
  condColor = conditionColors[c]
  ylimit = ifelse(cond == "HP", 25, 45)
  if(cond == "HP") AUC = AUCHP else AUC = AUCLP
  tempt = data.frame(paraComb, AUC = AUC)
  for(i in 1 : nPara){
    para = paras[i]
    tempt1 = tempt %>% group_by_at(vars(para)) %>% summarise(mu = mean(AUC))
    tempt2 = tempt %>% group_by_at(vars(para)) %>% summarise(std = sd(AUC))
    plotData = data.frame(tempt1, std = tempt2$std) 
    plotData[[para]]= as.factor(plotData[[para]])
    plotData$ymin = plotData$mu - plotData$std
    plotData$ymax = plotData$mu + plotData$std
    ggplot(plotData, aes_string(para, "mu")) +
      geom_bar(stat = "identity", color = condColor, fill = condColor) +
      saveTheme + xlab(capitalize(para)) + ylab("AUC / min") + ylim(c(-3, ylimit)) +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2)
    fileName = sprintf("figures/simDataAnalysis/%s/AUC_%s_%s.pdf", modelName, cond, para)
    ggsave(fileName, width = 3, height = 4)
  }
}



## look at corrletion
modelName = "full_model"
paras = getParas(modelName )
nPara = length(paras)
load(sprintf("genData/simulation/%s/simParas.RData", modelName))
for(c in 1:2){
  cond = conditions[c]
  condColor = conditionColors[c]
  simPara_ = loadSimPara_(modelName, paras, cond)
  simPara = matrix(simPara_, 810, 4*5)
  if(c == 1){
    simParaHP = simPara
  }else{
    simParaLP = simPara
  }
}
simPara = rbind(simParaHP, simParaLP)
simPara = as.data.frame(simPara)
colnames(simPara) = c(paras,  paste0(paras, "SD"), paste0(paras, "Effe"), paste0(paras, "Rhat"),
                      paste0(paras, "real"))
simPara$condition = rep(c("HP", "LP"), each = nrow(simParaHP))

# choose only useful data
RhatCols = which(str_detect(colnames(simPara), "hat"))[1 : length(paras)]
EffeCols = which(str_detect(colnames(simPara), "Effe"))[1 : length(paras)]
simParaUse = simPara[apply(simPara[,RhatCols] < 1.1, MARGIN = 1, sum) == length(paras) & 
                 apply(simPara[,EffeCols] >100, MARGIN = 1, sum) == length(paras), ]

plotParaPara = function(expPara, paraX, paraY){
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

for(i in 1 : (nPara-1)){
  paraX = paras[i]
  for(j in (i+1): nPara){
    paraY = paras[j]
    plotParaPara(simParaUse, paraX, paraY)
    fileName = sprintf("figures/simDataAnalysis/cor_%s_%s.pdf", paraX, paraY)
    ggsave(fileName, width = 6, height = 4)
  }
}

library("corrplot")
tempt = simParaHP[,1:4]
colnames(tempt) = paras
M = cor(tempt, method = "spearman")
pdf('figures/simDataAnalysis/corrplotHP.pdf')
corrplot(M, method = "color", type = "upper", outline = T, tl.col = "black",
          tl.cex = 2, cl.pos = "b", cl.cex = 1.2) 
dev.off()

tempt = simParaLP[,1:4]
colnames(tempt) = paras
M = cor(tempt, method = "spearman")
pdf('figures/simDataAnalysis/corrplotLP.pdf')
corrplot(M, method = "color", type = "upper", outline = T, tl.col = "black",
         tl.cex = 2, cl.pos = "b", cl.cex = 1.2) 
dev.off()

c= 1
cond = conditions[c]
condColor = conditionColors[c]
tempt = simParaUse[simParaUse$condition == "LP", ]
##### recaluclate the recover accuracy
for(i in 1 : nPara){
  para = paras[i]
  realName = sprintf("%sreal", para)
  tempt1 = tempt %>% group_by_at(vars(realName)) %>% summarise(mu = mean(get(para)))
  tempt2 = tempt %>% group_by_at(vars(realName)) %>% summarise(std = sd(get(para)))
  plotData = data.frame(tempt1, std = tempt2$std) 
  plotData$ymin = plotData$mu - plotData$std
  plotData$ymax = plotData$mu + plotData$std
  if(para == "phi" || para == "gamma"){
    wid = 0.01
  }else if(para == "tau" ){
    wid = 2
  }else wid = 0.5
  
  ggplot(plotData, aes_string(realName, "mu")) +
    geom_bar(stat = "identity", fill = condColor,color = condColor) +
    saveTheme + xlab(capitalize(para)) + ylab("Phi estimation")  +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = wid) +
    geom_point(aes_string(realName, realName), shape = 17, size = 5, fill = "black")
  fileName = sprintf("figures/simDataAnalysis/recovery_%s_%s.pdf", cond, para)
  ggsave(fileName, width = 4, height = 4)
}
