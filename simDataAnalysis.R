# 
library("ggplot2")
library("Hmisc")
source("subFxs/analysisFxs.R")
source("subFxs/helpFxs.R")
dir.create("figures/simDataAnalysis")
dirName = "genData/simulation/full_model"
load(sprintf("%s/trialHPData.RData", dirName))
load(sprintf("%s/trialLPData.RData", dirName))
load(sprintf("%s/simParas.RData", dirName))

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

paras = getParas("full_model")
nPara = length(paras)
for(c in 1 : 2)

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
    fileName = sprintf("figures/simDataAnalysis/AUC_%s_%s.pdf", cond, para)
    ggsave(fileName, width = 3, height = 4)
  }
}


