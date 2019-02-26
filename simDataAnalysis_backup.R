# this script analysized the simulation data on the group level
modelName = "monteSteep"
load(sprintf("genData/simulation/%s/realParas.RData", modelName)) # for nComb, paraNames, nPara, nBlock
nBlock = 3 
############ load data and functions #########
# generally loading 
library("ggplot2")
library("dplyr")
library("tidyr")
source(file = './subFxs/plotThemes.R')
load("wtwSettings.RData")

# specific loading 
library('scales')
source('subFxs/analysisFxs.R')
load(sprintf("genData/simulation/%s/trialHPData.RData", modelName)) # trialHPData
load(sprintf("genData/simulation/%s/trialLPData.RData", modelName)) # trialLPData
nRep = 5

# trial-level analysis 
plotTrialwiseData = F
plotKMSC = F
plotWTW = F

# loop over conditions 
for(cIdx in 1 : 2){
  cond = conditions[cIdx]
  if(cIdx == 1){
    trialData = trialHPData
  }else{
    trialData = trialLPData
  }
  # generate arguments for later analysis 
  tMax = tMaxs[cIdx]
  kmGrid = seq(0, tMax, by=0.1) 
  
  # initialize outputs
  totalEarnings_ = matrix(0, c(n,  nRep))
  AUC_ = array(0, c(n, nRep))
  
  # loop over combIdx
  for(combIdx in 1 : nComb){
    realPara = realParas[combIdx,]
    # loop over blocks
    for(bkIdx in 1 : nBlock){
      # loop over repitations
      for(rIdx in 1 : nRep){
        # select data
        thisTrialData = trialData[[simIdx[combIdx, rIdx]]]
        select = which(thisTrialData$sellTime <= bkIdx * blockSecs &
                         thisTrialData$sellTime > (bkIdx - 1) * blockSecs)
        thisTrialData$trialNum = thisTrialData$trialNum[select]
        thisTrialData$trialEarnings = thisTrialData$trialEarnings[select]
        thisTrialData$timeWaited = thisTrialData$timeWaited[select]
        thisTrialData$sellTime = thisTrialData$sellTime[select]
        thisTrialData$scheduledWait = thisTrialData$scheduledWait[select] 
        
        # generate arguments for later analysis 
        label = sprintf('%s %s', conditions[cIdx], paste0(paraNames, realPara, collapse = " "))
        
        # summarise totalEarnings 
        totalEarnings_[combIdx, bkIdx, rIdx] =  sum(thisTrialData$trialEarnings)
        
        # plot trial-by-trial data
        if (plotTrialwiseData) {
          trialPlots(thisTrialData,label)
        }
        # survival analysis
        kmscResults = kmsc(thisTrialData,tMax,label,plotKMSC,kmGrid)
        AUC_[combIdx, bkIdx, rIdx] = kmscResults[['auc']]
        
        # WTW time series
        wtwCeiling = tMax
        wtwtsResults = wtwTS(thisTrialData, tGrid, wtwCeiling, label, plotWTW)
      }# end of repetitions
    }# end of blocks
  }# end of para combinations 
  # organize and save blockData
  dim(AUC_) = c(nComb * nBlock, nRep) 
  dim(totalEarnings_) = c(nComb * nBlock, nRep) 
  AUC = rowSums(AUC_) / nRep
  totalEarnings = rowSums(totalEarnings_) / nRep
  blockData = data.frame(combIdx = rep(1 : combIdx, nBlock), blockNum = rep(1 : nBlock, each = nComb),
                         condition = factor(rep(cond, each = nComb * nBlock), levels = c("HP", "LP")),
                         AUC = AUC, 
                         totalEarnings = totalEarnings)
  if(cIdx == 1){
    blockHPData = blockData
  }else{
    blockLPData = blockData
  }
}

####### plot distribution of totalEarnings
realParaRanks = apply(realParas, MARGIN = 2, function(x) match(x, sort(unique(x))))
blockData = data.frame(rbind(blockHPData, blockLPData), rbind(realParaRanks, realParaRanks))
colnames(blockData)[(ncol(blockData) - nPara + 1) : ncol(blockData) ] = paraNames


############ summarise para effects on total earnings ###########
paraValues = 1:nValue
paraData = data.frame(paraName = rep(paraNames,each = nValue * 2 * nBlock),
                      paraRank = rep(rep(1 : nValue, each =  2 * nBlock), nPara),
                      condition = rep(rep(c("HP", "LP"), each = nBlock), nValue),
                      blockNum = rep(1 : nBlock, nValue * 2 * nPara)
                      ) # 2 is nCondition
paraData$paraName = factor(paraData$paraName, levels = paraNames)

# mu
tempt = vector("list", length = nPara)
for(i in 1 : nPara){
  tempt[[i]] = blockData %>%
    group_by_at(vars(paraNames[i], condition, blockNum)) %>%
    summarize_at(vars(AUC:totalEarnings), mean)
}
junk = ldply (tempt, data.frame)
paraData = data.frame(paraData, muAUC = junk$AUC, muEarn = junk$totalEarnings)
# std
tempt = vector("list", length = nPara)
for(i in 1 : nPara){
  tempt[[i]] = blockData %>%
    group_by_at(vars(condition, blockNum, paraNames[i],)) %>%
    summarize_at(vars(AUC:totalEarnings), sd)
}
junk = ldply (tempt, data.frame)
paraData = data.frame(paraData, stdAUC = junk$AUC, stdEarn = junk$totalEarnings)
# 
paraData$maxAUC = paraData$muAUC + paraData$stdAUC
paraData$minAUC = paraData$muAUC - paraData$stdAUC
paraData$maxEarn = paraData$muEarn + paraData$stdEarn
paraData$minEarn= paraData$muEarn - paraData$stdEarn
paraData$paraRank = as.factor(paraData$paraRank)
# plot earn
dir.create("figures")
dir.create("figures/simDataAnalysis")
dir.create(sprintf("figures/simDataAnalysis/%s", modelName))
for(c in 1:2){
  cond = conditions[c]
  ggplot(paraData[paraData$condition == cond,], aes(paraRank, muAUC)) +
    geom_bar(stat = "identity", width=0.5, fill = conditionColors[c]) + geom_errorbar(aes(ymin = minAUC, ymax = maxAUC), width=.2)+
    facet_wrap(paraName ~ blockNum, nrow = 2)+ saveTheme +
    xlab("Parameter value") + ylab("AUC / s") + ggtitle(cond) 
  fileName = sprintf("figures/simDataAnalysis/%s/paraAUCEffect%s.pdf", modelName, cond)
  ggsave(fileName, width = 6, height = 6) 
}

# how to define different 
realParaRanks = apply(realParas, MARGIN = 2, function(x) match(x, sort(unique(x))))
blockData = data.frame(rbind(blockHPData, blockLPData), rbind(realParas, realParas))
colnames(blockData)[(ncol(blockData) - nPara + 1) : ncol(blockData) ] = paraNames
blockData$steep = as.factor(blockData$steep)
blockData$phi = as.factor(blockData$phi)

ggplot(blockData, aes(steep, AUC)) + geom_bar(stat = "identity", width=0.5, fill = conditionColors[c])+
  facet_wrap(phi ~ blockNum, nrow = 3)

