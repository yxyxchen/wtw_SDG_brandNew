load("genData/expDataAnalysis/blockData.RData")

library("dplyr")
library("tidyr")
source("subFxs/loadFxs.R")
source("subFxs/plotThemes.R")

data = loadAllData()
hdrData = data$hdrData
trialData = data$trialData
idList = hdrData$ID
n = length(hdrData$ID)

# para follows location(mu), scale(sigma) and shape (-k)
library("lmomco")
gpaParas = vec2par(c(0, 2, -4), type = "gpa")
# when t = 40, CDF = 2/3
cdfgpa(40, gpaParas)
# calculate theoretic quatitls from 1/9 to 1
theoreticQs = quagpa(seq(1/9, 6/9, by = 1/9),
                     gpaParas, paracheck=TRUE)
for(sIdx in 1:n){
  thisID = idList[sIdx]
  thisTrialData = trialData[[thisID]]
  thisTrialData = thisTrialData[thisTrialData$blockNum == 1,]
  cond = unique(thisTrialData$condition)
  if(cond == "LP"){
    thisSeq = thisTrialData$scheduledWait
    empiricalQs = as.double(quantile(thisSeq, probs = seq(1/9, 6/9, by = 1/9)))
    plotData = data.frame(theoreticQs, empiricalQs)
    p = ggplot(plotData, aes(theoreticQs, empiricalQs)) + geom_point() +
      geom_abline(intercept = 0, slope = 1) + ylim(c(-2,42)) + xlim(c(-2,42)) +
      ggtitle(thisID)
    print(p)
    print(thisSeq)
    readline("continue")
  }
}


######## calculate waiting threshold using an ideal observer
######## estimated from the empirical sequences
#
isPlot = T
# inputs
tMaxs = c(20, 40) # trial durations
iti = 2 # iti duration in secs
tGrid = seq(0, blockSecs, 0.1)
stepDuration = 0.1
# time ticks within a trial for timeEarnings or wtw analysis
trialTicks = list(
  'HP' = round(seq(0, tMaxs[1], by = stepDuration), 1),
  'LP' = round(seq(0, tMaxs[2], by = stepDuration), 1)
)
trialGapValues = list(
  "HP" = round(seq(stepDuration, tMaxs[1], by = stepDuration), 2),
  "LP" = round(seq(stepDuration, tMaxs[2], by = stepDuration), 2)
)
trialGapIdxs = list(
  "HP" = 1 : length(trialGapValues$HP),
  "LP" = 1 : length(trialGapValues$LP)
)
# PDF and CDF
empOptimRewardRates = vector(length = n)
empOptimWaitTimes = vector(length = n)
for(sIdx in 1:n){
  thisID = idList[sIdx]
  thisTrialData = trialData[[thisID]]
  thisTrialData = thisTrialData[thisTrialData$blockNum == 1,]
  thisSeq = thisTrialData$scheduledWait
  cond = unique(thisTrialData$condition)
  
  thisTrialGapValues = trialGapValues[[cond]]
  tMax = tMaxs[conditions == cond]
  nTimeStep = tMax / stepDuration
  delayCDF = sapply(1 : nTimeStep, function(i) sum(thisSeq <= thisTrialGapValues[i])/
                      length(thisSeq)) 
  delayPDF =  diff(c(0, delayCDF))
  meanDelay = cumsum(thisTrialGapValues * delayPDF) / cumsum(delayPDF)
  rewardRate = tokenValue * delayCDF /
    ((meanDelay * delayCDF + thisTrialGapValues * (1 - delayCDF)) + iti)
  # save
  empOptimRewardRates[sIdx] = max(rewardRate)
  empOptimWaitTimes[sIdx] = thisTrialGapValues[which.max(rewardRate)]
  
  # plot
  if(isPlot){
    plotData = data.frame(time = thisTrialGapValues, emp = rewardRate, theo = rewardRates[[cond]])
    p = ggplot(plotData, aes(time, emp)) + geom_point() + geom_line(aes(time, theo)) +
      ggtitle(sprintf("%s, %s", sIdx, cond))
    print(p)
    readline("continue")
  }
}
plotData = data.frame(hdrData, empOptimRewardRates, empOptimWaitTimes = empOptimWaitTimes)
ggplot(plotData, aes(empOptimWaitTimes, fill = condition)) + geom_histogram() +
  scale_fill_manual(values = conditionColors) + saveTheme + xlab("waiting threshold / s")
ggsave("ideal_observer.png", width = 6, height = 4)


######## calculate waiting threshold using an asymptotical RL model
######## estimated from the empirical sequences
