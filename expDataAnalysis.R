# in this dataset, only trials within the 7 mins will be kept. Therefore, we don't need to delete any data
# load libraries
source('subFxs/loadFxs.R') # for loading data 
source('subFxs/analysisFxs.R') # for analysis 
source("subFxs/plotThemes.R")
library("ggplot2")
library('dplyr')
dir.create("genData")
dir.create("genData/expDataAnalysis")

# load setting parameters 
load("wtwSettings.RData")

# load all data
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
allIDs = hdrData$ID                   # column of subject IDs
n = length(allIDs)                    # n
cat('Analyzing data for',n,'subjects.\n')

# define nBlock
nBlock = 3

# control which individual-level plots to generate
plotTrialwiseData =F
plotKMSC = F
plotWTW = F

# initialize outputs, organised by block
tGrid = seq(0, blockSecs, by = 0.1)
AUC = numeric(length =n * nBlock)
totalEarnings =  numeric(length =n * nBlock)
nAction = numeric(length =n * nBlock)
wtwEarly = numeric(length =n * nBlock)
timeWTW_ = vector(mode = "list", length = n * nBlock)
trialWTW_ = vector(mode = "list", length = n * nBlock)
kmOnGrid_ = vector(mode = "list", length = n * nBlock)
stdQuitTime = numeric(length =n * nBlock)
cvQuitTime = numeric(length =n * nBlock)
muQuitTime = numeric(length =n * nBlock)
nQuit = numeric(length =n * nBlock)
nTrial = numeric(length =n * nBlock)
# descriptive statistics for individual subjects and blocks
for (sIdx in 1 : n) {
  thisID = allIDs[sIdx]
  #if(blockData[blockData$id == thisID, "AUC"] > 20 & blockData$condition[blockData$id == thisID] == "LP"){
  for (bkIdx in 1: nBlock){
    # select data 
    thisTrialData = trialData[[thisID]]
    thisBlockIdx = (thisTrialData$blockNum == bkIdx)
    thisTrialData = thisTrialData[thisBlockIdx,]
    # generate arguments for later analysis 
    label = sprintf('Subject %s, Cond %s, %s',thisID, unique(thisTrialData$condition), hdrData$stress[sIdx])
    noIdx = (sIdx - 1) * nBlock + bkIdx # 
    tMax = ifelse(unique(thisTrialData$condition) == conditions[1], tMaxs[1], tMaxs[2])
    kmGrid = seq(0, tMax, by=0.1) # grid on which to average survival curves.
    
    # calcualte totalEarnings
    totalEarnings[noIdx] =  sum(thisTrialData$trialEarnings)
    timeWaited = thisTrialData$timeWaited
    trialEarnings = thisTrialData$trialEarnings
    scheduledWait = thisTrialData$scheduledWait
    timeWaited[trialEarnings > loseValue] = scheduledWait[trialEarnings > loseValue]
    nAction[noIdx] = sum(round(ifelse(trialEarnings > loseValue, ceiling(timeWaited / stepDuration), floor(timeWaited / stepDuration) + 1)))
    nTrial[noIdx] = length(timeWaited)
    # calculate varQuitTime
    stdQuitTime[noIdx] = ifelse(totalEarnings[noIdx] == 0, NA, sd(timeWaited[trialEarnings == 0]))
    cvQuitTime[noIdx] = ifelse(totalEarnings[noIdx] == 0, NA, sd(timeWaited[trialEarnings == 0]) / mean(timeWaited[trialEarnings == 0]))
    muQuitTime[noIdx] = mean(timeWaited[trialEarnings == 0])
    nQuit[noIdx] = sum(trialEarnings == 0)
      
    # plot trial-by-trial data
    if (plotTrialwiseData) {
      trialPlots(thisTrialData,label)
      readline(prompt = paste('subject',thisID, "block", bkIdx, '(hit ENTER to continue)'))
      graphics.off()
    }
    
    # survival analysis
    kmscResults = kmsc(thisTrialData,tMax,label,plotKMSC,kmGrid)
    AUC[noIdx] = kmscResults[['auc']]
    kmOnGrid_[[noIdx]] = kmscResults$kmOnGrid
    if (plotKMSC) {
      readline(prompt = paste('subject',thisID, "block", bkIdx, '(hit ENTER to continue)'))
      graphics.off()
    }

    # WTW time series
    wtwCeiling = tMax
    wtwtsResults = wtwTS(thisTrialData, tGrid, wtwCeiling, label, plotWTW)
    timeWTW_[[noIdx]] = wtwtsResults$timeWTW
    trialWTW_[[noIdx]] = wtwtsResults$trialWTW
    wtwEarly[noIdx] =   wtwtsResults$trialWTW[1]
    # wait for input before continuing, if individual plots were requested
    if (any(plotTrialwiseData, plotKMSC, plotWTW)) {
      readline(prompt = paste('subject',thisID, "block", bkIdx, '(hit ENTER to continue)'))
      graphics.off()
    }
  } # loop over blocks
}

# save data
save(kmOnGrid_, file = 'genData/expDataAnalysis/kmOnGridBlock.RData')
blockData = data.frame(id = rep(allIDs, each = nBlock), blockNum = rep( t(1 : nBlock), n),
                       cbal = rep(hdrData$cbal, each = nBlock), condition = factor(rep(hdrData$condition, each = nBlock), levels = c("HP", "LP")),
                       stress = factor(rep(hdrData$stress, each = nBlock), levels = c("no stress", "stress")), AUC = AUC, wtwEarly = wtwEarly,
                       totalEarnings = totalEarnings, nAction = nAction, stdQuitTime = stdQuitTime, cvQuitTime = cvQuitTime,
                       muQuitTime = muQuitTime, nQuit = nQuit, nTrial = nTrial)
save(blockData, file = 'genData/expDataAnalysis/blockData.RData')

# I think I would like to know the cv values
hist(varQuitTime[blockData$blockNum == 1 & blockData$condition == "HP"])

idListHP =  as.numeric(row.names(hdrData)[hdrData$condition == "LP"])
plotData = data.frame(wtw = unlist(lapply(1:60, function(i) timeWTW_[[idListHP[i]*3-2]])), id = rep(idListHP, each = length(tGrid)),
                      time = rep(1: length(tGrid), 60))
summaryData =  dplyr::summarise(group_by(plotData, time), mu = mean(wtw), std = sd(wtw))
plotData2 = data.frame(mu = summaryData$mu, std = summaryData$std, time = tGrid)
ggplot(plotData2, aes(time, mu)) + geom_line(group = 1)

# get session data 
tGrid = seq(0, blockSecs * nBlock, by = 0.1)
AUC = numeric(length = n)
totalEarnings =  numeric(length = n)
nAction = numeric(length = n)
wtwEarly = numeric(length = n)
timeWTW_ = vector(mode = "list", length = n)
trialWTW_ = vector(mode = "list", length = n)
kmOnGrid_ = vector(mode = "list", length = n)
stdQuitTime = numeric(length = n)
cvQuitTime = numeric(length = n)
muQuitTime = numeric(length = n)
nQuit = numeric(length = n)
nTrial = numeric(length = n)
plotTrialwiseData =F
plotKMSC = F
plotWTW = F
for (sIdx in 1 : n) {
  thisID = allIDs[sIdx]
  tempt = trialData[[thisID]]
  # change trialNum, sellTime, trialS, totalEarnings
  thisTrialData = trialData[[thisID]]
  nTrials = sapply(1:nBlock, function(i) sum(tempt$blockNum == i))
  thisTrialData = within(tempt, {trialNum = trialNum + rep(c(0, cumsum(nTrials)[1:2]), time = nTrials);
                                  sellTime = sellTime + rep((1:3-1) * blockSecs, time = nTrials);
                                  trialStartTime = trialStartTime + rep((1:3-1) * blockSecs, time = nTrials);
                                  totalEarnings = totalEarnings +  rep(c(0, totalEarnings[cumsum(nTrials)[1:2]]),
                                                                         time = nTrials)
                                  })
  # generate arguments for later analysis 
  label = sprintf('Subject %s, Cond %s',thisID, unique(thisTrialData$condition))
  tMax = ifelse(unique(thisTrialData$condition) == conditions[1], tMaxs[1], tMaxs[2])
  kmGrid = seq(0, tMax, by=0.1) # grid on which to average survival curves.
  
  # calcualte totalEarnings
  totalEarnings[sIdx] =  sum(thisTrialData$trialEarnings)
  timeWaited = thisTrialData$timeWaited
  trialEarnings = thisTrialData$trialEarnings
  scheduledWait = thisTrialData$scheduledWait
  timeWaited[trialEarnings > loseValue] = scheduledWait[trialEarnings > loseValue]
  nAction[sIdx] = sum(round(ifelse(trialEarnings > loseValue, ceiling(timeWaited / stepDuration), floor(timeWaited / stepDuration) + 1)))
  nTrial[sIdx] = length(timeWaited)
  # calculate varQuitTime
  stdQuitTime[sIdx] = ifelse(totalEarnings[sIdx] == 0, NA, sd(timeWaited[trialEarnings == 0]))
  cvQuitTime[sIdx] = ifelse(totalEarnings[sIdx] == 0, NA, sd(timeWaited[trialEarnings == 0]) / mean(timeWaited[trialEarnings == 0]))
  muQuitTime[sIdx] = mean(timeWaited[trialEarnings == 0])
  nQuit[sIdx] = sum(trialEarnings == 0)
  
  # plot trial-by-trial data
  if (plotTrialwiseData) {
    trialPlots(thisTrialData,label)
    readline(prompt = paste('subject',thisID, '(hit ENTER to continue)'))
    graphics.off()
  }
  
  # survival analysis
  kmscResults = kmsc(thisTrialData, tMax, label, plotKMSC, kmGrid)
  AUC[sIdx] = kmscResults[['auc']]
  kmOnGrid_[[sIdx]] = kmscResults$kmOnGrid
  if (plotKMSC) {
    readline(prompt = paste('subject',thisID, '(hit ENTER to continue)'))
    graphics.off()
  }
  
  # WTW time series
  wtwCeiling = tMax
  wtwtsResults = wtwTS(thisTrialData, tGrid, wtwCeiling, label, plotWTW)
  timeWTW_[[sIdx]] = wtwtsResults$timeWTW
  trialWTW_[[sIdx]] = wtwtsResults$trialWTW
  wtwEarly[sIdx] =   wtwtsResults$trialWTW[1]
}
sessionData = data.frame(id = allIDs, condition = factor(hdrData$condition, levels = c("HP", "LP")), cbal = hdrData$cbal,
                       stress = factor(hdrData$stress, levels = c("no stress", "stress")), AUC = AUC, wtwEarly = wtwEarly,
                     totalEarnings = totalEarnings, nAction = nAction, stdQuitTime = stdQuitTime, cvQuitTime = cvQuitTime,
                     muQuitTime = muQuitTime, nQuit = nQuit, nTrial = nTrial)
save(sessionData, file = 'genData/expDataAnalysis/sessionData.RData')
save(kmOnGrid_, file = 'genData/expDataAnalysis/kmOnGridSess.RData')


## correlation between AUC and triats 
summaryData = read.csv("data/SDGdataset.csv", header = T)
summaryData$AUCblock1 = blockData$AUC[blockData$blockNum == 1]
summaryData$condition = ifelse(summaryData$Task..1...unif..2...gp. == 1, "HP", "LP")
predictors = c("Barratt.Impulsiveness",
               "Intolerance.of.Uncertainty",
               "Trait.Anxiety..STAIT.", 
               "Delay.of.Gratification",
               "BDI")
predictorNames = c("Impulsiveness", 
                   "Intolerence of uncertainty",
                   "Trait anxiety",
                   "Delay of gratification",
                   "Depression")
nPredictor = length(predictors)
for(i in 1 : nPredictor){
    predictor = predictors[i]
    predictorName = predictorNames[i]
    # calculate correlations
    corHP = cor.test(summaryData[summaryData$condition == "HP", predictor], summaryData[summaryData$condition == "HP",]$AUCblock1, method = "spearman")
    corLP = cor.test(summaryData[summaryData$condition == "LP", predictor], summaryData[summaryData$condition == "LP",]$AUCblock1, method = "spearman")
    rhoHP = round(corHP$estimate, 3)
    rhoLP= round(corLP$estimate, 3)
    pHP = round(corHP$p.value, 3)
    pLP = round(corLP$p.value, 3)
    textData = data.frame(label = c(paste(round(rhoHP,2), "(p=", pHP, ")"), paste(round(rhoLP,2), "(p=", pLP, ")")),
                          condition = c("HP", "LP"))
    ggplot(summaryData, aes_string(predictor, "AUCblock1")) + geom_point() + facet_grid(~condition) + saveTheme + ylim(c(0, 40)) +
      xlab(predictorName) + ylab("AUC / min") + geom_text(data = textData,aes(x = -Inf,y = -Inf, label = label),
                                                          hjust   = -0.2, vjust   = -0.7,color = "blue",size = 4, fontface = 2) 
    fileName = sprintf("figures/expDataAnalysis/%s_AUC.pdf", predictorName)
    ggsave(fileName, width = 6, height = 3)
}

for(i in 1 : nPredictor){
  predictor = predictors[i]
  predictorName = predictorNames[i]
  # calculate correlations
  cor = cor.test(summaryData[, predictor], summaryData$AUCblock1, method = "spearman")
  rho = round(cor$estimate, 3)
  p = round(cor$p.value, 3)
  label = paste(round(rho,2), "(p=", p, ")")
  ggplot(summaryData, aes_string(predictor, "AUCblock1")) + geom_point() + saveTheme + ylim(c(0, 40)) +
    xlab(predictorName) + ylab("AUC / min") + annotate("text", label = label, y = 0,
                                                       x = (max(summaryData[,predictor], na.rm = T) + min(summaryData[,predictor], na.rm = T)) / 2,
                                                       color = "blue",size = 6, fontface = 2) 
  fileName = sprintf("figures/expDataAnalysis/sm_%s_AUC_sm.pdf", predictorName)
  ggsave(fileName, width = 4, height = 6)
}


# analysis the effect of sequence
cutMins = 0.5
earlyScheduleMean = vector(length = n)
for (sIdx in 1 : n) {
  thisID = allIDs[sIdx]
  bkIdx =1
  # select data 
  thisTrialData = trialData[[thisID]]
  thisBlockIdx = (thisTrialData$blockNum == bkIdx)
  thisTrialData = thisTrialData[thisBlockIdx,]
  earlyScheduleMean[sIdx] = mean(thisTrialData$scheduledWait[1:5])
}
blockData = blockData[blockData$blockNum == 1,]
blockData$earlyScheduleMean = earlyScheduleMean
ggplot(blockData, aes(earlyScheduleMean, AUC)) + geom_point() + facet_grid(~condition)
cor.test(blockData$AUC[blockData$condition == "LP" & blockData$AUC < 20],
         blockData$earlyScheduleMean[blockData$condition == "LP" & blockData$AUC < 20], method = "spearman")