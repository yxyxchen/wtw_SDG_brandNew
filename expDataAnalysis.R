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
plotTrialwiseData =T
plotKMSC = F
plotWTW = F

# initialize outputs, organised by block
tGrid = seq(0, blockSecs, by = 0.1)
AUC = numeric(length =n * nBlock)
totalEarnings =  numeric(length =n * nBlock)
nAction = numeric(length =n * nBlock)
wtwEarly = numeric(length =n * nBlock)
kmOnGrid_ = vector(mode = "list", length = n * nBlock)
varQuitTime = numeric(length =n * nBlock)
cvQuitTime = numeric(length =n * nBlock)
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
    # delete the NA row
    if(any(is.na(thisTrialData[nrow(thisTrialData),]))) thisTrialData = thisTrialData[1 : (nrow(thisTrialData) -1),]
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
    varQuitTime[noIdx] = ifelse(totalEarnings == 0, NA, var(timeWaited[trialEarnings == 0]))
    cvQuitTime[noIdx] = ifelse(totalEarnings == 0, NA, sd(timeWaited[trialEarnings == 0]) / mean(timeWaited[trialEarnings == 0]))
    # plot trial-by-trial data
    if (plotTrialwiseData) {
      trialPlots(thisTrialData,label)
    }
    
    # survival analysis
    kmscResults = kmsc(thisTrialData,tMax,label,plotKMSC,kmGrid)
    AUC[noIdx] = kmscResults[['auc']]
    kmOnGrid_[[noIdx]] = kmscResults$kmOnGrid

    # WTW time series
    wtwCeiling = tMax
    wtwtsResults = wtwTS(thisTrialData, tGrid, wtwCeiling, label, plotWTW)
    wtwEarly[noIdx] = mean(wtwtsResults[1 : (1 * 60 * 10)])
    
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
                       totalEarnings = totalEarnings, nAction = nAction, varQuitTime = varQuitTime, cvQuitTime = cvQuitTime, nTrial = nTrial)
save(blockData, file = 'genData/expDataAnalysis/blockData.RData')

# get session data 
AUC = numeric(length =n)
totalEarnings =  numeric(length =n)
kmOnGrid_ = vector(mode = "list", length = n)
for (sIdx in 1 : n) {
  thisID = allIDs[sIdx]
  # select data 
  thisTrialData = trialData[[thisID]]
  
  # generate arguments for later analysis 
  label = sprintf('Subject %s, Cond %s, Stress %s)',thisID, hdrData$condition[sIdx], hdrData$stress[sIdx])
  tMax = ifelse(hdrData$condition[sIdx] == conditions[1], tMaxs[1], tMaxs[2])
  kmGrid = seq(0, tMax, by= 0.1) # grid on which to average survival curves.
  
  totalEarnings[sIdx] =  sum(blockData$totalEarnings[blockData$id == thisID])
  
  # survival analysis
  kmscResults = kmsc(thisTrialData,tMax,label,plotKMSC,kmGrid)
  AUC[sIdx] = kmscResults[['auc']]
  kmOnGrid_[[sIdx]] = kmscResults$kmOnGrid
}
subData = data.frame(id = allIDs, condition = factor(hdrData$condition, levels = c("HP", "LP")),
                       stress = factor(hdrData$stress, levels = c("no stress", "stress")), AUC = AUC, 
                       totalEarnings = totalEarnings)
save(subData, file = 'genData/expDataAnalysis/subData.RData')
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