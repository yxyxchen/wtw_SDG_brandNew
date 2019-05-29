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
stdWd = numeric(length =n * nBlock) # standard deviation from the survival curve for the whole block
cvWd =  numeric(length =n * nBlock)
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
    stdWd[noIdx] = kmscResults$stdWd
    cvWd[noIdx] = kmscResults$stdWd / kmscResults$auc
    
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
    if (plotWTW) {
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
                       muQuitTime = muQuitTime, nQuit = nQuit, nTrial = nTrial, stdWd = stdWd, cvWd = cvWd)
save(blockData, file = 'genData/expDataAnalysis/blockData.RData')



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
stdWd = numeric(length =n)
cvWd =  numeric(length =n)
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
  kmGrid = seq(0, tMax, by= 0.1) # grid on which to average survival curves.
  
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
  stdWd[[sIdx]] = kmscResults$stdWd
  cvWd[[sIdx]] =  kmscResults$stdWd / kmscResults$auc
  
  # WTW time series
  wtwCeiling = tMax
  wtwtsResults = wtwTS(thisTrialData, tGrid, wtwCeiling, label, plotWTW)
  timeWTW_[[sIdx]] = wtwtsResults$timeWTW
  trialWTW_[[sIdx]] = wtwtsResults$trialWTW
  wtwEarly[sIdx] =   wtwtsResults$trialWTW[1]
  if (plotWTW) {
    readline(prompt = paste('subject',thisID, '(hit ENTER to continue)'))
    graphics.off()
  }
  
  # # moving auc
  # window = 20
  # by = 10
  # tempt = kmscMoving(thisTrialData, tMax, label, plotKMSC, tGrid, window, by)
  # timeAUC_[[sIdx]] = tempt$timeAUCs
  # winAUC_[[sIdx]] = tempt$winAUCs
  # timeStdWd_[[sIdx]] = tempt$timeStdWds
}
sessionData = data.frame(id = allIDs, condition = factor(hdrData$condition, levels = c("HP", "LP")), cbal = hdrData$cbal,
                       stress = factor(hdrData$stress, levels = c("no stress", "stress")), AUC = AUC, wtwEarly = wtwEarly,
                     totalEarnings = totalEarnings, nAction = nAction, stdQuitTime = stdQuitTime, cvQuitTime = cvQuitTime,
                     muQuitTime = muQuitTime, nQuit = nQuit, nTrial = nTrial,
                     stdWd = stdWd, cvWd = cvWd)
save(sessionData, file = 'genData/expDataAnalysis/sessionData.RData')
save(kmOnGrid_, file = 'genData/expDataAnalysis/kmOnGridSess.RData')


# correlation between AUC and triats 
# determine summaryData
dataType = "sess"
if(dataType == "block"){
  load("genData/expDataAnalysis/blockData.RData")
  load("genData/expDataAnalysis/kmOnGridBlock.RData")
  summaryData = blockData[blockData$blockNum == 1, ] # so summaryData only have something relevant
}else{
  load("genData/expDataAnalysis/sessionData.RData")
  load("genData/expDataAnalysis/kmOnGridSess.RData")
  summaryData = sessionData
}

# load personality data
library("Hmisc")
personality = read.csv("data/SDGdataset.csv")
personality$id = personality$SubjectID
traits = c("Delay.of.Gratification", "Barratt.Impulsiveness","Intolerance.of.Uncertainty", "Trait.Anxiety..STAIT.")
traitNames = c("DG", "IMP", "UC", "AX")
nTrait = length(traits)
traitAUCCorr = list()
# plot separately for two conditions
for(i in 1 : nTrait){
  trait = traits[i];
  traitName = traitNames[i]
  input = data.frame(personality[summaryData$stress =="no stress",trait],
                     summaryData$AUC[summaryData$stress == "no stress"],
                     summaryData$condition[summaryData$stress == "no stress"])
  traitAUCCorr[[i]]= getCorrelation(input)
  p = plotCorrelation(input, isRank = T) 
  p + xlab(paste(capitalize(traitName), "(rank)")) + ylab("AUC (rank)") + myTheme
  fileName = sprintf("%s/AUC_%s_%s.png", "figures/expDataAnalysis", traitName, dataType)
  ggsave(fileName, width = 6, height = 3)
}
rhoTable = lapply(1:2, function(j) sapply(1: (nTrait), function(i) traitAUCCorr[[i]]$rhos[j]))
pTable = lapply(1:2, function(j) sapply(1: (nTrait), function(i) traitAUCCorr[[i]]$ps[j]))

# plot AUC in two conditions
library("ggpubr")
load("wtwSettings.RData")
summaryData[summaryData$stress == "no stress",]%>% ggplot(aes(condition, AUC)) + geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', aes(fill = condition)) +
  scale_fill_manual(values = conditionColors) + 
  xlab("") + ylab("AUC (s)") + myTheme +
  stat_compare_means(comparisons = list(c("HP", "LP")),
                     aes(label = ..p.signif..), label.x = 1.5, symnum.args= symnum.args,
                     bracket.size = 1, size = 6) + ylim(c(0, 47))
dir.create("figures/expDataAnalysis")
ggsave(sprintf("figures/expDataAnalysis/AUC_%s.png", dataType), width = 4, height = 3.5)

# plot wtw 
select = (sessionData$stress == "no stress")
plotData = data.frame(wtw = unlist(timeWTW_[select]), time = rep(tGrid, sum(select)),
           condition = rep(summaryData$condition[select], each = length(tGrid))) %>% group_by(condition, time) %>%
  summarise(mean = mean(wtw), se = sd(wtw) / sqrt(30), min = mean - se, max = mean + se) 

policy = data.frame(condition = c("HP", "LP"), wt = c(20, 2.2))
plotData %>% ggplot(aes(time, mean, color = condition)) +
  geom_ribbon(aes(ymin=min, ymax=max, fill = condition, colour=NA),alpha = 0.3) +
  geom_line(size = 1) + facet_wrap(~condition, scales = "free") +
  scale_color_manual(values = conditionColors) + scale_fill_manual(values = conditionColors) +
  ylim(c(3, 26)) + xlab("Cumulative task time (min)") +
  scale_x_continuous(breaks = seq(0, max(tGrid), by = 60*7),
                     labels = paste(seq(0, 21, by = 7))) + 
  ylab("Willingness to wait (s)") +
  geom_hline(data = policy, aes(yintercept = wt), linetype = "dotted", color = "#969696", size = 1.5) +
  myTheme + ylim(c(0, 25))
ggsave("figures/expDataAnalysis/wtw_timecourse.png", width = 6, height = 3)

# plot survival curve
select = (summaryData$stress == "no stress")
condition =  rep(summaryData$condition[select])
step = 0.1
len= sapply(1 : 60, function(i) ifelse(condition[i] == "HP", tMaxs[1] / step + 1, tMaxs[2] / step + 1))
ideal = data.frame(kmsc = c(rep(1, tMaxs[1] / step + 1),
                            c(rep(1, 2.2 / step + 1), rep(0, tMaxs[2] / step -  2.2 / step))),
                   time = c(seq(0, tMaxs[1], by = step), seq(0, tMaxs[2], by = step)),
                   condition = rep(c("HP", "LP"), time = (tMaxs / step) + 1))

data.frame(kmsc = unlist(kmOnGrid_[select]),
                      time = unlist(lapply(1 : 60, function(i) (1 : len[i]) * step)),
                      condition = rep(condition, time = len)) %>% group_by(condition, time) %>%
  summarise(mean = mean(kmsc), se = sd(kmsc) / sqrt(30), min = mean - se, max = mean + se) %>% 
  ggplot(aes(time, mean, color = condition)) + 
  geom_ribbon(aes(ymin=min, ymax=max, fill = condition, colour=NA),alpha = 0.3)+
  geom_line(size = 1.5) + myTheme + scale_fill_manual(values = conditionColors) + 
  xlab("Elapsed time (s)") + ylab("Survival rate") + scale_color_manual(values = conditionColors)
# +
#   geom_line(data = ideal, aes(time, kmsc, color = condition), linetype = "dotted", size = 1)
ggsave("figures/expDataAnalysis/kmsc_timecourse.png", width = 5, height = 4) 

