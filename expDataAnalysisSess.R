# in this dataset, only trials within the 7 mins will be kept. Therefore, we don't need to delete any data
# determine whether to truncate data
isTrun = T
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
nBlock = 3
if(isTrun){
  tGrid = seq(0, blockSecs * nBlock, by = 1) # here I use a truncated tGrid, according to max(sellTime) 
}else{
  tGrid = seq(0, blockSecs * nBlock, by = 1)
}

# load all data
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
allIDs = hdrData$ID                   # column of subject IDs
n = length(allIDs)                    # n
cat('Analyzing data for',n,'subjects.\n')
# control which individual-level plots to generate
plotTrialwiseData =T
plotKMSC = F
plotWTW = F

# get session data 
tGrid = seq(0, blockSecs * nBlock, by = 1)
AUC = numeric(length = n)
totalEarnings =  numeric(length = n)
nAction = numeric(length = n)
wtwEarly = numeric(length = n)
nExclude = numeric(length = n)
timeWTW_ = vector(mode = "list", length = n)
trialWTW_ = vector(mode = "list", length = n)
kmOnGrid_ = vector(mode = "list", length = n)
winAUC_ = vector(mode = "list", length = n)
timeAUC_ = vector(mode = "list", length = n)
trialReRate_ = vector(mode = "list", length = n)
trialEndTime_ = vector(mode = "list", length = n)
stdQuitTime = numeric(length = n)
cvQuitTime = numeric(length = n)
muQuitTime = numeric(length = n)
nQuit = numeric(length = n)
nTrial = numeric(length = n)
stdWd = numeric(length =n)
cvWd =  numeric(length =n)
AUCEarly = numeric(length =n)
for (sIdx in 1 : n) {
  thisID = allIDs[sIdx]
  thisTrialData = trialData[[thisID]]
  if(unique(thisTrialData$stress == "no stress")){
    # truncate the last min(tMaxs) seconds
    if(isTrun){
      cond = unique(thisTrialData$condition)
      cIdx = ifelse(cond == "HP", 1, 2)
      excludedTrials = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[cIdx]))
      thisTrialData = thisTrialData[! (1 : nrow(thisTrialData) %in% excludedTrials),]
      nExclude[sIdx] = length(excludedTrials)
    }
    thisTrialData = block2session(thisTrialData)
    # generate arguments for later analysis 
    label = sprintf('Subject %s, Cond %s',thisID, unique(thisTrialData$condition))
    tMax = ifelse(unique(thisTrialData$condition) == conditions[1], tMaxs[1], tMaxs[2])
    
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
    kmscResults = kmsc(thisTrialData, min(tMaxs), label, plotKMSC, kmGrid)
    AUC[sIdx] = kmscResults[['auc']]
    kmOnGrid_[[sIdx]] = kmscResults$kmOnGrid
    if (plotKMSC) {
      readline(prompt = paste('subject',thisID, '(hit ENTER to continue)'))
      graphics.off()
    }
    stdWd[[sIdx]] = kmscResults$stdWd
    cvWd[[sIdx]] =  kmscResults$stdWd / kmscResults$auc
    
    # WTW time series
    wtwCeiling = min(tMaxs)
    wtwtsResults = wtwTS(thisTrialData, tGrid, wtwCeiling, label, plotWTW)
    timeWTW_[[sIdx]] = wtwtsResults$timeWTW
    trialWTW_[[sIdx]] = wtwtsResults$trialWTW
    wtwEarly[sIdx] =   wtwtsResults$trialWTW[1]
    if (plotWTW) {
      readline(prompt = paste('subject',thisID, '(hit ENTER to continue)'))
      graphics.off()
    }
    
    # moving auc
    # window = 10
    # by = 10
    # tempt = kmscMoving(thisTrialData, tMax, label, plotKMSC, tGrid, window, by)
    # timeAUC_[[sIdx]] = tempt$timeAUCs
    # winAUC_[[sIdx]] = tempt$winAUCs
    AUCEarly[sIdx] =  kmsc(truncateTrials(thisTrialData, 1, 10), min(tMaxs), label, plotKMSC, kmGrid)$auc
    
    # calculate rewardRates
    trialReRate_[[sIdx]] = trialEarnings / (timeWaited + iti)
    trialEndTime_[[sIdx]] = thisTrialData$sellTime    
  }
}
sessionData = data.frame(id = allIDs, condition = factor(hdrData$condition, levels = c("HP", "LP")), cbal = hdrData$cbal,
                         stress = factor(hdrData$stress, levels = c("no stress", "stress")), AUC = AUC, wtwEarly = wtwEarly,
                         totalEarnings = totalEarnings, nAction = nAction, stdQuitTime = stdQuitTime, cvQuitTime = cvQuitTime,
                         muQuitTime = muQuitTime, nQuit = nQuit, nTrial = nTrial,
                         stdWd = stdWd, cvWd = cvWd, AUCEarly = AUCEarly, nExclude = nExclude)
save(sessionData, file = 'genData/expDataAnalysis/sessionData.RData')
save(kmOnGrid_, file = 'genData/expDataAnalysis/kmOnGridSess.RData')


# plot AUC in two conditions
library("ggpubr")
sessionData[sessionData$stress == "no stress",]%>% ggplot(aes(condition, AUC)) + geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', aes(fill = condition)) +
  scale_fill_manual(values = conditionColors) + 
  xlab("") + ylab("Average WTW(s)") + myTheme +
  stat_compare_means(comparisons = list(c("HP", "LP")),
                     aes(label = ..p.signif..), label.x = 1.5, symnum.args= symnum.args,
                     bracket.size = 1, size = 6) + ylim(c(0, 23))
dir.create("figures/expDataAnalysisSess")
ggsave("figures/expDataAnalysisSess/zTruc_AUC.png", width = 4, height = 3)

# plot stdWd in two conditions
# sessionData[sessionData$stress == "no stress",]%>% ggplot(aes(condition, stdWd)) + geom_boxplot() +
#   geom_dotplot(binaxis='y', stackdir='center', aes(fill = condition)) +
#   scale_fill_manual(values = conditionColors) + 
#   xlab("") + ylab(expression(bold(paste("WTW S.D.","(", "s"^2, ")")))) + myTheme +
#   stat_compare_means(comparisons = list(c("HP", "LP")),
#                      aes(label = ..p.signif..), label.x = 1.5, symnum.args= symnum.args,
#                      bracket.size = 1, size = 6) + ylim(c(0, 24))
# dir.create("figures/expDataAnalysis")
# ggsave(sprintf("figures/expDataAnalysis/stdWD_%s.png", dataType), width = 4, height = 3.5)

# plot correlations 
# sessionData[sessionData$stress == "no stress",]%>% ggplot(aes(AUC, stdWd, color = condition)) + geom_point() +
#   facet_grid(~condition, scales = "free") + xlab("WTW Average (s)") +
#   ylab(expression(bold(paste("WTW S.D.","(", "s"^2, ")")))) +  scale_color_manual(values = conditionColors) +
#   myTheme
# dir.create("figures/expDataAnalysis")
# ggsave(sprintf("figures/expDataAnalysis/stdWD_AUC_%s.png", dataType), width = 6, height = 3.5)

# plot wtw 
select = (sessionData$stress == "no stress")
plotData = data.frame(wtw = unlist(timeWTW_[select]), time = rep(tGrid, sum(select)),
                      condition = rep(sessionData$condition[select], each = length(tGrid))) %>% group_by(condition, time) %>%
  dplyr::summarise(mean = mean(wtw), se = sd(wtw) / sqrt(length(wtw)), min = mean - se, max = mean + se) 

policy = data.frame(condition = c("HP", "LP"), wt = c(20, 2.2))
blockEnds = cumsum(c(blockSecs, blockSecs))
plotData %>% ggplot(aes(time, mean, color = condition)) +
  geom_ribbon(aes(ymin=min, ymax=max, fill = condition, colour=NA),alpha = 0.3) +
  geom_line(size = 1) + facet_wrap(~condition, scales = "free") +
  scale_color_manual(values = conditionColors) + scale_fill_manual(values = conditionColors) + xlab("Cumulative task time (min)") +
  scale_x_continuous(breaks = seq(0, max(tGrid), by = 60*7),
                     labels = seq(0, max(tGrid), by = 60*7) / 60) + 
  ylab("Willingness to wait (s)") +
  geom_hline(data = policy, aes(yintercept = wt, color = condition), linetype = "dotted", size = 1.5) +
  myTheme + ylim(c(0, 22))  +
  geom_vline(xintercept = blockEnds, color = "#969696", linetype = 2)
ggsave("figures/expDataAnalysisSess/zTruc_wtw_timecourse.png", width = 6, height = 3)

# plot survival curve
select = (sessionData$stress == "no stress")
condition =  sessionData$condition[select]

data.frame(kmsc = unlist(kmOnGrid_[select]),
           time = rep(kmGrid, sum(select)),
           condition = rep(condition, each = length(kmGrid))) %>% group_by(condition, time) %>%
  dplyr::summarise(mean = mean(kmsc), se = sd(kmsc) / sqrt(length(kmsc)), min = mean - se, max = mean + se) %>% 
  ggplot(aes(time, mean, color = condition, fill = condition)) + 
  geom_ribbon(aes(ymin=min, ymax=max),alpha = 0.3,  colour=NA)+
  geom_line(size = 1.5) + myTheme + scale_fill_manual(values = conditionColors) + 
  xlab("Elapsed time (s)") + ylab("Survival rate") + scale_color_manual(values = conditionColors)
ggsave("figures/expDataAnalysisSess/zTruc_kmsc_timecourse.png", width = 5, height = 4) 

# learning curve
trialReRareMove_ = lapply(1 : n, function(i){
  movAve(trialReRate_[[i]], 11)
})
timeReRate_ = lapply(1 :  n,
                     function(i) trial2sec(trialReRareMove_[[i]], trialEndTime_[[i]], tGrid))
select = (sessionData$stress == "no stress")
data.frame(value = unlist(timeReRate_[select]), time = rep(tGrid, sum(select)* nBlock),
           condition = factor(rep(sessionData$condition[select], each = length(tGrid)),  levels = conditions)) %>%
  group_by(condition, time) %>%
  dplyr::summarise(mean = mean(value), se = sd(value) / sqrt(length(value)), min = mean - se, max = mean + se) %>% 
  ggplot(aes(time, mean, color = condition, fill = condition)) + 
  geom_ribbon(aes(ymin=min, ymax=max), colour=NA, alpha = 0.3)+
  geom_line(size = 1.5) + myTheme + scale_fill_manual(values = conditionColors) +
  xlab("Elapsed time (s)") + ylab(expression(bold("Reward rate","(", "s"^2, ")"))) +
  scale_color_manual(values = conditionColors)
ggsave("figures/expDataAnalysisSess/zTruc_reRate.png", width = 4, height = 3.5) 

