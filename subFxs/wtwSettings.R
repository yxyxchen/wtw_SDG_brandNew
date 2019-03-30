# not so good to generate varibale in scripts, since sometimes it invloves variables you don't wont
######## condition varibles #########
conditions = c("HP", "LP")
conditionColors = c("#008837", "#7b3294")

######## timing variables ########
tMaxs = c(20, 40) # trial durations
blockMins = 7 # block duration in mins
blockSecs = blockMins * 60 # block duration in secs
iti = 2 # iti duration in secs
tGrid = seq(0, blockSecs, stepDuration) # use stepDu 

######### reward variable ########
tokenValue = 10 #value of the token
loseValue = 0
stepDuration = 0.5
########## supporting vairbales ########
# time ticks within a trial for timeEarnings or wtw analysis
trialTicks = list(
  'HP' = round(seq(0, tMaxs[1], by = 0.1), 1),
  'LP' = round(seq(0, tMaxs[2], by = 0.1), 1)
)

########## additional  variables for optimal analysis ########
# CDF of reward delays: p(t_reward <= T)
k = 4
mu = 0
sigma = 2
pareto = list()
pareto[['k']] = k
pareto[['mu']] = mu
pareto[['sigma']] = sigma
HP = 1 / trialTicks$HP[length(trialTicks$HP)]  * trialTicks$HP
LP = 1 - (1 + k * (trialTicks$LP - mu) / sigma) ^ (-1 / k)
LP[length(trialTicks$LP)] = 1 
rewardDelayCDF = list(
  HP = HP,
  LP = LP
)

# library(ggplot2)
# source('plotThemes.R')
# plotData = data.frame(time = c(trialTicks$HP, trialTicks$LP),
#                       cdf = c(rewardDelayCDF$HP, rewardDelayCDF$LP),
#                       condition = c(rep('HP', length(trialTicks$HP)),
#                                     rep('LP', length(trialTicks$LP))))
# ggplot(plotData, aes(time, cdf, linetype = condition)) + geom_line() + xlab('Elapsed time / s') +
#   ylab('CDF') + ggtitle('Timing conditions') + saveTheme
# 
# ggsave('../outputs/exp_figures/timing_conditions.png', width = 3, height = 2)


#  PDF of reward delays: p(t_reward = T)
#  make it discrete 
HP = c(0, diff(rewardDelayCDF$HP))
LP = c(0, diff(rewardDelayCDF$LP))
rewardDelayPDF = list(
  "HP" = HP,
  "LP" = LP
)

# E(t_reward | t_reward <= T) 
HP = cumsum(trialTicks$HP * rewardDelayPDF$HP) / cumsum(rewardDelayPDF$HP)
HP[1] = NaN
LP = cumsum(trialTicks$LP * rewardDelayPDF$LP) / cumsum(rewardDelayPDF$LP)
LP[1] = NaN
# no reward arrives before the first reward timing, so points before that turn to NAN
meanRewardDelay = list('HP' = HP, 'LP' = LP)

# rewardRate
HP = tokenValue * rewardDelayCDF$HP /
  (meanRewardDelay$HP * rewardDelayCDF$HP + trialTicks$HP * (1 - rewardDelayCDF$HP) + iti)
LP = tokenValue * rewardDelayCDF$LP /
  (meanRewardDelay$LP * rewardDelayCDF$LP + trialTicks$LP * (1 - rewardDelayCDF$LP) + iti)
# quitting before the first reward timing get 0 reward
HP[which(is.nan(HP))] = 0 
LP[which(is.nan(LP))] = 0 
rewardRate = list('HP' = HP, 'LP' = LP)

optimWaitTimes = list()
optimWaitTimes$HP = trialTicks$HP[which.max(HP)]
optimWaitTimes$LP = trialTicks$LP[which.max(LP)]

optimRewardRates = list()
optimRewardRates$HP = max(HP)
optimRewardRates$LP = max(LP)

# calculate wInis 
wInisTheory = vector()
for(c in 1 : 2){
  cond = conditions[c];
  trialTick = trialTicks[[cond]]
  thisDelayPDF = rewardDelayPDF[[cond]]
  nTicks = length(trialTick)

  # assume gamma = 0.9
  gamma = 0.90
  r = - log(gamma) / stepDuration
  actionValueWaits = rep(0, nTicks)
  for(k in 1 : nTicks){
    actionValueWaits[k] = sum(tokenValue * exp(- (trialTick[k : nTicks] - trialTick[k]) * r)* thisDelayPDF[k : nTicks] / sum( thisDelayPDF[k : nTicks]))
  }
  junk = mean(actionValueWaits)
  wInisTheory[c] = junk
}
wInis = vector(length = 2)
wInis[1] = 3 # value of waiting, since participants didn't know differences in conditions 
wInis[2] = 2 # value of quitting, ensuring waiting first.

load("genData/expDataAnalysis/subData.RData")
load("genData/expDataAnalysis/blockData.RData")
modelName = "full_model"
source("subFxs/helpFxs.R")
source("subFxs/loadFxs.R")
paras = getParas(modelName)
expPara = loadExpPara(modelName, paras)

idList = unique(blockData$id)
useID = getUseID(blockData, expPara, paras)
# load expPara

wInisExp = vector(length = 2)
wInisExp[1] = median(expPara[subData$id %in% useID & subData$condition == "HP",4])
wInisExp[2] = median(expPara[subData$id %in% useID & subData$condition == "LP",4])


paraColors = list("phi" = "#78AB05","tau" = "#D9541A", "gamma" = "deepskyblue4", "QwaitIni" = "darkgoldenrod2") 

save("conditions", "conditionColors", "tMaxs", "blockMins", "blockSecs", "iti", "tGrid", 
     "tokenValue", "stepDuration", "trialTicks", "pareto", "rewardDelayCDF", 
     "rewardDelayPDF", "meanRewardDelay", "rewardRate", "optimRewardRates", 
     "optimWaitTimes", "wInis", "wInisTheory", "wInisExp", "paraColors", "loseValue", file = "wtwSettings.RData")