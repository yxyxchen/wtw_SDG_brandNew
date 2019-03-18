# not so good to generate varibale in scripts, since sometimes it invloves variables you don't wont
######## condition varibles #########
conditions = c("HP", "LP")
conditionColors = c("#008837", "#7b3294")

######## timing variables ########
tMaxs = c(20, 40) # trial durations
blockMins = 7 # block duration in mins
blockSecs = blockMins * 60 # block duration in secs
iti = 2 # iti duration in secs
tGrid = seq(0, blockSecs, 0.1)

######### reward variable ########
tokenValue = 10 #value of the token
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

# calculate action value in the Q learning 
getMeanReward= function(gapIdx, quitAfter, gammaPerStep){
  gamma = gammaPerStep^2 # discount per 1s
  nGap = length(trialTicks$HP) - 1
  nTick =  length(trialTicks$HP) 
  tickIdx = gapIdx + 1
  quitGap = quitAfter / 0.1
  quitTick = quitGap  + 1
  if(quitTick < tickIdx){
    meanRewardDiscount = NaN
    meanRewardTime = NaN
  }else{
    meanRewardDiscount =sum(gamma ^ (seq(0, quitGap - gapIdx) * 0.1 + 0.1) * rewardDelayPDF$HP[tickIdx : quitTick]) /
      sum(rewardDelayPDF$HP[tickIdx : quitTick])
    meanRewardTime = sum((seq(0, quitGap - gapIdx) * 0.1 + 0.1) * rewardDelayPDF$HP[tickIdx : quitTick]) /
      sum(rewardDelayPDF$HP[tickIdx : quitTick])
  }
  meanReward = list("discount" = meanRewardDiscount,
                    "time" = meanRewardTime)
  return(meanReward)
}

getMeanWait = function(gapIdx, quitAfter, gammaPerStep){
  gamma = gammaPerStep^2 # discount per 1s
  nGap = length(trialTicks$HP) - 1
  nTick =  length(trialTicks$HP) 
  tickIdx = gapIdx + 1
  quitGap = quitAfter / 0.1
  quitTick = quitGap  + 1
  if(quitTick < tickIdx){
    meanWaitDiscount = NaN
    meanWaitTime = NaN
  }else{
    meanWaitDiscount =sum(gamma ^ pmin(seq(0, nGap - gapIdx) * 0.1 + 0.1, quitAfter - gapIdx * 0.1 + 0.1) * rewardDelayPDF$HP[tickIdx : nTick]) /
      sum(rewardDelayPDF$HP[tickIdx : nTick]) 
    meanWaitTime = sum(pmin(seq(0, nGap - gapIdx) * 0.1 + 0.1, quitAfter - gapIdx * 0.1 + 0.1) * rewardDelayPDF$HP[tickIdx : nTick]) /
      sum(rewardDelayPDF$HP[tickIdx : nTick])
  }
  meanWait = list("discount" = meanWaitDiscount,
                  "time" = meanWaitTime)
  return(meanWait)
}

getPreward = function(gapIdx, quitAfter, gammaPerStep){
  gamma = gammaPerStep^2 # discount per 1s
  nGap = length(trialTicks$HP) - 1
  nTick =  length(trialTicks$HP) 
  tickIdx = gapIdx + 1
  quitGap = quitAfter / 0.1
  quitTick = quitGap  + 1 
  pReward = (rewardDelayCDF$HP[quitTick] - rewardDelayCDF$HP[tickIdx - 1]) / (1 - rewardDelayCDF$HP[tickIdx - 1])
  return(pReward)
}



# calculate Qwait for inf trials
gammaPerStep = 0.99 # discount per 0.5s 
nGap = length(trialTicks$HP) - 1
nTick =  length(trialTicks$HP) 
Qquit_ = vector(length = nGap)
Qwait_ = vector(mode = "list", length = nTick)
Qwait_[[1]] = NaN
# t means quitting after t
for(quitTick in 1 : nTick){
  quitAfter = 0.1 * (quitTick-1)
  meanReward = getMeanReward (1, quitAfter, gammaPerStep)
  meanRewardDiscount =  meanReward$discount * gamma ^ iti
  meanWait = getMeanWait(1, quitAfter, gammaPerStep) 
  meanWaitDiscount =  meanWait$discount * gamma ^ iti
  pReward = getPreward(1, quitAfter, gammaPerStep)
  Qquit =  pReward * tokenValue * meanRewardDiscount / (1 - meanWaitDiscount)
  # Qquit =  sum(pReward * tokenValue * meanRewardDiscount * meanWaitDiscount ^ (0 : 10))
  Qquit_[quitTick] = Qquit
  if(quitTick > 1) {
    Qwait = vector(length = (quitTick - 1))
    for(j in 1 :  (quitTick - 1)){
      meanReward = getMeanReward(j, quitAfter, gammaPerStep)
      meanWait = getMeanWait(j, quitAfter, gammaPerStep)
      meanRewardDiscount = meanReward$discount
      meanWaitDiscount = meanWait$discount
      pReward = getPreward(j, quitAfter, gammaPerStep)
      Qwait[j] = tokenValue * meanRewardDiscount * pReward + Qquit * meanWaitDiscount
    }
    Qwait_[[quitTick]] = Qwait
  }
}
plot(Qquit_)

plotData = data.frame(Qwait = unlist(Qwait_), Qquit = rep(Qquit_, times = sapply(1 : nTick, function(i) length(Qwait_[[i]]))),
                      threshold = rep(trialTicks$HP, times = sapply(1 : nTick, function(i) length(Qwait_[[i]]))),
                      time = c(NaN, unlist(sapply(1 : nGap, function(i) seq(0.1, i * 0.1, by = 0.1)))))
library("ggplot2")
library("dplyr")
library("tidyr")
plotData = gather(plotData, action, value, -c("time", "threshold"))
ggplot(plotData[plotData$threshold == 1,], aes(time, value, color = action))+
  geom_point()

quitAfter = 20
nGap = quitAfter / 0.1
meanRewardTime_ = vector(length = nGap)
pReward_ = vector(length = nGap)
for(i in 1 : nGap){
  junk = getMeanReward(i, quitAfter, gammaPerStep)
  meanRewardTime_[i] = junk$time
  pReward_[i]= getPreward(i, quitAfter, gammaPerStep)
}
plot(tempt)
plot(pReward_)
# calculate action value in the R learning, given the policy is always wait 
# remember, here the average reward is for per 0.1s, namely per tick, not per s
Qq






# LP                                                       
tList = c(2, 4, 10)
nT = length(tList)
LP_ = vector(mode = "list", length = nT)
LPQquit_ = vector(mode = "list", length = nT)
for(h in 1: nT){
  t =  tList[h]
  n = length(trialTicks$LP)
  LP = unlist(lapply(1 : (t / 0.1), function(i){
    values = rep(-rewardRate$LP[t / 0.1 + 1] *  (t -  0.1 * i + iti), n - i)# when get no rewards, wait for vein for (t - 0.1 * i) gap duration
    values[1 : (t / 0.1 + 1 - i)] = tokenValue - (0.1 * (0 : (t / 0.1 - i)) + iti)* rewardRate$LP[t / 0.1 + 1]
    
    weights = rewardDelayPDF$LP[(i + 1) : n] # rewardDelays for i : (n-1) gap = (i + 1) : n tick
    sum(values * weights)  / sum(weights)
  })) / 0.1
  meanTimeOutIti = sum(pmin(trialTicks$LP, t) * rewardDelayPDF$LP)
  LPQuit = rep(LP[1] * (meanTimeOutIti + iti) / (meanTimeOutIti + iti*2), length(LP))
  LP_[[h]] = LP
  LPQquit_[[h]] = LPQuit
}
plotData = data.frame(Qwait = unlist(LP_), Qquit = unlist(LPQquit_), threshold = factor(rep(tList, times = unlist(lapply(LP_, length)))),
                      time = unlist(lapply(1 : nT, function(i )seq(0.1, tList[i], by = 0.1))))
plotData = gather(plotData, action, value, -c("time", "threshold"))
ggplot(plotData, aes(time, value,  color = action)) + geom_point() + facet_grid(~threshold)


# calculate rewardRate from another way.
for(t in 0.1 : 20){
  n = length(trialTicks$HP)
  statProbHP = unlist(lapply(1 : (n-1), function(i) {
    sum(rewardDelayPDF$HP[(i+1): n]) / sum(rewardDelayPDF$HP[2:n] * (1 : (n-1)))
  }))
  thisStatProbHP = statProbHP[1 : (t / 0.1)] / sum(statProbHP[1 : (t / 0.1)] )
  
  rHP = unlist(lapply(1 : (t / 0.1), function(i) {
    rewardDelayPDF$HP[i+1] / sum(rewardDelayPDF$HP[(i+1): n]) 
  }))
  
  meanTimeOutIti = sum(pmin(trialGapValues$HP, t) * rewardDelayPDF$HP)
  (sum(rHP * thisStatProbHP * tokenValue) * meanTimeOutIti * 10) / (meanTimeOutIti + 2)
}




# library(ggplot2)
# source('plotThemes.R')
# plotData = data.frame(time = c(trialTicks$HP, trialTicks$LP),
#                       rewardRate = c(rewardRate$HP, rewardRate$LP),
#                       condition = c(rep('HP', length(trialTicks$HP)), 
#                                     rep('LP', length(trialTicks$LP))))
# ggplot(plotData, aes(time, rewardRate, linetype = condition)) + geom_line() + xlab('Giving-up time /s') +
#   ylab('Total Earnings / cent') + ggtitle('Expected payoff functions') + saveTheme
# ggsave('../outputs/exp_figures/expected_payoff.png', width = 3, height = 2)


#### calculate action value #######
# suppose the strategy is always wait
# using specific discount
# assume rewards all arrives at the right side 
# r = 0.5
# for(c in 1 : 2){
#   cond = conditions[c]
#   trialTick = trialTicks[[cond]]
#   nTicks = length(trialTick)
#   thisDelayPDF = rewardDelayPDF[[cond]]
#   # 
#   actionValueWaits = rep(0, nTicks)
#   for(i in 1 : nTicks){
#     actionValueWaits[i] = sum(tokenValue * exp(- (trialTick[i : nTicks] - trialTick[i]) * r)* thisDelayPDF[i : nTicks] / sum( thisDelayPDF[i : nTicks]))    
#   }
#   if(c == 1) HP = actionValueWaits else LP = actionValueWaits
# }

# ### calculate action values   
# # # suppose waiting to the last seconds
# condIdx = 2
# cond = conditions[[condIdx]]
# nTimePoint = tMaxs[[condIdx]] / 0.5
# thisRewardDelayPDF = rewardDelayPDF[[cond]]
# sparseRewardDelayPDF = c(0, rowSums(matrix(thisRewardDelayPDF[2:length(thisRewardDelayPDF)],
#                                       nrow = (length(thisRewardDelayPDF) - 1)/ 5)))
# 
# waitRate = 0.90
#   
# dv = matrix(NA, nTimePoint, nGamma )
# vaWaits = matrix(NA, nTimePoint, nGamma)
# vaQuits = matrix(NA, nTimePoint, nGamma )
# nGamma = 5
# gammaList = seq(0.60, 0.98, length.out = 5) # gamma for 0.5s
# for(h in 1 : nGamma){
#   gamma = gammaList[h]
#   Qwait = vector(length = nTimePoint)
#   for(i in 1 : nTimePoint){
#     discount = gamma ^ (0 : (nTimePoint - i)) 
#     waitRateSeq = waitRate^ (1 : (nTimePoint - i + 1))
#     Qwait[i] = sum(tokenValue * discount * thisRewardDelayPDF[(i + 1) : (nTimePoint+ 1)] * waitRateSeq)/
#       sum(thisRewardDelayPDF[(i + 1) : (nTimePoint + 1)])
#   }
#   Qquit = rep(Qwait[1] * gamma ^ 4, nTimePoint) * waitRate  + (1 - waitRate) * 0
#   vaWaits[,h] = Qwait
#   vaQuits[,h]= Qquit
#   dv[,h] = Qwait - Qquit
# }
# 
# tempt = data.frame(dv = dv, time = (1 : nTimePoint) * 0.1, Qwait = vaWaits, Qquit = vaQuits)
# plotData = data.frame(dv = as.vector(dv), Qwait = as.vector(vaWaits),
#                       Qquit = as.vector(vaQuits), gamma = as.factor(rep(gammaList, each = nTimePoint)),
#                       time = rep(1 : nTimePoint * 0.1, length(gammaList)))
# ggplot(plotData, aes(time, dv, color = gamma)) + geom_line(size = 1) +
#   ylab('Decision variable') + saveTheme+ggtitle(cond) + xlab('Time / s')

# library('dplyr')
# library('tidyr')
# plotData2 = gather(plotData, action, actionValue, -c(1,4,5))
# plotData2$action = ifelse(plotData2$action == 'Qwait', 'wait', 'quit')
# ggplot(plotData2, aes(time, actionValue, color = gamma, linetype= action)) + geom_line(size = 1) + saveTheme +ggtitle(cond) +
#   ylab('Action value') + xlab('Time / s')

load("genData/expDataAnalysis/subData.RData")
modelName = "full_model"
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
     "optimWaitTimes", "wInis", "wInisTheory", "wInisExp", "paraColors", file = "wtwSettings.RData")