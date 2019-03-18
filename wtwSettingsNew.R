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
stepDuration = 0.1
########## supporting vairbales ########
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


########## additional  variables for optimal analysis ########
# CDF of reward delays: p(t_reward <= T)
k = 4
mu = 0
sigma = 2
pareto = list()
pareto[['k']] = k
pareto[['mu']] = mu
pareto[['sigma']] = sigma
HP = 1 / length(trialGapIdxs$HP) * trialGapIdxs$HP
LP = 1 - (1 + k * (trialGapValues$LP - mu) / sigma) ^ (-1 / k)
LP[length(trialGapValues$LP)] = 1 
rewardDelayCDF = list(
  HP = HP,
  LP = LP
)

#  PDF of reward delays: p(t_reward = T)
#  there is no gap before the first tick, therefore the first element is NaN
#  remmember to use the middle value to the gap when calculating the mean
HP = diff(c(0, rewardDelayCDF$HP))
LP = diff(c(0, rewardDelayCDF$LP))
rewardDelayPDF = list(
  "HP" = HP,
  "LP" = LP
)

# E(t_reward | t_reward <= T) 
HP = cumsum((trialGapValues$HP - stepDuration / 2) * rewardDelayPDF$HP) / cumsum(rewardDelayPDF$HP)
LP = cumsum((trialGapValues$LP - stepDuration / 2) * rewardDelayPDF$LP) / cumsum(rewardDelayPDF$LP)
# no reward arrives before the first reward timing, so points before that turn to NAN
meanRewardDelay = list('HP' = HP, 'LP' = LP)

# rewardRate
HP = tokenValue * rewardDelayCDF$HP /
  (meanRewardDelay$HP * rewardDelayCDF$HP + trialGapValues$HP * (1 - rewardDelayCDF$HP) + iti)
LP = tokenValue * rewardDelayCDF$LP /
  (meanRewardDelay$LP * rewardDelayCDF$LP + trialGapValues$LP * (1 - rewardDelayCDF$LP) + iti)
rewardRate = list('HP' = HP, 'LP' = LP)

optimWaitTimes = list()
optimWaitTimes$HP = trialTicks$HP[which.max(HP)]
optimWaitTimes$LP = trialTicks$LP[which.max(LP)]

optimRewardRates = list()
optimRewardRates$HP = max(HP)
optimRewardRates$LP = max(LP)

# calculate action values in Qlearning
getMeanRewardDelay = function(gapIdx, quitAfter, gammaPerStep){
  gamma = gammaPerStep^2 # discount per 1s
  nGap = length(trialTicks$HP) - 1
  nTick =  length(trialTicks$HP) 
  tickIdx = gapIdx + 1
  quitGap = quitAfter / stepDuration
  quitTick = quitGap  + 1
  if(quitTick < tickIdx){
    discount = NaN
    discount = NaN
  }else{
    discount =sum(gamma ^ (seq(0, quitGap - gapIdx) * stepDuration + stepDuration / 2) * rewardDelayPDF$HP[gapIdx : quitGap]) /
      sum(rewardDelayPDF$HP[gapIdx : quitGap])
    time = sum((seq(0, quitGap - gapIdx) * stepDuration + stepDuration / 2) * rewardDelayPDF$HP[gapIdx : quitGap]) /
      sum(rewardDelayPDF$HP[gapIdx : quitGap])
  }
  meanRewardDelay = list(discount = discount, time = time)
  return(meanRewardDelay)
}

getMeanWaitDelay = function(gapIdx, quitAfter, gammaPerStep){
  gamma = gammaPerStep^2 # discount per 1s
  nGap = length(trialTicks$HP) - 1
  nTick =  length(trialTicks$HP) 
  tickIdx = gapIdx + 1
  quitGap = quitAfter / stepDuration
  quitTick = quitGap  + 1
  if(quitTick < tickIdx){
    discount = NaN
    time = NaN
  }else{
    discount =sum(gamma ^ pmin(seq(0, nGap - gapIdx) * stepDuration + stepDuration / 2, quitAfter - gapIdx * stepDuration + stepDuration)
                  * rewardDelayPDF$HP[gapIdx : nGap]) /
      sum(rewardDelayPDF$HP[gapIdx : nGap]) 
    time = sum(pmin(seq(0, nGap - gapIdx) * stepDuration + stepDuration / 2, quitAfter - gapIdx * stepDuration + stepDuration) * rewardDelayPDF$HP[gapIdx : nGap]) /
      sum(rewardDelayPDF$HP[gapIdx : nGap])
  }
  meanWaitDelay = list(discount = discount, time = time)
  return(meanWaitDelay)
}

getPreward = function(gapIdx, quitAfter, gammaPerStep){
  gamma = gammaPerStep^2 # discount per 1s
  nGap = length(trialTicks$HP) - 1
  nTick =  length(trialTicks$HP) 
  tickIdx = gapIdx + 1
  quitGap = quitAfter / stepDuration
  if(gapIdx == 1){
    pReward = rewardDelayCDF$HP[quitGap]  
  }else{
    pReward = (rewardDelayCDF$HP[quitGap] - rewardDelayCDF$HP[gapIdx - 1]) / (1 - rewardDelayCDF$HP[gapIdx - 1])
  }
  return(pReward)
}

gammaPerStep = 0.90 # discount per 0.5s 
nGap = length(trialGapValues$HP) 
nTick =  length(trialTicks$HP) 
Qquit_ = vector(length = nGap)
Qwait_ = vector(mode = "list", length = nGap)
# t means quitting after t
for(quitGap in 1 : nGap){
  quitAfter = stepDuration * quitGap
  meanRewardDelay = getMeanRewardDelay(1, quitAfter, gammaPerStep)
  rewardedDiscount = meanRewardDelay$discount * gamma ^ iti
  meanWaitDelay = getMeanWaitDelay(1, quitAfter, gammaPerStep) 
  discount =  meanWaitDelay$discount * gamma ^ iti
  pReward = getPreward(1, quitAfter, gammaPerStep)
  Qquit =  pReward * tokenValue * rewardedDiscount / (1 - discount)
  # Qquit =  sum(pReward * tokenValue * meanRewardDiscount * meanWaitDiscount ^ (0 : 10))
  Qquit_[quitGap] = Qquit
  Qwait = vector(length = quitGap)
  for(j in 1 : quitGap){
    meanRewardDelay = getMeanRewardDelay(j, quitAfter, gammaPerStep)
    meanWaitDelay = getMeanWaitDelay(j, quitAfter, gammaPerStep)
    rewardedDiscount = meanRewardDelay$discount
    discount = meanWaitDelay$discount
    pReward = getPreward(j, quitAfter, gammaPerStep)
    Qwait[j] = tokenValue * rewardedDiscount * pReward + Qquit * discount
  }
  Qwait_[[quitGap]] = Qwait
}
plot(Qquit_)

plotData = data.frame(Qwait = unlist(Qwait_), Qquit = rep(Qquit_, times = sapply(1 : nGap, function(i) length(Qwait_[[i]]))),
                      threshold = rep(trialGapValues$HP, times = sapply(1 : nGap, function(i) length(Qwait_[[i]]))),
                      time = c(unlist(sapply(1 : nGap, function(i) seq(stepDuration, i * stepDuration, by = stepDuration)))))
library("ggplot2")
library("dplyr")
library("tidyr")
plotData = gather(plotData, action, value, -c("time", "threshold"))
ggplot(plotData[plotData$threshold %in% c(2, 4, 10, 20),], aes(time, value, color = action))+
  geom_point(size = 0.5) + facet_grid(~threshold)


##
Rt = vector(length = length(trialGapValues$HP))
for(i in 1 : length(trialGapValues$HP)){
  t = i * stepDuration
  n = length(trialGapValues$HP)
  statProbHP = unlist(lapply(1 : n, function(i) {
    sum(rewardDelayPDF$HP[i : n]) / sum(rewardDelayPDF$HP[1:n] * (1 : (n)))
  }))
  thisStatProbHP = statProbHP[1 : (t / 0.1)] / sum(statProbHP[1 : (t / 0.1)] )
  
  rHP = unlist(lapply(1 : (t / 0.1), function(i) {
    rewardDelayPDF$HP[i] / sum(rewardDelayPDF$HP[i: n]) 
  }))
  meanWaitDelay = getMeanWaitDelay(1, t, gammaPerStep)
  # meanTimeOutIti = sum(pmin(trialGapValues$HP, t) * rewardDelayPDF$HP)
  meanTimeOutIti = meanWaitDelay$time
  Rt[i] = (sum(rHP * thisStatProbHP * tokenValue) * meanTimeOutIti * 10) / (meanTimeOutIti + 2)
}

rewardRate$HP[[200]] / (1 - gamma)


sum(thisStatProbHP *meanTimeOutIti / (meanTimeOutIti + 2) * Qwait_[[200]]) +
    mean(gamma ^ -seq(0, 1.9, by = 0.1)) * Qquit_[[200]] * 2 / (meanTimeOutIti + 2)


Rt = vector(length = length(trialGapValues$HP))
for(i in 1 : length(trialGapValues$HP)){
  t = i * stepDuration
  n = length(trialGapValues$HP)
  statProbHP = unlist(lapply(1 : n, function(i) {
    sum(rewardDelayPDF$HP[i : n]) / sum(rewardDelayPDF$HP[1:n] )
  }))
  thisStatProbHP = statProbHP[1 : (t / 0.1)] / sum(statProbHP[1 : (t / 0.1)] )
  
  rHP = unlist(lapply(1 : (t / 0.1), function(i) {
    rewardDelayPDF$HP[i] / sum(rewardDelayPDF$HP[i: n]) 
  }))
  meanWaitDelay = getMeanWaitDelay(1, t, gammaPerStep)
  # meanTimeOutIti = sum(pmin(trialGapValues$HP, t) * rewardDelayPDF$HP)
  meanTimeOutIti = meanWaitDelay$time
  Rt[i] = (sum(rHP * thisStatProbHP * tokenValue) * meanTimeOutIti * 10) / (meanTimeOutIti + 2)
}

plotData = data.frame(Qquit = Qquit_, time = trialGapValues$LP)
ggplot(plotData, aes(time, Qquit)) + geom_point()









