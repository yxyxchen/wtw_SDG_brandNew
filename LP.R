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

# assume the rewards happen at the end of the gap
# E(t_reward | t_reward <= T) 
HP = cumsum((trialGapValues$HP) * rewardDelayPDF$HP) / cumsum(rewardDelayPDF$HP)
LP = cumsum((trialGapValues$LP) * rewardDelayPDF$LP) / cumsum(rewardDelayPDF$LP)
# no reward arrives before the first reward timing, so points before that turn to NAN
meanRewardDelay = list('HP' = HP, 'LP' = LP)

# rewardRate
HP = tokenValue * rewardDelayCDF$HP /
  ((meanRewardDelay$HP * rewardDelayCDF$HP + trialGapValues$HP * (1 - rewardDelayCDF$HP)) + iti)
LP = tokenValue * rewardDelayCDF$LP /
  ((meanRewardDelay$LP * rewardDelayCDF$LP + trialGapValues$LP * (1 - rewardDelayCDF$LP)) + iti)
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
  nGap = length(trialTicks$LP) - 1
  nTick =  length(trialTicks$LP) 
  quitGap = round(quitAfter / stepDuration)
  if(quitGap < gapIdx){
    discount = NaN
    discount = NaN
  }else{
    discount =sum(gamma ^ (seq(0, quitGap - gapIdx) * stepDuration + stepDuration) * rewardDelayPDF$LP[gapIdx : quitGap]) /
      sum(rewardDelayPDF$LP[gapIdx : quitGap])
    time = sum((seq(0, quitGap - gapIdx) * stepDuration + stepDuration) * rewardDelayPDF$LP[gapIdx : quitGap]) /
      sum(rewardDelayPDF$LP[gapIdx : quitGap])
  }
  meanRewardDelay = list(discount = discount, time = time)
  return(meanRewardDelay)
}

getMeanWaitDelay = function(gapIdx, quitAfter, gammaPerStep){
  gamma = gammaPerStep^2 # discount per 1s
  nGap = length(trialTicks$LP) - 1
  nTick =  length(trialTicks$LP) 
  tickIdx = gapIdx + 1
  quitGap = round(quitAfter / stepDuration)
  quitTick = quitGap  + 1
  if(quitTick < tickIdx){
    discount = NaN
    time = NaN
  }else{
    discount =sum(gamma ^ pmin(seq(0, nGap - gapIdx) * stepDuration + stepDuration, quitAfter - gapIdx * stepDuration + stepDuration)
                  * rewardDelayPDF$LP[gapIdx : nGap]) /
      sum(rewardDelayPDF$LP[gapIdx : nGap]) 
    time = sum(pmin(seq(0, nGap - gapIdx) * stepDuration + stepDuration, quitAfter - gapIdx * stepDuration + stepDuration) * rewardDelayPDF$LP[gapIdx : nGap]) /
      sum(rewardDelayPDF$LP[gapIdx : nGap])
  }
  meanWaitDelay = list(discount = discount, time = time)
  return(meanWaitDelay)
}

getPreward = function(gapIdx, quitAfter, gammaPerStep){
  gamma = gammaPerStep^2 # discount per 1s
  nGap = length(trialTicks$LP) - 1
  nTick =  length(trialTicks$LP) 
  tickIdx = gapIdx + 1
  quitGap = round(quitAfter / stepDuration)
  if(gapIdx == 1){
    pReward = rewardDelayCDF$LP[quitGap]  
  }else{
    pReward = (rewardDelayCDF$LP[quitGap] - rewardDelayCDF$LP[gapIdx - 1]) / (1 - rewardDelayCDF$LP[gapIdx - 1])
  }
  return(pReward)
}

gammaPerStep = 0.90 # discount per 0.5s 
gamma = gammaPerStep ^ 2
nGap = length(trialGapValues$LP) 
nTick =  length(trialTicks$LP) 
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
library("ggplot2")
library("dplyr")
library("tidyr")
source("subFxs/plotThemes.R")

# ok there are the same
nGap = length(trialGapValues$LP)
QLP = vector(length = length(nGap))
for(lastWaitGap in 1 : nGap){
  # define the transition matrix
  trans = matrix(rep(0, length = (lastWaitGap + nQuitGap)^2), nrow = lastWaitGap + nQuitGap)
  if(lastWaitGap > 1) trans[1 : (lastWaitGap - 1), lastWaitGap + 1]  =
      unlist(lapply(1:(lastWaitGap-1), function(i) rewardDelayPDF$LP[i] / sum(rewardDelayPDF$LP[i : nGap])))
  trans[lastWaitGap, lastWaitGap + 1] = 1
  if(lastWaitGap > 1) for(i in 1 : (lastWaitGap - 1)) trans[i, i+1] = 1 - trans[i, lastWaitGap + 1]
  if(nQuitGap > 1) for(i in (lastWaitGap+1) : (lastWaitGap + nQuitGap-1)) trans[i, i + 1] = 1
  trans[lastWaitGap + nQuitGap, 1] = 1
  
  # calculate the steady state
  n = ncol(trans)
  A = t(trans - diag(n))
  A = rbind(A, rep(1, n))
  b = c(rep(0, n), 1)
  ULP = qr.solve(A, b)
  QquitList = Qquit_[[lastWaitGap]] * gamma ^ (-(0 : (nQuitGap - 1)) * stepDuration)
  QLP[[lastWaitGap]] =  sum(ULP[1:lastWaitGap] * Qwait_[[lastWaitGap]]) + sum(ULP[lastWaitGap + (1 : nQuitGap)] * QquitList)
}


nQuitGap = 20
Rt = vector(length = length(trialGapValues$LP))
Rstep= vector(length = length(trialGapValues$LP))
for(lastWaitGap in 1 : nGap){
  # define the transition matrix
  trans = matrix(rep(0, length = (lastWaitGap + nQuitGap)^2), nrow = lastWaitGap + nQuitGap)
  if(lastWaitGap > 1) trans[1 : (lastWaitGap - 1), lastWaitGap + 1]  =
      unlist(lapply(1:(lastWaitGap-1), function(i) rewardDelayPDF$LP[i] / sum(rewardDelayPDF$LP[i : nGap])))
  trans[lastWaitGap, lastWaitGap + 1] = 1
  if(lastWaitGap > 1) for(i in 1 : (lastWaitGap - 1)) trans[i, i+1] = 1 - trans[i, lastWaitGap + 1]
  if(nQuitGap > 1) for(i in (lastWaitGap+1) : (lastWaitGap + nQuitGap-1)) trans[i, i + 1] = 1
  trans[lastWaitGap + nQuitGap, 1] = 1
  # calculate the steady state
  n = ncol(trans)
  A = t(trans - diag(n))
  A = rbind(A, rep(1, n))
  b = c(rep(0, n), 1)
  ULP = qr.solve(A, b)
  rLP = c(unlist(lapply(1 : lastWaitGap, function(i) {
    rewardDelayPDF$LP[i] / sum(rewardDelayPDF$LP[i: nGap]) 
  })), rep(0, nQuitGap))
  timeInGaps = c(unlist(lapply(1 : lastWaitGap, function(i) {
    junk = rewardDelayPDF$LP[i] / sum(rewardDelayPDF$LP[i : nGap])
    junk * stepDuration + (1 - junk) * stepDuration
  })), rep(iti / nQuitGap, nQuitGap))
  meanTimePerAction = sum(ULP * timeInGaps)
  Rt[lastWaitGap] = sum(ULP * rLP * 10) / meanTimePerAction
  Rstep[lastWaitGap] = sum(ULP * rLP * 10)
}

plot(Rt / rewardRate$LP)
Rstep[200] * gamma ^ stepDuration / (1 - gamma ^ stepDuration)

# plot state values
tau = 10
nGap = length(trialGapValues$LP)
waitRate_ = vector(mode = "list", length = nGap)
stateValueLP_ = vector(mode = "list", length = nGap)
for(lastWaitGap in 1 : nGap){
  Qwait = Qwait_[[lastWaitGap]]
  Qquit = Qquit_[[lastWaitGap]]
  waitRate = unlist(lapply(1 : lastWaitGap, function(i) exp(tau * Qwait[i]) / sum(exp(tau * Qwait[i]), exp(tau * Qquit))))
  stateValueLP = unlist(lapply(1 : lastWaitGap, function(i) waitRate[i] * Qwait[i] + (1 - waitRate[i]) * Qquit))
  stateValueLP_[[lastWaitGap]] = stateValueLP
}




