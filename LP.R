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

gammaPerStep = 0.5 # discount per 0.5s 
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

# find the cross point for every policy(where Qquit > Qquit)
crossPoints = vector(length = nGap)
for(lastWaitGap in 1 : nGap){
  tempt = which(Qwait_[[lastWaitGap]] <= Qquit_[[lastWaitGap]])
  if(length(tempt) == 0) crossPoints[lastWaitGap] = min(nGap, lastWaitGap)* stepDuration 
  else crossPoints[lastWaitGap] = min(tempt * stepDuration) 
}
plot(crossPoints)
max(crossPoints)
which.max(crossPoints)

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

# the same
# always true, even the rewards occurs at the middle of the gap
plot(Rstep / 0.1, rewardRate$LP)

# plot state values
# only true when all action have the same time
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
plot(Rstep * gamma ^ stepDuration / (1 - gamma ^ stepDuration),
     QLP)


plot(rewardRate$LP*stepDuration * gamma^stepDuration / (1 - gamma^stepDuration),
     QLP)

stepDuration = 0.001
QLPAp = rewardRate$LP * stepDuration * gamma ^ stepDuration / (1 - gamma ^ stepDuration)
QLPAp[[22]]