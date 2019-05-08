# this script is used to unify all algorithms. 
# here we assume the rewards arrive at the end
# if we need stats for decrete steps, we can use wtwSettings.
# for the optimal statistics, don't use both of them!
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
k = 4 # shape
mu = 0 # location
sigma = 2 # scale
pareto = list()
pareto[['k']] = k # shape
pareto[['mu']] = mu # location
pareto[['sigma']] = sigma # scale
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
rewardRates = list('HP' = HP, 'LP' = LP)

optimWaitTimes = list()
optimWaitTimes$HP = trialGapValues$HP[which.max(HP)]
optimWaitTimes$LP = trialGapValues$LP[which.max(LP)]

optimRewardRates = list()
optimRewardRates$HP = max(HP)
optimRewardRates$LP = max(LP)

# the mean waited time without ITI, counting from the every begining of the gap.
getMeanRewardDelay = function(gapIdx, quitAfter, gammaPerStep, cond){
  gamma = gammaPerStep^2 # discount per 1s
  nGap = length(trialTicks[[cond]]) - 1
  nTick =  length(trialTicks[[cond]]) 
  quitGap = round(quitAfter / stepDuration)
  if(quitGap < gapIdx){
    discount = NaN
    discount = NaN
  }else{
    discount =sum(gamma ^ (seq(0, quitGap - gapIdx) * stepDuration) * rewardDelayPDF[[cond]][gapIdx : quitGap]) /
      sum(rewardDelayPDF[[cond]][gapIdx : quitGap])
    time = sum((seq(1, quitGap - gapIdx + 1) * stepDuration) * rewardDelayPDF[[cond]][gapIdx : quitGap]) /
      sum(rewardDelayPDF[[cond]][gapIdx : quitGap])
  }
  # meanRewardDelay = list(discount = discount, time = time)
  meanRewardDelay = list(discount = discount * gamma ^ stepDuration, time = time)
  return(meanRewardDelay)
}

# what it is this?
getMeanWaitDelay = function(gapIdx, quitAfter, gammaPerStep, cond){
  gamma = gammaPerStep^2 # discount per 1s
  nGap = length(trialTicks[[cond]]) - 1
  nTick =  length(trialTicks[[cond]]) 
  tickIdx = gapIdx + 1
  quitGap = round(quitAfter / stepDuration)
  quitTick = quitGap  + 1
  if(quitTick < tickIdx){
    discount = NaN
    time = NaN
  }else{
    discount =sum(gamma ^ pmin(seq(0, nGap - gapIdx) * stepDuration, quitAfter - gapIdx * stepDuration)
                  * rewardDelayPDF[[cond]][gapIdx : nGap]) /
      sum(rewardDelayPDF[[cond]][gapIdx : nGap]) 
    time = sum(pmin(seq(1, 1 + nGap - gapIdx) * stepDuration, quitAfter - gapIdx * stepDuration ) * rewardDelayPDF[[cond]][gapIdx : nGap]) /
      sum(rewardDelayPDF[[cond]][gapIdx : nGap]) + iti
  }
  # meanWaitDelay = list(discount = discount, time = time)
  meanWaitDelay = list(discount = discount * gamma ^ stepDuration, time = time)
  return(meanWaitDelay)
}

# the prob of getting rewards at the very beggining of the gap. If gapIdx = 1, quitAfter = 10, p = 0.5
getPreward = function(gapIdx, quitAfter, gammaPerStep, cond){
  gamma = gammaPerStep^2 # discount per 1s
  nGap = length(trialTicks[[cond]]) - 1
  nTick =  length(trialTicks[[cond]]) 
  tickIdx = gapIdx + 1
  quitGap = round(quitAfter / stepDuration)
  if(gapIdx == 1){
    pReward = rewardDelayCDF[[cond]][quitGap]  
  }else{
    pReward = (rewardDelayCDF[[cond]][quitGap] - rewardDelayCDF[[cond]][gapIdx - 1]) / (1 - rewardDelayCDF[[cond]][gapIdx - 1])
  }
  return(pReward)
}

cond = "LP"
gammaPerStep = 0.80 # discount per 0.5s 
gamma = gammaPerStep ^ 2
kd = -log(gamma)
nGap = length(trialGapValues[[cond]]) 
nTick =  length(trialTicks[[cond]]) 
Qquit_ = vector(length = nGap)
Qwait_ = vector(mode = "list", length = nGap)
# calculate rewardRate rr_
rr_ = vector(mode = "list", length = nGap)
potential_ = vector(mode = "list", length = nGap)
rewardTimes_ = vector(mode = "list", length = nGap)
pRewards_ = vector(mode = "list", length = nGap)
rrBar_ = vector(length = nGap)
for(quitGap in 1 : nGap){
  quitAfter = stepDuration * quitGap
  meanWaitDelay = getMeanWaitDelay(1, quitAfter, gammaPerStep, cond)
  pReward = getPreward(1, quitAfter, gammaPerStep, cond)
  rrBar = pReward * tokenValue / (meanWaitDelay$time- 0.5 * stepDuration)
  rrBar_[[quitGap]] = rrBar
  
  rr = vector(length = quitGap)
  potential = vector(length = quitGap)
  pRewards = vector(length = quitGap)
  rewardTimes  = vector(length = quitGap)
  for(j in 1 : quitGap){
    meanRewardDelay = getMeanRewardDelay(j, quitAfter, gammaPerStep, cond)
    rewardTime = meanRewardDelay$time
    pReward = getPreward(j, quitAfter, gammaPerStep, cond)
    
    pRewards[j] = pReward
    rewardTimes[j] = meanRewardDelay$time
    rr[j] = tokenValue * pReward  / rewardTime
    potential[j] = tokenValue * pReward - rrBar * (rewardTime - 0.5 * stepDuration)
    # potential[j] = tokenValue * pReward - rrBar * rewardTime 
  }
  rr_[[quitGap]] = rr
  potential_[[quitGap]] = potential
  pRewards_[[quitGap]] = pRewards
  rewardTimes_[[quitGap]] = rewardTimes
}
# rr_HP = rr_
# rr_LP = rr_
exits = sapply(1 : nGap, function(i) which(potential_[[i]] < 0)[1])
plotData = data.frame(policy = trialGapValues$LP, quitTime = exits * stepDuration)
ggplot(plotData, aes(policy, quitTime)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) + xlim(c(0, 5)) + ylim(c(0, 5))

potentialMatrix = matrix(NA, nrow = 400, ncol = 400)
for(i in 1 : 400){
  potentialMatrix[1:i, i] = potential_[[i]]
}
policyList = c(140, 220, 300, 380)
nPolicy = length(policyList)
plotData = data.frame(potential = unlist(lapply(1:nPolicy,
                                         function(i) potentialMatrix[,policyList[i]])),
                      exist_threshold = as.factor(rep(policyList, each = 400)*stepDuration),
                      time = rep(1:400, nPolicy))
ggplot(plotData, aes(time, potential, color = exist_threshold)) + geom_point(size = 1) +
  ylab("Potential") + xlab("Time / s") + saveTheme + geom_hline(yintercept = 0)
ggsave("potenTialLP.png",width = 8, height = 4)

