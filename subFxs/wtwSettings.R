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
loseValue = 0
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
optimWaitTimes$HP = trialGapValues$HP[which.max(HP)]
optimWaitTimes$LP = trialGapValues$LP[which.max(LP)]

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
wInisArbitrary = vector(length = 2)
wInisArbitrary[1] = 3 
wInisArbitrary[2] = 2 

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