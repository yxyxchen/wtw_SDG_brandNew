# library 
library('ggplot2')
library('dplyr')
library('tidyr')
load("wtwSettings.RData")
source('subFxs/repetitionFxs.R') # 
source("subFxs/helpFxs.R")
source("subFxs/plotThemes.R")
source("subFxs/loadFxs.R") 
source("subFxs/analysisFxs.R")
source("subFxs/taskFxs.R")
###### simulate using the empirical schedualedWait
# load raw data 
allData = loadAllData()
hdrData = allData$hdrData  
allIDs = hdrData$ID     
expTrialData = allData$trialData   
idList = hdrData$ID
n = length(idList)

##### get a sense of the model ######
# simluation for one para and one scheduledWait from simulation or ...
paras = c(0.05, 30)
modelName = "functionLinear"
repModelFun = getRepModelFun(modelName)
sIdx = 6
id = idList[[sIdx]]
cond = hdrData$cond[hdrData$ID == id]
thisExpTrialData = expTrialData[[id]]
thisExpTrialData = thisExpTrialData[thisExpTrialData$blockNum ==1, ]
scheduledWait = thisExpTrialData$scheduledWait
set.seed(123)
# cond = "LP"
# scheduledWait = unlist(lapply(1:1000, function(x) drawSample(cond)))
tempt = repModelFun(paras, cond, scheduledWait)
trialPlots(tempt, cond)

## try to make some predictions
scheduledWait = thisExpTrialData$scheduledWait
trialEarnings = thisExpTrialData$trialEarnings
timeWaited = thisExpTrialData$timeWaited
nTrial = length(scheduledWait)
timeWaited[trialEarnings > 0] = scheduledWait[trialEarnings > 0]
Ts = round(ceiling(timeWaited / stepDuration) + 1)

yHat = matrix(0, nTimeStep, 42)
beta0Hat = vector(length = 42)
beta1Hat = vector(length = 42)
quitTime = vector(length = 42)

beta0Hat[,1] = 1.2
beta1Hat[,1] = 0
yHat[,1] = 1.2
quitTime[,1] = nTimeStep
for(k in 1 : 41){
  nTrial = k
  x = unlist(lapply(1 : nTrial, function(i) 1 : (Ts[i] - 1)))
  y = unlist(lapply(1 : nTrial,
                    function(i) rev(  (trialEarnings[i] + 1) * 2 / (1 : (Ts[i] - 1)) )))
  fit = lm(y ~ x)
  yHat[, k+1] = predict(fit, newdata=data.frame(x=1:nTimeStep))
  beta0Hat[k+1] = fit$coefficients[1]
  beta1Hat[k+1] = fit$coefficients[2]
  if(sum(tempt$rrBars[k] >= yHat[, k]) == nTimeStep || sum(tempt$rrBars[k] < yHat[, k]) == nTimeStep){
    quitTime[k+1] = NA
  }else{
    quitTime[k+1] = which.min(abs(tempt$rrBars[k] - yHat[, k]))
  }
}


nTimeStep = tMaxs[2] / stepDuration
beta.prior = matrix(c(1.2, 0), ncol = 1)
sigmaSq.prior = matrix(c(0.1, 0, 0, 0.1), nrow = 2, ncol = 2)
x.star =  cbind(rep(1, nTimeStep), (1 : nTimeStep - 1) * stepDuration)s
yHat = matrix(0, nTimeStep, 42)
quitTime = vector(length = 42)
yHat[, 1] =  x.star %*% beta.prior
quitTime[1] = nTimeStep
for(k in 1 : 41){
  nTrial = k
  x = unlist(lapply(1 : nTrial, function(i) 1 : (Ts[i] - 1)))
  y = unlist(lapply(1 : nTrial,
                    function(i) rev(  (trialEarnings[i] + 1) * 2 / (1 : (Ts[i] - 1)) )))
  X = cbind(rep(1, length(x)), x)
  Y = matrix(y, ncol = 1)
  sigmaSq = sum((y - X %*% solve(t(X) %*% X) %*% t(X) %*% Y)^2) / (length(y) - 2)
  A = solve(sigmaSq.prior) + t(X) %*% X / sigmaSq 
  betaHat.mu = solve(A)  %*% (t(X) %*% Y/sigmaSq + solve(sigmaSq.prior) %*% beta.prior)
  yHat.mu = x.star %*% betaHat.mu
  yHat[, k+1] = yHat.mu
  
  if(sum(tempt$rrBars[k] >= yHat[, k]) == nTimeStep || sum(tempt$rrBars[k] < yHat[, k]) == nTimeStep){
    quitTime[k+1] = NA
  }else{
    quitTime[k+1] = which.min(abs(tempt$rrBars[k] - yHat[, k]))
  }
}
plot(yHat[,20], ylim = c(-5, 5))

# simulate 
source('subFxs/simulationFxs.R') # 
paraTable = list("phi" = c(0.025, 0.05, 0.1), "tau" = c(5, 25, 45))
cond = "LP"
nTimeStep = tMaxs[2] / stepDuration
nComb = 3^2
nRep = 10
nSeq = 20

# simulates
seqLens = c(5, 25, 100)
rrs.mu = vector(length = nTimeStep * 3)
rrs.std = vector(length = nTimeStep * 3)
set.seed(123)
for(lIdx in 1 : 3){
  seqLen = seqLens[lIdx]
  scheduledWaitList = replicate(nSeq, replicate(seqLen, drawSample(cond)), simplify = F)
  trialData = simulate("functionRL", nRep, paraTable, scheduledWaitList, cond)
  # extract rrs
  rrs_ = array(t(seq(1 : (nComb * nSeq * nRep * nTimeStep))), dim = c(nRep, nSeq, nComb, nTimeStep)) 
  load("genData/simulation/functionRL/simParas.RData")
  for(cIdx in 1 : nComb){
    for(sIdx in 1 : nSeq){
      for(rIdx in 1 : nRep){
        thisTrialData = trialData[[simNo[rIdx, sIdx, cIdx]]]
        rrs_[rIdx, sIdx, cIdx, ] = thisTrialData$rrs[,seqLen]
      }
    }
  }
  rrs.mu[ (nTimeStep * (lIdx - 1) + 1) : (nTimeStep * lIdx)] = apply(rrs_, 4, mean)
  rrs.std[ (nTimeStep * (lIdx - 1) + 1) : (nTimeStep * lIdx)] = apply(rrs_, 4, sd)
}

plotData = data.frame(rrs.mu, rrs.std, rrs.max = rrs.mu + rrs.std, rrs.min = rrs.mu - rrs.std,
                      trial.number = rep(seqLens, each = nTimeStep), 
                      time = rep((1 : nTimeStep) * stepDuration, 3))

# not smooth enough... shall I show quitting time underlying different tau and different time threshold
ggplot(plotData, aes(time, rrs.mu)) + facet_grid(~trial.number) +
  geom_errorbar(aes(ymin = rrs.min, ymax = rrs.max), width = 0.2, color = "grey") + saveTheme + 
  geom_line(size= 1) + xlab("Time / s") + ylab("Reward rate")
ggsave("rrs_RL.png", width = 8, height = 4)


##### get a sense of the model from a lot scheduledWait######
# simluation for one para and a lot scheduledWait
paras = c(0.02, 10, 0.001) 
modelName = "curiosityTrialR"
repModelFun = getRepModelFun(modelName)
nRep = 1 # number of repetitions
trialData = vector(length = n * nRep, mode ='list')
repNo = matrix(1 : (n * nRep), nrow = n, ncol = nRep)
set.seed(231)
for(sIdx in 1 : n){
  id = idList[[sIdx]] 
  cond = hdrData$cond[hdrData$ID == id]
  thisExpTrialData = expTrialData[[id]]
  thisExpTrialData = thisExpTrialData[thisExpTrialData$blockNum ==1, ]
  scheduledWait = thisExpTrialData$scheduledWait
  for(rIdx in 1 : nRep){
    tempt = repModelFun(paras, cond, scheduledWait)
    trialData[[repNo[sIdx, rIdx]]] = tempt
  }
}

for(i in 1:n){
  sIdx = repNo[i, 1]
  id = idList[sIdx]
  label = sprintf("%s,%s", hdrData$condition[hdrData$ID== i], id)
  thisTrialData = trialData[[sIdx]]
  trialPlots(thisTrialData, label)
  readline("continue")
  tMax = ifelse(label == "HP", 20, 40)
  kmGrid = seq(0, tMax, by=0.1) 
  thisTrialData = trialData[[sIdx]]
  kmscResults = kmsc(thisTrialData,tMax,label,T,kmGrid)
  readline("continue")
  # actionValueViewer(thisTrialData)
}


# one sequences example
paras = c(0.01, 15, 1)
set.seed(123)
modelName = "curiosityTrial"
repModelFun = getRepModelFun(modelName)
longScheduledWaitLP =  unlist(lapply(1:5000, function(x) drawSample("LP")))
tempt = repModelFun(paras, "LP", longScheduledWaitLP)
label = sprintf("tau = %.2f", paras[2])
trialPlots(tempt, label)
fileName = sprintf("zoomIn_%s.png", label)
trialPlots(truncateTrials(tempt, 4800, 5000), label)
ggsave(fileName, width = 6, height = 12)
actionValueViewer(truncateTrials(tempt, 4000, 4100))
plot(tempt$Qwaits[,4000])
which.min(tempt$Qwaits[,2000])

# simulated scheduledWait
set.seed(123)
for(c in 1 : 2){
  cond  = conditions[c]
  if(cond == "HP"){
    scheduledWaitHPlist = lapply(1 : nRep, function(i) unlist(lapply(1:50, function(x) drawSample(cond))))
  }else{
    scheduledWaitLPlist = lapply(1 : nRep, function(i) unlist(lapply(1:50, function(x) drawSample(cond))))
  }
}

# simulate
modelName = "full_model" 
nBlock = 1
nRep = 10
paraTable = data.frame(phi = c(0.02, 0.05, 0.08), tau = c(5, 10, 15),
                       gamma = c(0.85, 0.90, 0.95), QwaitIni = c(2, 3, 4))
simulate(modelName, nBlock, nRep, paraTable)


modelName = "R_learning2" 
nBlock = 1
nRep = 10
paraTable = data.frame(phi1 = c(0.02, 0.05, 0.08), phi2 = c(0.02, 0.05, 0.08),tau = c(5, 10, 15))
simulate(modelName, nBlock, nRep, paraTable)



