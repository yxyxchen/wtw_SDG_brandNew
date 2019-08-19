# libraries and scripts
library("stringr")
library("ggplot2")
library("dplyr") 
library("tidyr")
source("subFxs/helpFxs.R")
source("subFxs/loadFxs.R")
source("subFxs/modelComparisonFxs.R")
source("subFxs/plotThemes.R")
load("wtwSettings.RData")
# load model names
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
load("genData/expDataAnalysis/sessionData.RData")
# select common useID
ids = hdrData$ID[hdrData$stress == "no stress"]
n  = length(ids)
modelNames = factor(c("QL1", "QL2", "RL1", "RL2"),
                    levels = c("QL1", "QL2", "RL1", "RL2", "BL"))
nModel = length(modelNames)
useID_ = vector(mode = "list", length = nModel)
source("subFxs/loadFxs.R")
for(i in 1 : nModel){
  modelName = modelNames[i]
  paraNames = getParaNames(modelName)
  expPara = loadExpPara(paraNames, sprintf("genData/expModelFitting/%sdb", modelName))
  useID_[[i]] = getUseID(expPara, paraNames)
}
useID = ids[apply(sapply(1 : nModel, function(i )ids %in% useID_[[i]]), MARGIN = 1,
              all)]
nUse = length(useID)

# extract logEvidece_ from loo 
logEvidence_ = matrix(NA, nUse, nModel)
logLik_ = matrix(NA, nUse, nModel)
pWaic_ = matrix(NA, nUse, nModel)
for(m in 1 : nModel){
  modelName = modelNames[m]
  for(sIdx in 1 : nUse ){
    id = useID[sIdx]
    fileName = sprintf("genData/expModelFitting/%sdb/s%s_waic.RData", modelName, id)
    load(fileName)
    logEvidence_[sIdx, m] = WAIC$elpd_waic # here is like loglikelyhood, larger the better 
    logLik_[sIdx, m] = WAIC$elpd_waic  + WAIC$p_waic / 2
    pWaic_[sIdx, m] = WAIC$p_waic
  }
}
# save output for modelComparision 
output = data.frame(logEvidence_,
                    condition = ifelse(hdrData$condition[hdrData$ID %in% useID] == "HP", 1, 2),
                    AUC = sessionData$AUC[sessionData$id %in% useID])
f= "genData/expModelFitting/logEvidenceList.csv"
write.table(file = f, output, sep = ",", col.names = F, row.names = F)

# participants best desribed by 
library("ggpubr")
bestNums = sapply(1 : nModel, function(i) sum(apply(logEvidence_[,1:nModel], MARGIN = 1, FUN = function(x) which.max(x) == i)))
data.frame(model = modelNames, bestNums = bestNums) %>%  ggplot(aes(x="", y=bestNums, fill=model)) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) + ylab("") + xlab("") + ggtitle(sprintf("Best described (n = %d)", nUse))+ 
  myTheme
dir.create("figures/expModelComparison")
ggsave("figures/expModelComparison/loo_nBest.png", width = 5, height = 3.5)

data.frame(pwaic = as.vector(pWaic_), model = rep(modelNames, each = nUse)) %>%
  group_by(model) %>% ggplot(aes(model, pwaic)) + geom_boxplot() + myTheme
ggsave("figures/expModelComparison/loo_pwaic.png", width = 5, height =3.5)

data.frame(pwaic = as.vector(pWaic_), model = rep(modelNames, each = nUse)) %>%
  group_by(model) %>% 
  summarise(muData = mean(pwaic), seData = sd(pwaic) / sqrt(length(pwaic)),
            minData = muData - seData, maxData = muData + seData) 

# one and two learning rates
ratioQL2 = sum(logEvidence_[,modelNames == 'QL2'] >= logEvidence_[,modelNames == 'QL1']) / n
ratioRL2 = sum(logEvidence_[,modelNames == 'RL2'] >= logEvidence_[,modelNames == 'RL1']) / n
ratioRL2QL2= sum(logEvidence_[,modelNames == 'RL2'] > logEvidence_[, modelNames == 'QL2']) / n
data.frame(ratio = c(ratioQL2, ratioRL2, ratioRL2QL2),
           model = factor(c('QL', 'RL', 'QLRL'), levels = c('QL', 'RL', 'QLRL')))%>%
  ggplot(aes(model, ratio)) + geom_bar(stat = 'identity') + xlab(' ') + ylab('Ratio') + 
  scale_x_discrete(labels = c('LL_QL2 > LL_QL1', 'LL_RL2 > LL_RL1', 'LL_RL2 > LL_QL2')) + myTheme 
ggsave('figures/expModelComparison/loo_ratio.png', width = 8, height = 3)


# extract logEvidence, cross validation
ids = hdrData$ID[hdrData$stress == "no stress"]; nSub = length(ids)
modelNames = factor(c("QL1", "QL2", "RL1", "RL2", "BL"),
                    levels = c("QL1", "QL2", "RL1", "RL2", "BL"))
nModel = length(modelNames)
nFold = 10
logEvidence = matrix(nrow = length(ids), ncol= nModel) 
logEvidenceTrain = list(length = nModel)
for(mIdx in 1 : nModel){
  modelName = modelNames[mIdx]
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  likFun = getLikFun(modelName)
  for(sIdx in 1 : nSub){
    id = ids[sIdx]
    load(sprintf("genData/expModelFittingCV/split/s%s.RData", id))
    thisTrialData = trialData[[id]]
    cond = unique(thisTrialData$condition)
    cIdx = ifelse(cond == "HP", 1, 2)
    excludedTrials = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[cIdx]))
    # prepare
    nTrial = length(thisTrialData$trialEarnings)
    tMax = ifelse(cond == "HP", tMaxs[1], tMaxs[2])
    trialEarnings = thisTrialData$trialEarnings
    scheduledWait = thisTrialData$scheduledWait
    timeWaited = thisTrialData$timeWaited
    timeWaited[trialEarnings != 0] = scheduledWait[trialEarnings != 0]
    Ts = round(ceiling(timeWaited / stepDuration) + 1)
    cvPara = loadCVPara(paraNames,
                      sprintf("genData/expModelFittingCV/%sdb",modelName),
                      pattern = sprintf("s%s_f[0-9]{1,2}_summary.txt", id))
    # initialize 
    LL_ = vector(length = nFold)
    if(length(getUseID(cvPara, paraNames)) == 10){
      for(f in 1 : nFold){
        # determine training end testing trials
        trials = partTable[f,]
        trials = trials[trials < nTrial]
        junk = 1 : nTrial
        paras = as.double(cvPara[f,1:nPara])
        lik_ = likFun(paras, cond, trialEarnings, timeWaited)$lik_
        LL_[f] = sum(sapply(1 : length(trials), function(i){
          trial = trials[i]
          if(trialEarnings[trial] > 0){
            junk = log(lik_[1 : max(Ts[trial]-1, 1), trial])
            junk[is.infinite(junk)] = -10000
            sum(junk)
          }else{
            if(Ts[trial] > 2){
              junk = c(log(lik_[1:max(Ts[trial] - 2,1), trial]), log(1-lik_[Ts[trial] - 1, trial]))
              junk[is.infinite(junk)] = -10000
              sum(junk)
            }else{
              junk = log(1-lik_[Ts[trial] - 1, trial])
              junk
            }
          }
        }))
        logEvidence[sIdx, mIdx] = sum(LL_)
      }
    }
  }
}
select = apply(sapply(1 : nModel, function(i) !is.na(logEvidence[,i])), MARGIN = 1, FUN = all)
useID = ids[select]

output = data.frame(cvLik = logEvidence[select,],
                    cond = ifelse(sessionData$condition[sessionData$id %in% useID] == "HP", 1, 2),
                    AUC = sessionData$AUC[sessionData$id %in% useID], id = useID)
f= "genData/expModelFitting/logEvidenceListCV.csv"
write.table(file = f, output, sep = ",", col.names = F, row.names = F)

bestNums = sapply(1 : nModel, function(i) sum(apply(logEvidence[,1:nModel], MARGIN = 1, FUN = function(x) which.max(x) == i)))
data.frame(model = modelNames, bestNums = bestNums) %>%  ggplot(aes(x="", y=bestNums, fill=model)) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) + ylab("") + xlab("") +
  ggtitle(sprintf("Participants best described (n = %d)", length(useID)))+ 
  myTheme
dir.create("figures/expModelComparison")
ggsave("figures/expModelComparison/CV_nBest.png", width = 5, height = 3.5)

