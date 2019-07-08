# libraries and scripts
library("stringr")
library("ggplot2")
source("subFxs/helpFxs.R")
source("subFxs/loadFxs.R")
source("subFxs/modelComparisonFxs.R")
source("subFxs/plotThemes.R")
# load model names
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
allIDs = hdrData$ID                   # column of subject IDs
n = length(allIDs) 

# select common useID
idList = hdrData$ID
modelNames = c("PRbs", "PRbsNC", "Rlearn", "RlearnL")
nModel = length(modelNames)
useID_ = vector(mode = "list", length = nModel)
source("subFxs/loadFxs.R")
for(i in 1 : nModel){
  modelName = modelNames[i]
  paras = getParas(modelName)
  expPara = loadExpPara(paras, sprintf("genData/expModelFitting/%s", modelName))
  useID_[[i]] = getUseID(expPara, paras)
}
useID = idList[apply(sapply(1 : nModel, function(i )idList %in% useID_[[i]]), MARGIN = 1,
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
    fileName = sprintf("genData/expModelFitting/%s/s%d_waic.RData", modelName, id)
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
  coord_polar("y", start=0) + ylab("") + xlab("") + ggtitle(sprintf("Participants best described (n = %d)", nUse))+ 
  myTheme


# extract logEvidence, cross validation
ids = hdrData$ID[hdrData$stress == "no stress"]
modelName = "PRbs"
paras = getParas(modelName)
nPara = length(paras)
nFold = 10
logLikFun = getLogLikFun(modelName)
logEvidence = vector(length = length(ids))
nCore = parallel::detectCores() -1 # only for the local computer
registerDoMC(nCore)
for(sIdx in 1 : length(ids)){
  id = ids[sIdx]
  load(sprintf("genData/expModelFittingCV/%s/s%d_split.RData",modelName, id))
  thisTrialData = trialData[[id]]
  nTrial = length(thisTrialData$trialEarnings)
  cond = unique(thisTrialData$condition)
  tMax = ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  trialEarnings = thisTrialData$trialEarnings
  timeWaited = pmin(thisTrialData$timeWaited, tMax)
  Ts = round(ceiling(timeWaited / stepDuration) + 1)
  scheduledWait = thisTrialData$scheduledWait
  cvPara = loadCVPara(paras, sprintf("genData/expModelFittingCV/%s",modelName),
                      id)
  # initialize 
  LL_ = vector(length = nFold)
  if(length(getUseID(cvPara, paras)) == 10){
    for(f in 1 : nFold){
      trials = partTable[f,]
      trials = trials[trials < nTrial]
      thisParas = tempt[1:nPara,1]
      logLik_ = logLikFun(thisParas, cond, trialEarnings, timeWaited)$logLik_
      LL_[f] = sum(sapply(1 : length(trials), function(i){
        trial = trials[i]
        if(trialEarnings[trial] > 0){
          sum(logLik_[1 : max(Ts[trial]-1, 1), trial])
        }else{
          sum(logLik_[1:max(Ts[trial] - 2,1), trial]) + log(1-exp(logLik_[Ts[trial] - 1, trial]))
        }
      }))
    }
    logEvidence[sIdx] = sum(LL_)
  }
}


