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
modelNames = c("PRbs", "PRbsNC", "Rlearn", "RlearnL", "reduce_gamma")
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
modelNames = c("PRbs", "PRbsNC", "Rlearn", "RlearnL", "reduce_gamma")
nModel = length(modelNames)
ids = hdrData$ID[hdrData$stress == "no stress"]
nSub = length(ids)
nFold = 10
logEvidence = matrix(nrow = length(ids), ncol= nModel) 
logEvidenceTrain = list(length = nModel)
for(mIdx in 4 : nModel){
  modelName = modelNames[mIdx]
  paras = getParas(modelName)
  nPara = length(paras)
  logLikFun = getLogLikFun(modelName)
  thisLogEvidenceTrain = matrix(nrow = nFold, ncol = nSub)
  for(sIdx in 1 : nSub){
    id = ids[sIdx]
    load(sprintf("genData/expModelFittingCV/split/s%d.RData", id))
    thisTrialData = trialData[[id]]
    nTrial = length(thisTrialData$trialEarnings)
    cond = unique(thisTrialData$condition)
    tMax = ifelse(cond == "HP", tMaxs[1], tMaxs[2])
    trialEarnings = thisTrialData$trialEarnings
    timeWaited = pmin(thisTrialData$timeWaited, tMax)
    Ts = round(ceiling(timeWaited / stepDuration) + 1)
    scheduledWait = thisTrialData$scheduledWait
    cvPara = loadCVPara(paras,
                      sprintf("genData/expModelFittingCV/%s",modelName),
                      pattern = sprintf("s%d_f[0-9]{1,2}_summary.txt", id))
    # initialize 
    LL_ = vector(length = nFold)
    if(length(getUseID(cvPara, paras)) == 10){
      for(f in 1 : nFold){
        # determine training end testing trials
        trials = partTable[f,]
        trials = trials[trials < nTrial]
        junk = 1 : nTrial
        trialsTrain = junk[!junk %in% trials]
        
        thisParas = as.double(cvPara[f,1:nPara])
        lik_ = logLikFun(thisParas, cond, trialEarnings, timeWaited)$lik_
        LL_[f] = sum(sapply(1 : length(trials), function(i){
          trial = trials[i]
          if(trialEarnings[trial] > 0){
            sum(log(lik_[1 : max(Ts[trial]-1, 1), trial]))
          }else{
            sum(log(lik_[1:max(Ts[trial] - 2,1), trial])) + log(1-lik_[Ts[trial] - 1, trial])
          }
        }))
        thisLogEvidenceTrain[f, sIdx] = sum(sapply(1 : length(trialsTrain), function(i){
          trial = trialsTrain[i]
          if(trialEarnings[trial] > 0){
            sum(log(lik_[1 : max(Ts[trial]-1, 1), trial]))
          }else{
            sum(log(lik_[1:max(Ts[trial] - 2,1), trial])) + log(1-lik_[Ts[trial] - 1, trial])
          }
        }))
          
      }
      logEvidence[sIdx, mIdx] = sum(LL_)
      logEvidenceTrain[[mIdx]] = thisLogEvidenceTrain
    }
  }
}
select = apply(sapply(1 : nModel, function(i) !is.na(logEvidence[,i])), MARGIN = 1, FUN = all)
useID = ids[select]

output = data.frame(cvLik = logEvidence[select,],
                    cond = ifelse(sessionData$condition[hdrData$ID %in% useID] == "HP", 1, 2),
                    AUC = sessionData$AUC[hdrData$ID %in% useID], id = useID)
f= "genData/expModelFitting/logEvidenceListCV.csv"
write.table(file = f, output, sep = ",", col.names = F, row.names = F)
nZeros = sapply(1 : nFold, function(i) sum(trialEarnings[(1 : nTrial) %in% as.vector(partTable[-i,])] == 0))
