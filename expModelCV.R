
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


