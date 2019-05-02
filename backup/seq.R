paraTable = list("phi" = 0.02, "tau" = 20)
set.seed(123)
scheduledWait = unlist(lapply(1:100, function(i) drawSample(cond)))
scheduledWaitList = lapply(1 : 10, function(i) scheduledWait[sample(100)])
trialData = simulate(modelName, nRep, paraTable, scheduledWaitList, "LP")
nRep = 10
modelName = "functionRL"

load(sprintf("genData/simulation/%s/simParas.RData", modelName))
plotKMSC = F
simNo = simNo[,,1]
# initialize 
totalEarningsRep_ = matrix(0, n, nRep)
AUCRep_ = matrix(0, 10, nRep)
timeWaitedRep_ = vector(mode = "list", length = n)
for(sIdx in 1 : 10){
  tMax = tMaxs[2]
  kmGrid = seq(0, tMax, by=0.1) 
  label = "sa"
  # initialize
  for(rIdx in 1 : nRep){
    thisTrialData = trialData[[simNo[rIdx, sIdx]]]
    junk = thisTrialData$timeWaited
    kmscResults = kmsc(thisTrialData,tMax,label,plotKMSC,kmGrid)
    AUCRep_[rIdx, sIdx] = kmscResults[['auc']]
  }
}

AUCRep_ = as.data.frame(AUCRep_)
names(AUCRep_) = paste0(1:10)

plotData = gather(AUCRep_, "seq", "AUC")
ggplot(plotData, aes(seq, AUC)  ) + geom_boxplot() +
  xlab("Seq No") + ylab("AUC / secs") + saveTheme
ggsave("seq_effect.png", width = 6, height = 4)
