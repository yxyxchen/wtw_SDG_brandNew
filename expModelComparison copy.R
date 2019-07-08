# libraries and scripts
library("stringr")
library("ggplot2")
source("subFxs/helpFxs.R")
source("subFxs/loadFxs.R")
# load model names
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
allIDs = hdrData$ID                   # column of subject IDs
n = length(allIDs) 


likRatioTest("PR", "PRNC", 1, "all")
likRatioTest("fullModel", "PR", 1, "all")
likRatioTest("baseline", "PR", 4, "all")

likRatioTest = function(modelName1, modelName2, df, group = "all"){
  paras1 = getParas(modelName1)
  expPara1 = loadExpPara(paras1,
                         sprintf("genData/expModelFitting/%s", modelName1))
  useID1 = getUseID(expPara1, paras1)
  paras2 = getParas(modelName2)
  expPara2 = loadExpPara(paras2,
                         sprintf("genData/expModelFitting/%s", modelName2))
  useID2 = getUseID(expPara2, paras2)
  if(group == "HP" || group == "LP"){
    useID = allIDs[allIDs %in% useID1 & allIDs %in% useID2 & hdrData$condition == group]
  }
  
  sumLogEvi1 = filter(expPara1, id %in% useID) %>% summarise(sum(LL_all))
  sumLogEvi2 = filter(expPara2, id %in% useID) %>% summarise(sum(LL_all))
  
  delta = -2 * as.double(sumLogEvi1 - sumLogEvi2)
  p = pchisq(delta, df * length(useID), lower.tail=FALSE)
  return(p)
}

# define a function for convinence
# select useID
idList = hdrData$ID
modelNames = c("PRbs", "PRbsNC")
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

logEvidence_ = matrix(NA, nUse, nModel)
logLik_ = matrix(NA, nUse, nModel)
pWaic_ = matrix(NA, nUse, nModel)
for(m in 1 : nModel){
  modelName = modelNames[m]
  for(sIdx in 1 : nUse ){
    id = useID[sIdx]
    fileName = sprintf("genData/expModelFitting/%s/s%d_waic.RData", modelName, id)
    load(fileName)
    logEvidence_[sIdx, m] = WAIC$elpd_waic
    logLik_[sIdx, m] = WAIC$elpd_waic  + WAIC$p_waic / 2
    pWaic_[sIdx, m] = WAIC$p_waic
  }
}
output = data.frame(logEvidence_,
                    condition = ifelse(hdrData$condition[hdrData$ID %in% useID] == "HP", 1, 2))
f= "genData/expModelFittingSub/logEvidenceList.csv"
write.table(file = f, output, sep = ",", col.names = F, row.names = F)

library("ggpubr")
condition = hdrData$condition[hdrData$ID %in% useID]
nQuit = sessionData$nQuit[hdrData$ID %in% useID]
data.frame(logLik = as.vector(logLik_), condition = rep(condition, nModel),
           model = as.factor(rep(modelNames, each = nrow(logLik_))),
           nQuit =  rep(nQuit, nModel)) %>% 
  filter(nQuit < 5) %>% ggplot(aes(model, logLik)) +
  geom_boxplot() + stat_compare_means()



# is there any logLik differences within a subset of data?
zoomInID = unique(blockData$id[blockData$AUC> 11 & blockData$AUC < 25 & blockData$condition == "LP"])
a = medianLogLik_[useID %in% zoomInID, ]
arbMinusThe = a[,1] - a[,2]
hist(arbMinusThe)


