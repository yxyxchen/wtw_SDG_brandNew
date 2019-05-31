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

# likRatioTest("para4", "PR", 1, "HP")
# likRatioTest("PR", "PR_cost", 1, "HP")
# likRatioTest("para4", "PR", 1, "LP")
# likRatioTest("PR", "PR_cost", 1, "LP")
likRatioTest("para4", "PR", 1, "all")
likRatioTest("fullModel", "PR", 1, "all")
likRatioTest("baseline", "PR", 4, "all")

likRatioTest = function(modelName1, modelName2, df, group = "all"){
  paras1 = getParas(modelName1)
  expPara1 = loadExpPara(paras1,
                         sprintf("genData/expModelFittingSub/%s", modelName1))
  useID1 = getUseID(expPara1, paras1)
  paras2 = getParas(modelName2)
  expPara2 = loadExpPara(paras2,
                         sprintf("genData/expModelFittingSub/%s", modelName2))
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
modelNames = c("PR", "baseline", "para4", "fullModel")
nModel = length(modelNames)
useID_ = vector(mode = "list", length = nModel)
useID = idList
source("subFxs/loadFxs.R")
for(i in 1 : nModel){
  modelName = modelNames[i]
  paras = getParas(modelName)
  expPara = loadExpPara(paras, sprintf("genData/expModelFittingSub/%s", modelName))
  useID_[[i]] = getUseID(expPara, paras)
}
useID = idList[idList %in% useID_[[i]] & idList %in% useID]
nUse = length(useID)

logEvidence_ = matrix(NA, nUse, nModel)
pWaic_ = matrix(NA, nUse, nModel)
for(m in 1 : nModel){
  modelName = modelNames[m]
  for(sIdx in 1 : nUse ){
    id = useID[sIdx]
    fileName = sprintf("genData/expModelFittingSub/%s/s%d_waic.RData", modelName, id)
    load(fileName)
    logEvidence_[sIdx, m] = WAIC$elpd_waic
    pWaic_[sIdx, m] = WAIC$p_waic
  }
}
output = data.frame(logEvidence_,
                    condition = ifelse(hdrData$condition[hdrData$ID %in% useID] == "HP", 1, 2))
f= "genData/expModelFittingSub/logEvidenceList.csv"
write.table(file = f, output, sep = ",", col.names = F, row.names = F)

meanWAIC = (-2 * logEvidence_) %>% apply(MARGIN  = 2, FUN = mean)
meanWAIC[2:4] - meanWAIC[1] 

png('diffBICHP.png', width = 400, height = 250,
    units = "px")
hist((waic_[output$condition == "HP",1] - waic_[output$condition == "HP",2]),
     xlab="BIC_Qlearn - BIC_Rlearn",
     xlim=c(-450, 40), main = "HP", breaks = 3)
dev.off()
png('diffBICLP.png', width = 400, height = 250,
    units = "px")
hist((waic_[output$condition == "LP",1] - waic_[output$condition == "LP",2]),
     xlab="BIC_Qlearn - BIC_Rlearn",
     xlim=c(-450, 40), main = "LP", breaks = 30)
dev.off()
# is there any logLik differences within a subset of data?
zoomInID = unique(blockData$id[blockData$AUC> 11 & blockData$AUC < 25 & blockData$condition == "LP"])
a = medianLogLik_[useID %in% zoomInID, ]
arbMinusThe = a[,1] - a[,2]
hist(arbMinusThe)


