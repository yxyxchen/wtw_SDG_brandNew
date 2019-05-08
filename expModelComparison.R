# this script extract waic for all models. 

# libraries and scripts
library("stringr")
library("ggplot2")
source("subFxs/helpFxs.R")

# load model names
modelNames = c("curiosityTrial", "curiosityTrialNC")
nModel = length(modelNames)

# load experimental data
load("genData/expDataAnalysis/blockData.RData")
idList = unique(blockData$id) 
n = length(idList)

# define a function for convinence
# select useID
useID_ = vector(mode = "list", length = nModel)
useID = idList
source("subFxs/loadFxs.R")
for(i in 1 : nModel){
  modelName = modelNames[i]
  expPara = loadExpPara(modelName, getParas(modelName))
  pars = getParas(modelName)
  RhatCols = which(str_detect(colnames(expPara), "hat"))[1 : length(pars)]
  EffeCols = which(str_detect(colnames(expPara), "Effe"))[1 : length(pars)]
  if(length(pars) <= 1){
    useID_[[i]]  = idList[expPara[,RhatCols]<1.1 & expPara[,EffeCols] > 100]
  }else{
    useID_[[i]] = idList[apply(expPara[,RhatCols] < 1.1, MARGIN = 1, sum) == length(pars) &
                           apply(expPara[,EffeCols] >100, MARGIN = 1, sum) == length(pars)]
    
  }
  useID = idList[idList %in% useID_[[i]] & idList %in% useID]
}
nUse = length(useID)

# use logLikelyhood of the median paramater
medianLogLik_ = matrix(NA, n, nModel)
for(m in 1 : nModel){
  modelName = modelNames[m]
  pars = getPars(modelName)
  tempt = loadExpParaExtra(modelName, pars)
  expParaMedian = tempt$expParaMedian
  medianLogLik_[,m] = expParaMedian[,(length(pars) +1)]
}
medianLogLik_ = medianLogLik_[expPara$id %in% useID, ]
f= "genData/expModelFitting/logEvidenceList.csv"
write.table(file = f, medianLogLik_, sep = ",", col.names = F, row.names = F)
# extract logEvidence per trial
# here logEvidence is on the loglikelyhood scale, so it is negtive 
# also, it accounts for model complexity 
# logLik_ doesn't account for model complexity
# since p_waic is very noise, therefore pay more attention to logLik_
logEvidence_ = matrix(NA, nUse, nModel)
logLik_ = matrix(NA, nUse, nModel)
logLikPerAct_ = matrix(NA, nUse, nModel)
logEvidencePerAct_ = matrix(NA, nUse, nModel) # save later for model comparison 
pWaic_ = matrix(NA, nUse, nModel)
for(m in 1 : nModel){
  modelName = modelNames[m]
  for(sIdx in 1 : nUse ){
    id = useID[sIdx]
    nAction = blockData$nAction[blockData$id == id & blockData$blockNum == 1]
    fileName = sprintf("genData/expModelFitting/%s/s%d_waic.RData", modelName, id)
    load(fileName)
    logEvidence_[sIdx, m] = WAIC$elpd_waic
    logEvidencePerAct_[sIdx, m] = WAIC$elpd_waic/ nAction
    pWaic_[sIdx, m] = WAIC$p_waic
    paraSummary = read.csv(sprintf("genData/expModelFitting/%s/s%d_summary.txt", modelName, id), header = F)
    logLik_[sIdx, m] = paraSummary[ nrow(paraSummary) - 1, 1]
  }
}
waic_ = -2 * logEvidence_
mean(waic_[,1] - waic_[,2])
f= "genData/expModelFitting/logEvidenceList.csv"
write.table(file = f, logEvidence_, sep = ",", col.names = F, row.names = F)

# is there any logLik differences within a subset of data?
zoomInID = unique(blockData$id[blockData$AUC> 11 & blockData$AUC < 25 & blockData$condition == "LP"])
a = medianLogLik_[useID %in% zoomInID, ]
arbMinusThe = a[,1] - a[,2]
hist(arbMinusThe)


