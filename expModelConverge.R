source("subFxs/loadFxs.R")

# parameters
nFold = 10

# load sessionData             
load("genData/expDataAnalysis/sessionData.RData")
nSub = sum(sessionData$stress == "no stress")
ids = 1 : (nSub * nFold)

# initialize outputs
nFits = matrix(rep(1, nrow * ncol), nrow = nFold, ncol = nSub)

# loop over models
modelNames = c("PRbs", "PRbsNC", "Rlearn", "RlearnL", "reduce_gamma")

mIdx = 1
modelName = modelNames[mIdx]
paras = getParas(modelName)
cvPara = loadCVPara(paras,
                  sprintf("genData/expModelFittingCV/%s", modelName),
                  "*.txt")
useID = getUseID(cvPara, paras)
excID = ids[!ids %in% useID]
for(i in 1 : length(excID))


