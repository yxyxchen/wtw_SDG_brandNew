# model-free analysis

# inputs:
# isTrct : logical variable determining whether the last portion in each block is truncated 

# outputs:
# sumStats = {
  # id : [nSubx1 id]
  # condition : [nSubx1 fac]
  # nExcl : [nSubx1 int] # total number of excluded trials 
  # muWTW : [nSubx1 num] # average willingness to wait (WTW), measured by area under the Kaplan-Meier survival curve
  # stdWTW : [nSubx1 num] # standard deviation of WTW, measured in Kaplan-Meier survival analysis
  # totalEarnings :  [nSubx1 num] 
# }
# timeWTW_ : list(nSubx1) # wtw timecourse, each element is a vector
# trialWTW_ : list(nSubx1) # trial-wise WTW, each element is a vector
# survCurve_ : list(nSubx1) # Kaplan-Meier survival curve, each element is a vector

MFAnalysis = function(isTrct){
  # load libraries
  source('subFxs/loadFxs.R') 
  source('subFxs/analysisFxs.R') 
  library('dplyr')
  library("tidyr")
  
  # create the output directory
  dir.create("genData")
  dir.create("genData/MFAnalysis")
  
  # load experiment parameters
  load("expParas.RData")
  
  # load exp data
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  ids = hdrData$id[hdrData$stress == "no_stress"]  
  nSub = length(ids)                    # n
  cat('Analyzing data for',nSub,'subjects.\n')
  
  # initialize output variables 
  nExcl = numeric(length = nSub)
  muWTW = numeric(length = nSub) 
  stdWTW = numeric(length = nSub) 
  totalEarnings =  numeric(length = nSub) 
  timeWTW_ = vector(mode = "list", length = nSub) 
  trialWTW_ = vector(mode = "list", length = nSub) 
  survCurve_ = vector(mode = "list", length = nSub) 
  
  # loop over inidviduals
  for (sIdx in 1 : nSub) {
    # load (and truncate) trialData for this individual
    id = ids[sIdx]
    thisTrialData = trialData[[id]]
    if(isTrct){
      trctLine = blockSec - max(tMaxs)
      # truncate trials completed after tractline in each block
      nExcl[sIdx] = sum(thisTrialData$trialStartTime > trctLine)
      thisTrialData = thisTrialData %>% filter(trialStartTime <=  trctLine )
    }else{
      nExcl[sIdx] = 0
    }
    # calcualte totalEarnings
    totalEarnings[sIdx] =  sum(thisTrialData$trialEarnings)
    
    # intergrate data of 3 blocks
    thisTrialData = block2session(thisTrialData) 
    
    # survival analysis
    kmscResults = kmsc(thisTrialData, min(tMaxs), F, kmGrid)

    muWTW[sIdx] = kmscResults[['auc']]
    survCurve_[[sIdx]] = kmscResults$kmOnGrid
    stdWTW[[sIdx]] = kmscResults$stdWTW
    
    # WTW timecourse
    wtwtsResults = wtwTS(thisTrialData, tGrid, min(tMaxs), F)
    timeWTW_[[sIdx]] = wtwtsResults$timeWTW
    trialWTW_[[sIdx]] = wtwtsResults$trialWTW
  }
  # return outputs
  sumStats = data.frame(
    id = ids,
    condition = factor(hdrData$condition[hdrData$stress == "no_stress"], levels = c("HP", "LP")),
    nExcl = nExcl,
    totalEarnings = totalEarnings,
    muWTW = muWTW,
    stdWTW = stdWTW
  )
  outputs = list(
    sumStats = sumStats,
    survCurve_ = survCurve_,
    trialWTW_ = trialWTW_,
    timeWTW_ = timeWTW_ 
  )
  return(outputs)
}
