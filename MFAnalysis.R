# model-free analysis

# inputs:
# isTrct : logical variable determining whether the last portion in each block is truncated 

# outputs:
# nExcludes : [nSubx1 int] # total number of excluded trials 
# muWTWs : [nSubx1 num] # average willingness to wait (WTW), measured by area under the Kaplan-Meier survival curve
# stdWTWs : [nSubx1 num] # standard deviation of WTW, measured in Kaplan-Meier survival analysis
# totalEarnings_s :  [nSubx1 num] 
# timeWTW_ : list(nSubx1) # wtw timecourse, each element is a vector
# trialWTW_ : list(nSubx1) # trial-wise WTW, each element is a vector
# survCurve_ : list(nSubx1) # Kaplan-Meier survival curve, each element is a vector


# load libraries
source('subFxs/loadFxs.R') 
source('subFxs/analysisFxs.R') 
library('dplyr')

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

# analyses parameters
tGrid = seq(0, blockSec * nBlock, by = 3) # time grid for wtw time courses
kmGrid = seq(0, min(tMax_), by = 0.1) # time grid for Kaplan-Meier survival curves

# initialize output variables 
nExcls = numeric(length = nSub)
muWTWs = numeric(length = nSub) 
stdWTWs = numeric(length = nSub) 
totalEarnings_s =  numeric(length = nSub) 
timeWTW_ = vector(mode = "list", length = nSub) 
trialWTW_ = vector(mode = "list", length = nSub) 
survCurve_ = vector(mode = "list", length = nSub) 

# if True, perfom all analyses and plot inividual figures one by one and press any key to proceed
# if False, perform all analyses without plotting
plotTrialwiseData = F # plot trial-wise data
plotKMSC = F # plot Kaplan-Meier survival curves
plotWTW = F # plot williness to wait (WTW) time courses

# loop over inidviduals
for (sIdx in 1 : nSub) {
  # load (and truncate) trialData for this individual
  id = ids[sIdx]
  thisTrialData = trialData[[id]]
  if(isTrct){
    trctLine = blockSec - max(tMaxs)
    # truncate trials completed after tractline in each block
    nExcls[sIdx] = sum(thisTrialData$sellTime > trctLine))
    thisTrialData = thisTrialData %>% filter(sellTime <=  trctLine )
  }else{
    nExcls[sIdx] = 0
  }
  # calcualte totalEarnings
  totalEarnings_s[sIdx] =  sum(thisTrialData$trialEarnings)
  
  #  
  plotTitle = sprintf('Subject %s, Cond %s',id, unique(thisTrialData$condition))
  # in
  thisTrialData = block2session(thisTrialData) 

  # plot trial-wise data
  if (plotTrialwiseData) {
    trialPlots(thisTrialData,label)
    readline(prompt = paste( plotTitle, '(hit ENTER to continue)'))
    graphics.off()
  }
  
  # survival analysis
  kmscResults = kmsc(thisTrialData, min(tMaxs), label, plotKMSC, kmGrid)
  muWTWs[sIdx] = kmscResults[['auc']]
  survCurves_[[sIdx]] = kmscResults$kmOnGrid
  if (plotKMSC) {
    readline(prompt = paste( plotTitle, '(hit ENTER to continue)'))
    graphics.off()
  }
  stdWTWs[[sIdx]] = kmscResults$stdWd

  
  # WTW timecourse
  wtwtsResults = wtwTS(thisTrialData, tGrid, min(tMaxs), label, plotWTW)
  timeWTW_[[sIdx]] = wtwtsResults$timeWTW
  trialWTW_[[sIdx]] = wtwtsResults$trialWTW
  if (plotWTW) {
    readline(prompt = paste('subject',id, '(hit ENTER to continue)'))
    graphics.off()
  }
}