# replicate behavioral data by sumulating with individual fitted parameters
expModelRepitation = function(modelName){
  # create output directories
  dir.create("figures/expModelRep/")
  dir.create(sprintf("figures/expModelRep/%s",modelName))
  
  # load experiment parameters
  load('expParas.RData')
  
  # load packages and sub functions 
  library("ggplot2") 
  library("dplyr")
  library("tidyr")
  source("subFxs/plotThemes.R")
  source("subFxs/helpFxs.R") 
  source("subFxs/loadFxs.R") # 
  source("subFxs/analysisFxs.R") 
  
  
  # get the generative model 
  source(sprintf("subFxs/%s.R", modelName))
  gnrModel = get(modelName)
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  
  # num of repetitions 
  nRep = 10
  
  # load empirical data 
  allData = loadAllData()
  hdrData = allData$hdrData 
  trialData = allData$trialData       
  ids = hdrData$id[hdrData$stress == "no_stress"]
  nSub = length(ids)

  # initialize outputs
  repTrialData = vector(length = nSub * nRep, mode ='list')
  repNo = matrix(1 : (nSub * nRep), nrow = nRep, ncol = nSub)
  
  # loop over participants
  for(sIdx in 1 : nSub){
    # prepare empirical data 
    id = ids[[sIdx]]
    thisTrialData = trialData[[id]] 
    # excluded trials at the end of blocks 
    excluedTrials = which(thisTrialData$trialStartTime > (blockSec - max(tMaxs)))
    thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excluedTrials,]
    # load individually fitted paramters 
    fitSummary = read.table(sprintf("genData/expModelFit/%s/s%s_summary.txt",  modelName, id),sep = ",", row.names = NULL)
    paras =  fitSummary[1 : nPara,1]
    # simulate nRep times
    for(rIdx in 1 : nRep){
      tempt = gnrModel(paras, thisTrialData$condition, thisTrialData$scheduledWait)
      repTrialData[[repNo[rIdx, sIdx]]] = tempt
    }
  }
  
  # compare willingness to wait (WTW) from empirical and replicated data
  ## WTW from empirical data 
  source("MFAnalysis.R")
  MFResults = MFAnalysis(isTrct = T)
  sumStats = MFResults[['sumStats']]
  muWTWEmp = sumStats$muWTW
  stdWTWEmp = sumStats$stdWTW
  ## WTW from empirical data 
  muWTWRep_ = matrix(NA, nrow = nRep , ncol = nSub)
  stdWTWRep_ = matrix(NA, nrow = nRep, ncol = nSub)
  for(sIdx in 1 : nSub){
    id = ids[[sIdx]]
    for(rIdx in 1 : nRep){
      thisRepTrialData = repTrialData[[repNo[rIdx, sIdx]]]
      kmscResults = kmsc(thisRepTrialData, min(tMaxs), F, kmGrid)
      muWTWRep_[rIdx,sIdx] = kmscResults$auc
      stdWTWRep_[rIdx, sIdx] = kmscResults$stdWTW
    }
  }
  ## summarise WTW across simulations for replicated data 
  muWTWRep_mu = apply(muWTWRep_, MARGIN = 2, mean) # mean of average willingness to wait
  muWTWRep_std = apply(muWTWRep_, MARGIN = 2, sd) # std of average willingess to wait
  stdWTWRep_mu = apply(stdWTWRep_, MARGIN = 2, mean) # mean of std willingness to wait
  stdWTWRep_std = apply(stdWTWRep_, MARGIN = 2, sd) # std of std willingess to wait
  ## check fit
  expPara = loadExpPara(paraNames, sprintf("genData/expModelFit/%s", modelName))
  passCheck = checkFit(paraNames, expPara)
  
  ## plot to compare average willingess to wait
  data.frame(mu =  muWTWRep_mu, std = muWTWRep_std,
             empMu = muWTWEmp, passCheck,
             condition = sumStats$condition) %>%
    mutate(min = mu - std, max = mu + std) %>%
    filter(passCheck == T) %>%
    ggplot(aes(empMu, mu)) +  geom_errorbar(aes(ymin = min, ymax = max), color = "grey") +
    geom_point(size = 2) + facet_grid(~condition) + 
    geom_abline(slope = 1, intercept = 0) + xlim(c(-2, 22)) + ylim(c(-2, 22)) +
    ylab("Model-generated (s)") + xlab("Observed (s)") + ggtitle(sprintf("Average WTW, n = %d", sum(passCheck))) +
    myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5))
  fileName = sprintf("figures/expModelRep/%s/AUC_AUCRep.png", modelName)
  ggsave(filename = fileName,  width = 6, height = 4)


  
}
