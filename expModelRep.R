# replicate behavioral data by sumulating with individual fitted parameters
expModelRep = function(modelName){
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
  source(sprintf("subFxs/gnrModels/%s.R", modelName))
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
  
  
  ## save Results
  plotData = data.frame(mu =  muWTWRep_mu, std = stdWTWRep_mu,
                        empMu = muWTWEmp, empStd = stdWTWEmp,
                          passCheck = rep(passCheck, each = 2), 
                        condition = sumStats$condition) %>% filter(passCheck == T)
  save(plotData, file = "genData/rscRep.RData")
  
  cor.test(plotData$mu, plotData$empMu, methods = "spearman")
  ## plot to compare average willingess to wait
    plotData %>%
    ggplot(aes(empMu, mu, shape = condition)) + 
    geom_point(size = 2, color = themeColor, fill = "pink", stroke = 1) + 
    geom_abline(slope = 1, intercept = 0) +
    scale_shape_manual(values = c(21, 23)) + 
    ylab("Model-generated (s)") + xlab("Observed (s)") + ggtitle(sprintf("Average WTW, n = %d", sum(passCheck))) +
    myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    scale_x_continuous(breaks = c(0, 20), limits = c(-1, 21)) + 
    scale_y_continuous(breaks = c(0, 20), limits = c(-1, 21))
  
  fileName = sprintf("figures/expModelRep/%s/muWTW_muWTWRep.eps", modelName) 
  ggsave(filename = fileName,  width = 4, height = 4)
  
  ## plot to compare std willingess to wait
    plotData %>%
    ggplot(aes(empStd, std, shape = condition)) + 
    geom_point(size = 2, color = themeColor, fill = "pink", stroke = 1)  + 
    geom_abline(slope = 1, intercept = 0) +
    ylab(expression(bold(paste("Model-generated (s"^2,")")))) + xlab(expression(bold(paste("Observed (s"^2,")")))) + ggtitle(sprintf("Std WTW, n = %d", sum(passCheck))) +
    myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    scale_shape_manual(values = c(21, 23)) + 
    scale_x_continuous(breaks = c(0, 8), limits = c(0, 10)) + 
    scale_y_continuous(breaks = c(0, 8), limits = c(0, 10))
  fileName = sprintf("figures/expModelRep/%s/stdWTW_stdWTWRep.eps", modelName)
  ggsave(filename = fileName,  width = 4, height = 4)
  
  
  # plot example participants 
  # model generated 
  sIdx = 1
  thisRepTrialData = repTrialData[[repNo[1, sIdx]]]
  thisRepTrialData$timeWaited = apply(sapply(1 : nRep, FUN = function(i) 
    repTrialData[[repNo[i, sIdx]]]$timeWaited), MARGIN = 1, FUN = mean)
  thisRepTrialData = data.frame(thisRepTrialData[1:6])
  trialPlots(thisRepTrialData) 
  ggsave("figures/expModelRep/sim1.eps", width = 6, height = 4)
  
  # 
  thisTrialData = trialData[[sIdx]]
  thisTrialData  = thisTrialData %>% filter(trialStartTime <=  blockSec - max(tMaxs))
  thisTrialData = block2session(thisTrialData)
  trialPlots(thisTrialData)
  ggsave("figures/expModelRep/exp1.eps", width = 6, height = 4)
}
