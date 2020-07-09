# replicate behavioral data by sumulating with individual fitted parameters
expModelRep = function(modelName){
  set.seed(123)
  # create output directories
  dir.create("figures/expModelRep/")
  dir.create("genData/expModelRep/")
  dir.create(sprintf("figures/expModelRep/%s",modelName))
  
  # load experiment parameters
  load('expParas.RData')
  
  # load packages and sub functions 
  library("tidyverse")
  source("expSchematics.R")
  source("subFxs/plotThemes.R")
  source("subFxs/helpFxs.R") 
  source("subFxs/loadFxs.R") # 
  source("subFxs/analysisFxs.R") 
  source("subFxs/modelRep.R")
  
  
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

  ## check fit
  expPara = loadExpPara(paraNames, sprintf("genData/expModelFit/%s", modelName))
  passCheck = checkFit(paraNames, expPara)
  
  # compare willingness to wait (WTW) from empirical and replicated data
  ## WTW from empirical data 
  source("MFAnalysis.R")
  MFResults = MFAnalysis(isTrct = T)
  sumStats = MFResults[['sumStats']]
  muWTWEmp = sumStats$muWTW
  stdWTWEmp = sumStats$stdWTW
  
  ## replicate data
  repOutputs =  modelRep(trialData, ids, nRep, T)
  save(repOutputs, file = sprintf("genData/expModelRep/%s_trct.RData", modelName))
  plotData = data.frame(mu =  repOutputs$muWTWRep_mu, std = repOutputs$stdWTWRep_mu,
                        empMu = muWTWEmp, empStd = stdWTWEmp,
                          passCheck = passCheck, 
                        condition = sumStats$condition) %>% filter(passCheck)
  
  ## plot to compare average willingess to wait
    plotData %>%
    ggplot(aes(empMu, mu)) + 
    geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21) + 
    geom_abline(slope = 1, intercept = 0)  + 
    ylab("Model-generated (s)") + xlab("Observed (s)") + ggtitle(sprintf("AUC, n = %d", sum(passCheck))) +
    myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    scale_x_continuous(breaks = c(0, 20), limits = c(-1, 21)) + 
    scale_y_continuous(breaks = c(0, 20), limits = c(-1, 21)) +
    scale_color_manual(values = conditionColors) +
      theme(legend.position = "none")
  fileName = sprintf("figures/expModelRep/%s/muWTW_muWTWRep.eps", modelName)
  fileName = sprintf("figures/expModelRep/%s/muWTW_muWTWRep.png", modelName) 
  ggsave(filename = fileName,  width = 4, height = 4)

  ## plot to compare std willingess to wait
    plotData %>%
    ggplot(aes(empStd, std, shape = condition)) + 
    geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21)  + 
    geom_abline(slope = 1, intercept = 0) +
    ylab(expression(bold(paste("Model-generated (s"^2,")")))) +
      xlab(expression(bold(paste("Observed (s"^2,")")))) + ggtitle(sprintf("CIP, n = %d", sum(passCheck))) +
    myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) + 
    scale_x_continuous(breaks = c(0, 8), limits = c(0, 10)) + 
    scale_y_continuous(breaks = c(0, 8), limits = c(0, 10)) +
      scale_color_manual(values = conditionColors) +
      theme(legend.position = "none")
  fileName = sprintf("figures/expModelRep/%s/stdWTW_stdWTWRep.eps", modelName)
  fileName = sprintf("figures/expModelRep/%s/stdWTW_stdWTWRep.png", modelName)
  ggsave(filename = fileName,  width = 4, height = 4)
  
  # replicate data again 
  repOutputs =  modelRep(trialData, ids, nRep, F)
  save(repOutputs, file = sprintf("genData/expModelRep/%s.RData", modelName))
  timeWTW_ = repOutputs$timeWTW_
  
  # plot to check time course of learning
  data.frame(
    wtw = as.vector(timeWTW_),
    t = rep(tGrid, nSub),
    condition = rep(sumStats$condition, each = length(tGrid))
  ) %>% group_by(condition, t) %>%
    summarise(mu = mean(wtw, na.rm = F),
              se = sd(wtw, na.rm = F) / sqrt(length(wtw)),
              min = mu - se,
              max = mu + se) %>%
    ggplot(aes(t, mu)) + geom_line(aes(color = condition)) +
    geom_ribbon(aes(ymin=min, ymax=max, fill = condition), alpha = 0.5) +
    scale_color_manual(values = conditionColors) +
    scale_fill_manual(values = conditionColors) + myTheme +
    theme(legend.position = "none") +
    ylab("WTW (s)") + xlab("Task time (min)")  +
    scale_x_continuous(breaks = 0:3 * 420, labels = seq(0, blockSec * nBlock, blockSec) / 60) + 
    theme(legend.position = "none") + ylim(0, 16) 
  fileName = sprintf("figures/expModelRep/%s/timecourse.eps", modelName)
  ggsave(filename = fileName,  width = 5, height = 3)
  fileName = sprintf("figures/expModelRep/%s/timecourse.png", modelName)
  ggsave(filename = fileName,  width = 5, height = 3)
  
  # plot example participants 
  # model generated 
  repTrialData = repOutputs$repTrialData
  repNo = repOutputs$repNo
  sIdx = 1
  thisRepTrialData = repTrialData[[repNo[1, sIdx]]]
  
  thisRepTrialData = data.frame(thisRepTrialData[1:6])
  trialPlots(thisRepTrialData) +  
    ggtitle("Model-generated") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
  ggsave(sprintf("figures/expModelRep/%s/sim1.eps", modelName), width = 6, height = 4)
  ggsave(sprintf("figures/expModelRep/%s/sim1.png", modelName), width = 6, height = 4)
  
  # observed
  thisTrialData = trialData[[ids[sIdx]]]
  thisTrialData  = thisTrialData %>% filter(trialStartTime <=  blockSec - max(delayMaxs))
  thisTrialData = block2session(thisTrialData)
  trialPlots(thisTrialData) + ggtitle("Observed") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
  ggsave(sprintf("figures/expModelRep/%s/exp1.eps", modelName), width = 6, height = 4)
  ggsave(sprintf("figures/expModelRep/%s/exp1.png", modelName),, width = 6, height = 4)
}

