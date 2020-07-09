simExample = function(){
  set.seed(123)
  # default settings 
  modelName = "QL2"
  smallReward = 0 
  iti = 2
  
  # random seed
  set.seed(123)
  
  # create output directories
  dir.create("figures/simulation/")
  dir.create("figures/simulation/example")
  
  # load experiment parameters
  load('expParas.RData')
  
  # normative analysis 
  normResults = expSchematics(smallReward, iti, F)
  optimRewardRates = normResults$optimRewardRates
  optimWaitThresholds = normResults$optimWaitThresholds
  
  # load packages and sub functions 
  library("tidyverse")
  source("subFxs/plotThemes.R")
  source("subFxs/helpFxs.R") 
  source("subFxs/loadFxs.R") # 
  source("subFxs/analysisFxs.R") 
  
  
  # get the generative model 
  source(sprintf("subFxs/simModels/default/%s.R", modelName))
  simModel = get(modelName)
  paraNames = getParaNames(modelName)
  paraNames = factor(paraNames, levels = paraNames)
  nPara = length(paraNames)
  
  # num of repetitions 
  nRep = 10
  duration = 120 * 60
  periodBreaks = seq(0, 16 * 60, length.out = 5)
  
  HPSim_ = list()
  LPSim_ = list()
  HPauc_ = matrix(NA, nrow =  5, ncol = nRep)
  LPauc_ = matrix(NA, nrow =  5, ncol = nRep)
  for(i in 1 : nRep){
    # HP 
    tempt = simModel(c(0.05, 0.05, 5, 0.85, 6), "HP", duration, normResults)
    HPSim_[[i]] = tempt
    tempt$Qwaits_ = NULL
    tempt = as.data.frame(tempt)
    HPauc_[1 : 4,i] = sapply(1 : 4, function(x) kmsc(tempt[(i-1) * 4 * 60 <= tempt$sellTime & tempt$sellTime < i * 4 * 60,], min(delayMaxs), F, kmGrid)$auc)
    HPauc_[5,i] = kmsc(tempt[116 * 60 <= tempt$sellTime,], min(delayMaxs), F, kmGrid)$auc
    
    # LP
    tempt  = simModel(c(0.05, 0.05, 5, 0.85, 6), "LP", duration, normResults)
    LPSim_[[i]] = tempt
    tempt$Qwaits_ = NULL
    tempt = as.data.frame(tempt)
    LPauc_[1 : 4,i] = sapply(1 : 4, function(x) kmsc(tempt[(i-1) * 4 * 60 <= tempt$sellTime & tempt$sellTime < i * 4 * 60,], min(delayMaxs), F, kmGrid)$auc)
    LPauc_[5,i] = kmsc(tempt[116 * 60 <= tempt$sellTime,], min(delayMaxs), F, kmGrid)$auc
  }
  
  
  # asymptotic AUCs
  asyHPAUC = mean(HPauc_[5,])
  asyLPAUC = mean(LPauc_[5,])
  
  # plot AUCs in HP and LP
  data.frame(
    auc = c(as.vector(HPauc_[1:4,]), as.vector(LPauc_[1:4,])),
    condition = c(rep("HP", nRep * 4), rep("LP", nRep * 4)),
    t= rep(rep(seq(0,12, by = 4) + 2,  each = nRep), 2)
  ) %>% group_by(condition, t) %>%
    summarise(mu = mean(auc),
              se = sd(auc) / sqrt(length(auc)),
              min = mu - se,
              max = mu + se) %>% 
    ggplot(aes(t, mu, color = condition)) + geom_point() + geom_line() +
    myTheme + scale_color_manual(values = conditionColors) +
    scale_x_continuous(breaks = seq(0, 16, by = 4),
                       limits = c(0,16)) +
    xlab("Time (s)") +
    ylab("AUC (s)") + theme(legend.position = "none")  +
    ylim(c(0, 20)) + 
    geom_hline(yintercept = 2.2, color = conditionColors[2], linetype = "dashed") +
    geom_hline(yintercept = 20, color = conditionColors[1], linetype = "dashed") +
    geom_point(x = 16, y = asyHPAUC, color = conditionColors[1], shape = 8) + 
    geom_point(x = 16, y = asyLPAUC, color = conditionColors[2], shape = 8)
    ggsave("figures/simulation/example/auc.png", width = 3, height = 3)

  # Qwaits in HP and LP
  nCut = 5 # five obervations
  nHPSteps = nrow(HPSim_[[1]]$Qwaits_) - 1
  HPQwaits_ = matrix(NA, nrow = nHPSteps, ncol = nRep * nCut)
  for(i in 1 : nRep){
    junk = sapply(c(240, 480, 720), function(x) which.min(abs(HPSim_[[i]]$sellTime - x)))
    HPQwaits_[, (i-1) * nCut + 1 : nCut] = HPSim_[[i]]$Qwaits_[1 :nHPSteps ,c(1, junk,length(HPSim_[[i]]$sellTime))]
  }
  nLPSteps = nrow(LPSim_[[1]]$Qwaits_) - 1
  LPQwaits_ = matrix(NA, nrow = nLPSteps, ncol = nRep * nCut)
  for(i in 1 : nRep){
    junk = sapply(c(240, 480, 720), function(x) which.min(abs(LPSim_[[i]]$sellTime - x)))
    LPQwaits_[, (i-1) * nCut + 1 : nCut] = LPSim_[[i]]$Qwaits_[1 : nLPSteps,c(1, junk,length(LPSim_[[i]]$sellTime))]
  }
  HPdata = data.frame(
    Qwait = as.vector(HPQwaits_),
    id = rep(1 : nRep, each = nCut * nHPSteps),
    period = rep(rep(1 : nCut, nRep), each = nHPSteps),
    t = rep(1:nHPSteps, nRep * nCut),
    color = rep(rep(1 : nCut, nRep), each = nHPSteps),
    condition = "HP"
  ) 
  LPdata =  data.frame(
    Qwait = as.vector(LPQwaits_),
    id = rep(1 : nRep, each = nCut * nLPSteps),
    period = rep(rep(1 : nCut, nRep), each = nLPSteps),
    t = rep(1:nLPSteps, nRep * nCut),
    color = rep(rep(1 : nCut + nCut, nRep), each = nLPSteps),
    condition = "LP"
  ) 
  # asymptotic Q(wait, t)
  asyQwaitHP = HPdata %>% filter(period == 5) %>% group_by(t, condition) %>% summarise(mean(Qwait))
  asyQwaitLP = LPdata %>% filter(period == 5) %>% group_by(t, condition) %>% summarise(mean(Qwait))
  asyDf = rbind(asyQwaitHP, asyQwaitLP)
  # plot Qwaits 
  plotData = rbind(HPdata, LPdata)
  plotData %>% filter(period < 5) %>% group_by(period, t, condition) %>%
    summarise(mu = mean(Qwait), color = mean(color)) %>%
    ggplot(aes(t, mu, color = as.factor(color))) + geom_point() + 
    geom_line() + myTheme +
    facet_grid(~condition) + xlab("t (s)") + ylab("Q(wait, t)") +
    scale_color_manual(values = c("#c7e9c0", "#41ab5d", "#006d2c", "#00441b",
                                  "#bcbddc", "#807dba", "#6a51a3", "#3f007d")) +
    theme(legend.position = "none") + 
    geom_point(data = asyDf, aes(x = t, y = `mean(Qwait)`), inherit.aes = F) 
  ggsave("figures/simulation/example/Qwait.png", width = 4, height = 3) 
}





