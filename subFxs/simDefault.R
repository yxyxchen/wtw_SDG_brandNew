simDefault = function(){
  # default settings 
  modelName = "QL2"
  smallReward = 0
  iti = 2
  
  # random seed
  set.seed(123)
  
  # create output directories
  dir.create("figures/simulation/")
  dir.create("figures/simulation/default")
  
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
  paraLabels = c("alpha[r]","alpha[u]", "tau", "gamma", "eta")
  nPara = length(paraNames)
  
  # num of repetitions 
  nRep = 10
  duration = 16 * 60
  nPeriod = 4
  periodBreaks = seq(0, duration, length.out = nPeriod + 1)
  ######################### simulate with multiple parameter combinations #####################
  # parameter combinations 
  nCut = 4
  paraSampleHPs = data.frame(
    alphaR = seq(0.01, 0.05, length.out = nCut),
    alphaU = seq(0.01, 0.05, length.out = nCut),
    tau = exp(seq(log(0.5), log(8), length.out = nCut)),
    gamma = seq(0.75, 0.95, length.out = nCut),
    prior = seq(0, 1, length.out = nCut)
  )
  paraSampleLPs = data.frame(
    alphaR = seq(0.01, 0.05, length.out = nCut),
    alphaU = seq(0.01, 0.05, length.out = nCut),
    tau = exp(seq(log(0.5), log(8), length.out = nCut)),
    gamma = seq(0.75, 0.95, length.out = nCut),
    prior = seq(2, 6, length.out = nCut)
  )
  
  
  for(condition in conditions){
    if(condition == 'HP'){
      paraSamples = paraSampleHPs
    }else{
      paraSamples = paraSampleLPs
    }
    paraCombs = apply(paraSamples[3, ], 2, rep, nCut * nPara)
    for(pIdx in 1 : nPara){
      paraCombs[(nCut * (pIdx - 1) + 1) : (nCut * pIdx) ,pIdx] = paraSamples[,pIdx] 
    }
    nComb = nrow(paraCombs)
    
    # initialize outputs
    auc_ = matrix(NA, nrow = nComb, ncol = nRep)
    periodAuc_ = matrix(NA, nrow = nComb * nPeriod, ncol = nRep)
    # loop
    for(i in 1 : nComb){
      paras = paraCombs[i, ]
      for(j in 1 : nRep){
        tempt = simModel(as.numeric(paras), condition, duration, normResults)
        kmscResults = kmsc(tempt, min(delayMaxs), F, kmGrid)
        auc_[i, j] = kmscResults$auc
        tempt$Qwaits_ = NULL
        tempt = do.call(cbind.data.frame, tempt)
        periodAuc_[((i - 1) * nPeriod + 1) : (i * nPeriod), j] = 
          sapply(1 : nPeriod, function(x) kmsc(tempt[periodBreaks[x] <= tempt$sellTime & tempt$sellTime < periodBreaks[x+1],], min(delayMaxs), F, kmGrid)$auc)
      }
    }
    df = data.frame(
      auc = apply(periodAuc_, 1, mean),
      period = rep(1 : nPeriod, nComb),
      t = rep((periodBreaks[1 : nPeriod] + periodBreaks[2 : (nPeriod+1)] )/ 120, nComb),
      paraName = rep(paraNames, each = nCut * nPeriod),
      paraLabel = rep(paraLabels, each = nCut * nPeriod),
      paraRank = factor(rep(rep(1 : nCut, each = nPeriod), nPara)),
      condition = condition
    ) 
    if(condition == "HP"){
      HPdf = df
    }else{
      LPdf = df
    }
  }
  
  # plot
  cutValues = c(
    "#bdbdbd",
    "#737373",
    "#525252",
    "#252525"
  )
  optimDf = data.frame(
    t = rep(seq(0, 16, by = 4), 10),
    condition = rep(conditions, each = 5 * 5),
    optimThreshold = rep(unlist(optimWaitThresholds), each = 5 * 5),
    paraLabel= rep(rep(paraLabels, each = 5), 2)
  )
  rbind(HPdf, LPdf) %>% 
    mutate(paraLabel = factor(paraLabel, levels = paraLabels)) %>%
    ggplot(aes(t, auc, color = paraRank)) +
    geom_line() + geom_point() +
    facet_grid(condition~paraLabel, labeller = label_parsed) +
    geom_line(data = optimDf, aes(x = t, y = optimThreshold),
              inherit.aes = F, linetype = "dashed", color = rep(conditionColors, each = 25)) +
    scale_color_manual(values = cutValues) +
    myTheme + xlab("Task time (min)") + ylab("AUC (s)") +
    scale_x_continuous(breaks = periodBreaks / 60, limits = c(0, duration / 60)) +
    ylim(c(0, 20)) 
  ggsave("figures/simulation/default/paraEffect.png", width = 7, height = 4)
}





