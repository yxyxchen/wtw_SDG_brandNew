# this script plots delay distributions and reward rates in two environments
expSchematics = function(smallReward, iti, isPlot){
  # load libararies
  library("tidyverse")
  library(latex2exp)
  source('subFxs/plotThemes.R')
  if(!dir.exists('figures/expSchematics')){
    dir.create('figures/expSchematics')
  }
  
  # load experiment parameters 
  load("expParas.RData")
  
  # for display purposes, all variables on the continous time scale
  # are discretized into 0.1 second time bins
  bin = 0.1 # width of a time bin
  time = list(
    HP = seq(bin, delayMaxs[1], by = bin),
    LP = seq(bin, delayMaxs[2], by = bin)
  ) 
  
  # delay CDFS
  ### delays in HP follow unif(0, 20)
  ### delays in LP follow pareto(mu = 0, k =4, sigma = 2), truncated at 40
  HP = 1 / length(time[['HP']]) * (1 : length(time[['HP']]))
  mu = 0; k = 4; sigma = 2; # parameters for the pareto distribution
  LP = 1 - (1 + k * (time[['LP']]- mu) / sigma) ^ (-1 / k)
  LP[length(LP)] = 1 # truncated at 40
  rewardDelayCDFs = list(
    HP = HP,
    LP = LP
  )
  
  # delay PDFs
  HP = diff(c(0, rewardDelayCDFs$HP))
  LP = diff(c(0, rewardDelayCDFs$LP))
  rewardDelayPDFs = list(
    "HP" = HP,
    "LP" = LP
  )
  
  # average waiting durations given different policies
  # Here we assume rewards occur at the middle of each time bin
  HP = cumsum((time[['HP']] - 0.5 * bin) * rewardDelayPDFs$HP) / cumsum(rewardDelayPDFs$HP)
  LP = cumsum((time[['LP']] - 0.5 * bin) * rewardDelayPDFs$LP) / cumsum(rewardDelayPDFs$LP)
  meanRewardDelays = list('HP' = HP, 'LP' = LP)
  
  # rewardRates given different policies
  ## might be different from the values used in expParas.R, 
  ## which are calcuated with a higher temporal resoluation
  HP = (tokenValue * rewardDelayCDFs$HP + smallReward * (1 - rewardDelayCDFs$HP)) /
    ((meanRewardDelays$HP * rewardDelayCDFs$HP + time[['HP']] * (1 - rewardDelayCDFs$HP)) + iti)
  LP = (tokenValue * rewardDelayCDFs$LP + smallReward * (1 - rewardDelayCDFs$LP))/
    ((meanRewardDelays$LP * rewardDelayCDFs$LP + time[['LP']] * (1 - rewardDelayCDFs$LP)) + iti)
  rewardRates = list('HP' = HP, 'LP' = LP)
  
  # optimal raward rates and optimal policies
  optimWaitThresholds = list()
  optimWaitThresholds$HP = time$HP[which.max(HP)]
  optimWaitThresholds$LP = time$LP[which.max(LP)]
  optimRewardRates = list()
  optimRewardRates$HP = max(HP)
  optimRewardRates$LP = max(LP)
  
  # calculate subjective value as a function of elapsed time 
  subjectValues = list()
  for(cIdx in 1 : 2){
    condition = conditions[cIdx]
    delayMax = delayMaxs[cIdx]
    pdf = rewardDelayPDFs[[cIdx]]
    cdf = rewardDelayCDFs[[cIdx]]
    thisTime = time[[cIdx]]
    Rstar = optimRewardRates[[cIdx]]
    ts = seq(0, delayMax, by = 0.1) # elapsed time
    
    # initialize 
    thisSubjectValues = vector(length = length(ts))
    Tstars = vector(length = length(ts))
    # loop over different elapsed time
    for(i in 1 : length(ts)){
      t = ts[i] # this elapsed time
      trctTime = thisTime[thisTime > t]
      
      if(t == delayMax){
        Tstar = t
        gt_max = tokenValue 
      }else{
        # loop over different waiting policies 
        Tstar = t
        gt_max = -100
        for(T in seq(t, delayMax, by = 0.1)){
          trctPDF = pdf[thisTime > t] / sum(pdf[thisTime > t])
          at = tokenValue * sum(trctPDF[trctTime <= T])
          trctPDF[trctTime == T] = trctPDF[trctTime == T]
          bt = sum((trctTime[trctTime <= T] - 0.5 * bin - t) * trctPDF[trctTime <= T]) + 
            + (T - t) * sum(trctPDF[trctTime > T]) 
          gt = at - bt * Rstar
          if(gt > gt_max){
            gt_max = gt
            Tstar = T
          }
        }
      }
      
      thisSubjectValues[i] = gt_max 
      Tstars[i] = Tstar
    }
    subjectValues[[condition]] = thisSubjectValues
  }
  
  if(isPlot){
    # plot CDFs 
    ## here we extend the HP CDF to 32s for display purposes
    data.frame(CDF = c(0,c(rewardDelayCDFs$HP, rep(1, length(time$LP) - length(time$HP))), 0, rewardDelayCDFs$LP),
               time = c(0, time$LP, 0, time$LP),
               condition =  rep(c('HP', 'LP'), c(length(time$LP) + 1, length(time$LP) + 1))) %>%
      ggplot(aes(time, CDF)) + geom_line(size = 3, aes(color = condition)) + facet_grid(~condition) +
      ylim(c(0,1)) + scale_y_continuous(breaks = c(0,0.5,1)) + 
      scale_x_continuous(breaks = c(0, max(delayMaxs)/ 2, max(delayMaxs)), labels = c("0", max(delayMaxs)/2, max(delayMaxs)),
                         limits = c(0, max(delayMaxs) * 1.1)) + 
      myTheme + xlab('Delay duration (s)') + ylab('CDF') + 
      theme(plot.title = element_text(hjust = 0.5, color = themeColor)) + 
      annotate("text", x = max(delayMaxs)/ 2, y = 0.4, label = "+10¢", size = 6) +
      scale_color_manual(values = conditionColors) +
      theme(legend.position = "none")
    
    ggsave('figures/expSchematics/CDF.eps', width =4, height = 3)
    ggsave('figures/expSchematics/CDF.png', width =4, height = 3)
    
    # plot reward rates
    optimData = data.frame(condition = c("HP", "LP"), waitThreshold = as.double(optimWaitThresholds))
    data.frame(rewardRate = c(0, rewardRates[[1]], 0, rewardRates[[2]]),
               time = c(0, time[[1]], 0, time[[2]]),
               condition = rep(c("HP", "LP"), c(length(time$HP) + 1, length(time$LP) + 1))) %>%
      ggplot(aes(time, rewardRate)) +
      geom_line(size = 2, aes(color = condition))  + myTheme + 
      ylab(expression(bold("Reward rate (¢ s"^"-1"*")"))) + xlab("Waiting policy (s)") + 
      theme(plot.title = element_text(hjust = 0.5, color = themeColor)) + 
      scale_y_continuous(breaks = c(0, 0.4, 0.8, 1.2), limits = c(0, 1.5)) +
      scale_x_continuous(breaks = c(0, max(delayMaxs)/ 2, max(delayMaxs)),
                         limits = c(0, max(delayMaxs) * 1.1)) +
      scale_color_manual(values = conditionColors) +
      theme(legend.position = "none") + facet_grid(~condition)
    ggsave("figures/expSchematics/reward_rate.eps", width = 4, height = 3)
    ggsave("figures/expSchematics/reward_rate.png", width = 4, height = 3)
    
    # plot prior belief
    eta = 20
    gamma = 0.6
    V0 = (optimRewardRates$HP + optimRewardRates$LP) * 0.5 / (1/6) 
    Qquit = V0 * gamma
    ts = seq(0, max(delayMaxs), by = 0.5)
    Qwaits = (eta - ts) * 0.1 + Qquit
    data = data.frame(x = ts, y = Qwaits, Qquit = rep(Qquit, length = length(ts))) 
    data$z = ifelse(data$y > Qquit, data$y, Qquit)
    
    ggplot(data, aes(x, y)) +
      geom_line(size = 2, color = "#2166ac") + geom_line(aes(x, Qquit), size = 2, color = "#b2182b") +
      geom_segment(aes(x = eta, xend = eta, y = 0, yend = Qquit), linetype = "dashed") +
      scale_x_continuous(breaks = c(20), labels = c(expression(10~eta)), limits = c(0, max(delayMaxs))) +
      scale_y_continuous(breaks = NULL) + xlab("t") + ylab("Initial Value") +
      myTheme +
      theme(text = element_text(size=20))
    ggsave('figures/expSchematics/prior.eps', width = 3, height = 3)
    ggsave('figures/expSchematics/prior.png', width = 3, height = 3)
    
    
    # plot subjective value of waiting 
    data.frame(
      value =  c(subjectValues$HP, subjectValues$LP),
      t = c(seq(0, delayMaxs[1], by = 0.1), seq(0, delayMaxs[2], by = 0.1)),
      condition = rep(conditions, c(length(subjectValues$HP), length(subjectValues$LP))))%>%
      ggplot(aes(t, value)) +
      geom_line(aes(color = condition), size = 2) +
      myTheme+
      scale_color_manual(values = conditionColors) +
      scale_linetype_manual(values = c(1, 2)) +
      xlab("Elapsed time (s)") + ylab("Subjective value (¢)")  + 
      theme(legend.position = "none") + facet_grid(~condition)
    ggsave("figures/expSchematics/subjective.eps", width = 4, height = 3)
    ggsave("figures/expSchematics/subjective.png", width = 4, height = 3)      
  }
  
  # return outputs 
  outputs = list(
    "optimWaitThresholds" = optimWaitThresholds,
    "optimRewardRates" = optimRewardRates,
    "rewardRates" = rewardRates,
    "subjectValues" = subjectValues,
    "time" = time
  )
  return(outputs)
}


  


