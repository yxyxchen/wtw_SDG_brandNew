# fit a reinforcement learning model for a single participant in Rstan 

# inputs:
# thisTrialData: behavioral data for this participant
# fileName: the name of the output file
# modelName: the name of   the model 
# paraNames: parameters for the model
# model: the Bayesian model 
# config: a list containing the Rstab configuration 
modelFitSingle = function(id, thisTrialData, modelName, paraNames, model, config, outputFile){
    # load experiment paras
    load('expParas.RData')
  
    # parse the stan configuration
    nChain = config[['nChain']] # number of MCMC chains
    nIter = config[['nIter']] # number of iterations on each chain
    controlList = list(adapt_delta = config[['adapt_delta']],
                       max_treedepth = config[['max_treedepth']] )
    warningFile = config[['warningFile']] # output file for stan warnings and errors
    
    # duration of one time step (namely one temporal state) 
    stepSec = 1
    
    # prepare inputs for fitting the model
    condition = unique(thisTrialData$condition)
    ## maximal number of steps in a trial
    nStepMax = ifelse(condition == "HP", tMaxs[1] / stepSec, tMaxs[2] / stepSec)
    ## ensure timeWaited = scheduledWait on rewarded trials
    thisTrialData = within(thisTrialData, {timeWaited[trialEarnings!= 0] = scheduledWait[trialEarnings!= 0]})
    ## terminal state in each trial
    ## Noticeably, ceiling(timeWaited / stepSec) gives the num of steps
    ## and adding 1 upon it gives the terminal state, namely the state following the final action. 
    Ts = with(thisTrialData, {round(ceiling(timeWaited / stepSec) + 1)}) 
    ## orgianze inputs into a list
    inputs <- list(
      iti = iti,
      stepSec = stepSec,
      nStepMax = nStepMax,
      N = length(thisTrialData$blockNum), # number of trials
      Rs = thisTrialData$trialEarnings, # rewards on each trial
      Ts = Ts)
    if(modelName %in% c("QL1", "QL2")){
      ## in Q-learning, the initial value of the iti state is proportional to 
      ## the discounted total rewards averaged across two conditions
      ## discount factor for one step is 0.85
      VitiIni = 0.9 * mean(unlist(optimRewardRates) * stepSec / (1 - 0.85))
      inputs$VitiIni = VitiIni
    }else{
      ## in R-learning, the initial reward rate is proportional to
      ## the optimal reward rates averaged across two conditions
      reRateIni = 0.9 * mean(unlist(optimRewardRates)) * stepSec;
      inputs$reRateIni = reRateIni     
    }
   
   # fit the model
    withCallingHandlers({
      fit = sampling(object = model, data = inputs, cores = 1, chains = nChain,
                     iter = nIter, control = controlList) 
      print(sprintf("Finish %s !",id))
      write(sprintf("Finish %s !", id), warningFile, append = T, sep = "\n")
    }, warning = function(w){
      warnText = paste(modelName, id, w)
      write(warnText, warningFile, append = T, sep = "\n")
    })
            
  # extract posterior samples
  samples = fit %>%
    rstan::extract(permuted = F, pars = c(paraNames, "LL_all"))
  # save posterior samples
  samples = samples %>% adply(2, function(x) x) %>% dplyr::select(-chains) 
  write.table(matrix(unlist(samples), ncol = length(paraNames) + 1), file = sprintf("%s.txt", outputFile), sep = ",",
              col.names = F, row.names=FALSE) 
  
  # calculate WAIC and Efficient approximate leave-one-out cross-validation (LOO)
  log_lik = extract_log_lik(fit) 
  WAIC = waic(log_lik)
  LOO = loo(log_lik)
  save("WAIC", "LOO", file = sprintf("%s_waic.RData", outputFile))
  
  # summarise posterior parameters and LL_all
  fitSummary <- summary(fit, pars = c(paraNames, "LL_all"), use_cache = F)$summary
  write.table(matrix(fitSummary, nrow = length(paraNames) + 1), file = sprintf("%s_summary.txt", outputFile), 
              sep = ",", col.names = F, row.names=FALSE)
}

modelFittingCV = function(thisTrialData, fileName, paraNames, model, modelName){
  #
  load("wtwSettings.RData")
  # simulation parameters
  nChain = 4
  nIter = 5000
  
  # determine wIni
  # since the participants' initial strategies are unlikely optimal
  # we multiple the optimal opportunity cost by subOptimalRatio
  subOptimalRatio = 0.9 
  if(any(paraNames  == "gamma") || modelName == "BL" ){
    wIni = (5/6 + 0.93) / 2 * stepDuration / (1 - 0.9) * subOptimalRatio
  }else{
    wIni = (5/6 + 0.93) / 2 * stepDuration * subOptimalRatio
  }
  
  # prepare input
  timeWaited = thisTrialData$timeWaited
  scheduledWait = thisTrialData$scheduledWait
  trialEarnings = thisTrialData$trialEarnings
  timeWaited[trialEarnings > 0] = scheduledWait[trialEarnings > 0]
  cond = unique(thisTrialData$condition)
  tMax = ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  nTimeSteps = tMax / stepDuration
  Ts = round(ceiling(timeWaited / stepDuration) + 1)
  data_list <- list(wIni = wIni,
                    nTimeSteps = nTimeSteps,
                    N = length(timeWaited),
                    trialEarnings = trialEarnings,
                    Ts = Ts,
                    iti = iti,
                    stepDuration = stepDuration)
  rm(last.warning)
  fit = sampling(object = model, data = data_list, cores = 1, chains = nChain,
                 iter = nIter) 
  if(exists("last.warning")){
    fileNameShort = str_extract(fileName, pattern = "s[0-9]*_f[0-9]*")
    sapply(1 : length(last.warning), function(i) print(paste(modelName, fileNameShort, paste(as.character(names(last.warning)[i], collapse  = " ")))))
  }

  # save
  fitSummary <- summary(fit,pars = c(paraNames, "LL_all"), use_cache = F)$summary
  write.table(matrix(fitSummary, nrow = length(paraNames) + 1), file = sprintf("%s_summary.txt", fileName),  sep = ",",
              col.names = F, row.names=FALSE)
}


modelFittingdb = function(thisTrialData, fileName, paraNames, model, modelName,nPara, low, up){
  #
  load("wtwSettings.RData")
  # simulation parameters
  nChain = 4
  nIter = 5000
  
  # determine wIni
  # since the participants' initial strategies are unlikely optimal
  # we multiple the optimal opportunity cost by subOptimalRatio
  subOptimalRatio = 0.9 
  if(any(paraNames  == "gamma") || modelName == "BL" ){
    wIni = (5/6 + 0.93) / 2 * stepDuration / (1 - 0.9) * subOptimalRatio
  }else{
    wIni = (5/6 + 0.93) / 2 * stepDuration * subOptimalRatio
  }
  
  # prepare input
  timeWaited = thisTrialData$timeWaited
  scheduledWait = thisTrialData$scheduledWait
  trialEarnings = thisTrialData$trialEarnings
  timeWaited[trialEarnings > 0] = scheduledWait[trialEarnings > 0]
  cond = unique(thisTrialData$condition)
  tMax = ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  nTimeSteps = tMax / stepDuration
  Ts = round(ceiling(timeWaited / stepDuration) + 1)
  data_list <- list(wIni = wIni,
                    nTimeSteps = nTimeSteps,
                    nPara = nPara,
                    N = length(timeWaited),
                    trialEarnings = trialEarnings,
                    Ts = Ts,
                    iti = iti,
                    stepDuration = stepDuration,
                    low =low,
                    up = up)
  fit = sampling(object = model, data = data_list, cores = 1, chains = nChain,
                 iter = nIter) 
  # extract parameters
  extractedPara = fit %>%
    rstan::extract(permuted = F, pars = c(paraNames, "LL_all"))
  # save sampling sequences
  tempt = extractedPara %>%
    adply(2, function(x) x) %>%  # change arrays into 2-d dataframe 
    dplyr::select(-chains) 
  write.table(matrix(unlist(tempt), ncol = length(paraNames) + 1), file = sprintf("%s.txt", fileName), sep = ",",
              col.names = F, row.names=FALSE) 
  # calculate and save WAIC
  log_lik = extract_log_lik(fit) # quit time consuming
  WAIC = waic(log_lik)
  looStat = loo(log_lik)
  save("WAIC", "looStat", file = sprintf("%s_waic.RData", fileName))
  fitSummary <- summary(fit,pars = c(paraNames, "LL_all"), use_cache = F)$summary
  write.table(matrix(fitSummary, nrow = length(paraNames) + 1), file = sprintf("%s_summary.txt", fileName),  sep = ",",
              col.names = F, row.names=FALSE)
  
  # detmerine converge
  converge = all(fitSummary[,"Rhat"] < 1.1) & all(fitSummary[, "n_eff"] >100)
  return(converge)
}

modelFittingCVdb = function(thisTrialData, fileName, paraNames, model, modelName,nPara, low, up){
  #
  load("wtwSettings.RData")
  # simulation parameters
  nChain = 4
  nIter = 5000
  
  # determine wIni
  # since the participants' initial strategies are unlikely optimal
  # we multiple the optimal opportunity cost by subOptimalRatio
  subOptimalRatio = 0.9 
  if(any(paraNames  == "gamma") || modelName == "BL" ){
    wIni = (5/6 + 0.93) / 2 * stepDuration / (1 - 0.9) * subOptimalRatio
  }else{
    wIni = (5/6 + 0.93) / 2 * stepDuration * subOptimalRatio
  }
  
  
  # prepare input
  timeWaited = thisTrialData$timeWaited
  scheduledWait = thisTrialData$scheduledWait
  trialEarnings = thisTrialData$trialEarnings
  timeWaited[trialEarnings > 0] = scheduledWait[trialEarnings > 0]
  cond = unique(thisTrialData$condition)
  tMax = ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  nTimeSteps = tMax / stepDuration
  Ts = round(ceiling(timeWaited / stepDuration) + 1)
  data_list <- list(wIni = wIni,
                    nTimeSteps = nTimeSteps,
                    nPara = nPara,
                    N = length(timeWaited),
                    trialEarnings = trialEarnings,
                    Ts = Ts,
                    iti = iti,
                    stepDuration = stepDuration,
                    low =low,
                    up = up)
  fit = sampling(object = model, data = data_list, cores = 1, chains = nChain,
                 iter = nIter) 
  write.table(get_elapsed_time(fit), file = sprintf("%s_time.txt", fileName), sep = ",",
              col.names = F, row.names=FALSE)
  fitSummary <- summary(fit,pars = c(paraNames, "LL_all"), use_cache = F)$summary
  write.table(matrix(fitSummary, nrow = length(paraNames) + 1), file = sprintf("%s_summary.txt", fileName),  sep = ",",
              col.names = F, row.names=FALSE)
  # detmerine converge
  converge = all(fitSummary[,"Rhat"] < 1.1) & all(fitSummary[, "n_eff"] >100)
  return(converge)
}