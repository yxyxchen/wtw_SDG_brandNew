# we don't need lp__ which is a sum loglikelyhood scaled by a constant, something like a dispersion 
# we don't need lp__ for model comparison 
# we save LL_all 
# we use log_like to calculate WAIC and looStat
# but we don't save log_like
modelFitting = function(cond, wIni, timeWaited, trialEarnings, scheduledWait, fileName, pars, model){
  load("wtwSettings.RData")
  tMax = ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  condIdx = ifelse(cond =="HP", 1, 2)
  nChain = 4
  nIter = 5000
  nTimeSteps = tMax / stepDuration
  Ts = round(ifelse(trialEarnings >0, ceiling(timeWaited / stepDuration) + 1, floor(timeWaited / stepDuration) + 1))
  # If using the least distance loss function
  # nScheduledWaitPoints = ceiling(scheduledWait / stepDuration)
  data_list <- list(tMax = tMax,
                    # nScheduledWaitPoints = nScheduledWaitPoints,
                    wIni = wIni,# maybe change later
                    nTimeSteps = nTimeSteps,
                    N = length(timeWaited),
                    timeWaited = timeWaited,
                    trialEarnings = trialEarnings,
                    Ts = Ts)
  fit = sampling(object = model, data = data_list, cores = min(nChain, 3), chains = nChain,
               iter = nIter) 
  # extract parameters
  extractedPara = fit %>%
    rstan::extract(permuted = F, pars = c("Qwait[1]"))
  # save sampling sequences
  tempt = extractedPara %>%
    adply(2, function(x) x) %>%  # change arrays into 2-d dataframe 
    select(-chains) 
  write.table(matrix(unlist(tempt), ncol = length(pars) + 1), file = sprintf("%s.txt", fileName), sep = ",",
              col.names = F, row.names=FALSE) 
  # calculate and save WAIC
  log_lik = extract_log_lik(fit) # quit time consuming
  WAIC = waic(log_lik)
  looStat = loo(log_lik)
  save("WAIC", "looStat", file = sprintf("%s_waic.RData", fileName))
  fitSumary <- summary(fit,pars = c(pars, "LL_all"), use_cache = F)$summary
  write.table(matrix(fitSumary, nrow = length(pars) + 1), file = sprintf("%s_summary.txt", fileName),  sep = ",",
            col.names = F, row.names=FALSE)
}

# I think this is used to fited the Nolearning models, using the expParaMedian
modelFittingNo = function(cond, wIni, timeWaited, trialEarnings, scheduledWait, fileName, pars, model, expParaMedian){
  load("wtwSettings.RData")
  tMax = ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  condIdx = ifelse(cond =="HP", 1, 2)
  nChain = 4
  nIter = 5000
  nTimeStep = tMax / stepDuration
  nTimePoints = round(ifelse(trialEarnings >0, ceiling(timeWaited / stepDuration), floor(timeWaited / stepDuration) + 1))
  nScheduledWaitPoints = ceiling(scheduledWait / stepDuration)
  data_list <- list(tMax = tMax,
                    nScheduledWaitPoints = nScheduledWaitPoints,
                    wIni = wIni,
                    wInis = wInis,
                    nTimeStep = nTimeStep,
                    N = length(timeWaited),
                    timeWaited = timeWaited,
                    trialEarnings = trialEarnings,
                    nTimePoints = nTimePoints,
                    expParaMedian = expParaMedian)
  # init = list(list('phi' = 0.5, 'tau' = 15, 'gamma' = 0.5),
  #             list('phi' = 0.1, 'tau' = 5, 'gamma' = 0.1))
  fit = sampling(object = model, data = data_list, cores = min(nChain, 3), chains = nChain,
                 iter = nIter) 
  # extract parameters
  extractedPara = fit %>%
    rstan::extract(permuted = F, pars = c(pars, "LL_all"))
  # save sampling sequences
  tempt = extractedPara %>%
    adply(2, function(x) x) %>%  # change arrays into 2-d dataframe 
    select(-chains) 
  write.table(matrix(unlist(tempt), ncol = length(pars) + 1), file = sprintf("%s.txt", fileName), sep = ",",
              col.names = F, row.names=FALSE)
  # calculate and save WAIC
  log_lik = extract_log_lik(fit) # quit time consuming
  WAIC = waic(log_lik)
  looStat = loo(log_lik)
  save("WAIC", "looStat", file = sprintf("%s_waic.RData", fileName))
  # save dist and timeWaitedProbLog, not the same as simulation but ...
  # junk =  fit %>%
  #   rstan::extract(permuted = F, pars =  "timeWaitedLogProb")
  # timeWaitedProb = as.vector(apply(junk, MARGIN = 3, function(x) mean(exp(x))))
  # timeWaitedProbSd = as.vector(apply(junk, MARGIN = 3, function(x) sd(exp(x))))
  # junk = fit %>%
  #   rstan::extract(permuted = F, pars =  "dist") %>% 
  #   adply(2, function(x) x) %>%  select(-chains)
  # dist = as.vector(apply(junk, 2, mean))
  # distSd = as.vector(apply(junk, 2, sd))
  # save("timeWaitedProb", "timeWaitedProbSd", "dist", "distSd", file = sprintf("%s_prediction.RData", fileName))
  # save summarized fit 
  fitSumary <- summary(fit,pars = c(pars, "lp__", "LL_all"), use_cache = F)$summary
  write.table(matrix(fitSumary, nrow = length(pars) + 2), file = sprintf("%s_summary.txt", fileName),  sep = ",",
              col.names = F, row.names=FALSE)
}