
likRatioTest = function(modelName1, modelName2, df, group = "all"){
  paras1 = getParas(modelName1)
  expPara1 = loadExpPara(paras1,
                         sprintf("genData/expModelFitting/%s", modelName1))
  useID1 = getUseID(expPara1, paras1)
  paras2 = getParas(modelName2)
  expPara2 = loadExpPara(paras2,
                         sprintf("genData/expModelFitting/%s", modelName2))
  useID2 = getUseID(expPara2, paras2)
  if(group == "HP" || group == "LP"){
    useID = allIDs[allIDs %in% useID1 & allIDs %in% useID2 & hdrData$condition == group]
  }
  
  sumLogEvi1 = filter(expPara1, id %in% useID) %>% summarise(sum(LL_all))
  sumLogEvi2 = filter(expPara2, id %in% useID) %>% summarise(sum(LL_all))
  
  delta = -2 * as.double(sumLogEvi1 - sumLogEvi2)
  p = pchisq(delta, df * length(useID), lower.tail=FALSE)
  return(p)
}

# select modelFun by modelName
getLogLikFun = function(modelName){
  if(modelName == "para4"){
    logLikFun = para4
  }else if(modelName == "baseline"){
    logLikFun = baseLine
  }else if(modelName == "uniPrior"){
    logLikFun = uniPrior
  }else if(modelName == "uniPriorNC"){
    logLikFun = uniPriorNC
  }else if(modelName == "hyper"){
    logLikFun = hyper
  }else if(modelName == "PRbs"){
    logLikFun = PRbs
  }else if(modelName == "PRbsNC"){
    logLikFun = PRbsNC
  }else if(modelName == "MVT"){
    logLikFun = MVT
  }else if(modelName %in% c("Rlearn", "Rlearndb")){
    logLikFun = Rlearn
  }else if(modelName == "RlearnL"){
    logLikFun = RlearnL
  }else if(modelName == "reduce_gamma"){
    logLikFun = reduce_gamma
  }else{
    return("wrong model name!")
  }
  return(logLikFun)
}


PRbs = function(paras, cond, trialEarnings, timeWaited){
  # parse para
  phi = paras[1]; phiP = paras[2]; tau = paras[3]; gamma = paras[4]; zeroPoint = paras[5]
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  nTimeStep = tMax / stepDuration
  Ts = round(ceiling(timeWaited / stepDuration) + 1)
  
  # initialize actionValues
  subOptimalRatio = 0.9
  QHPApOptim = 5 / 6 * stepDuration / (1 - 0.9) * subOptimalRatio
  QLPApOptim = 0.93 * stepDuration / (1 - 0.9) * subOptimalRatio
  wIni = (QHPApOptim + QLPApOptim)/ 2 
  
  Qquit = wIni; Viti = wIni 
  Qwait = zeroPoint*0.1 - 0.1*(0 : (nTimeStep - 1)) + Qquit
  
  # initialize varibles for recording action values
  Qwaits = matrix(NA, nTimeStep, nTrial); Qwaits[,1] = Qwait
  Qquits = vector(length = nTrial); Qquits[1] = Qquit
  Vitis = vector(length = nTrial); Vitis[1] = Viti
  deltas = matrix(NA, nTimeStep, nTrial)
  Gs = matrix(NA, nTimeStep, nTrial)
  
  # initialize outputs 
  logLik_ = matrix(nrow = nTimeStep, ncol = nTrial)
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    nextReward = trialEarnings[tIdx]
    getReward = ifelse(nextReward == tokenValue, T, F)
    T = Ts[tIdx]
    # calculate logLik
    logLik_[,tIdx] =  sapply(1 : nTimeStep, function(i) 
      log(1 / sum(1  + exp((Qquit- Qwait[i])* tau))))
    # update values 
    if(tIdx < nTrial){
      returns = sapply(1 : (T-1), function(t) gamma^(T-t-1) *nextReward + gamma^(T-t) * Viti)
      if(getReward){
        Gs[1 : (T-1), tIdx] = returns[1 : (T-1)];
        deltas[1 : (T-1), tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
        Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
      }else{
        if(T > 2){
          Gs[1 : (T-2), tIdx] = returns[1 : (T-2)]
          deltas[1 : (T-2), tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phiP*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
        }
      }
      # update Viti
      if(getReward){
        Viti = Viti + phi*(gamma^(iti / stepDuration) * returns[1] - Viti)
      }else{
        Viti = Viti + phiP*(gamma^(iti / stepDuration) * returns[1] - Viti)
      }
      
      # update Qquit by counterfactual learning
      if(getReward){
        Qquit = Qquit + phi*(gamma^(iti / stepDuration + 1) * returns[1] - Qquit)
      }else{
        Qquit = Qquit + phiP*(gamma^(iti / stepDuration + 1) * returns[1] - Qquit)
      }
      
      # record updated values
      Qwaits[,tIdx + 1] = Qwait
      Qquits[tIdx + 1] = Qquit
      Vitis[tIdx + 1] = Viti
    }# end of the value update section
    
  } # end of the trial loop
  
  outputs = list( 
    "logLik_" = logLik_,
    "Qwaits" = Qwaits, "Qquits" = Qquits, "Gs" = Gs, "deltas" = deltas, "Vitis" = Vitis
  )
  return(outputs)
}

PRbsNC = function(paras, cond, trialEarnings, timeWaited){
  # parse para
  phi = paras[1]; phiP = paras[2]; tau = paras[3]; gamma = paras[4]; zeroPoint = paras[5]
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  nTimeStep = tMax / stepDuration
  Ts = round(ceiling(timeWaited / stepDuration) + 1)
  # initialize actionValues
  subOptimalRatio = 0.9
  QHPApOptim = 5 / 6 * stepDuration / (1 - 0.9) * subOptimalRatio
  QLPApOptim = 0.93 * stepDuration / (1 - 0.9) * subOptimalRatio
  wIni = (QHPApOptim + QLPApOptim)/ 2 
  
  Qquit = wIni; Viti = wIni 
  Qwait = zeroPoint*0.1 - 0.1*(0 : (nTimeStep - 1)) + Qquit
  
  # initialize varibles for recording action values
  Qwaits = matrix(NA, nTimeStep, nTrial); Qwaits[,1] = Qwait
  Qquits = vector(length = nTrial); Qquits[1] = Qquit
  Vitis = vector(length = nTrial); Vitis[1] = Viti
  deltas = matrix(NA, nTimeStep, nTrial)
  Gs = matrix(NA, nTimeStep, nTrial)
  
  # initialize outputs 
  logLik_ = matrix(nTimeStep, nTrial)
    
  # loop over trials
  for(tIdx in 1 : nTrial) {
    thisScheduledWait = scheduledWait[tIdx]
    # loop for each timestep t and determine At
    logLik_[,tIdx] =  log(1 / sum(1  + exp((rep(Qquit, nTimeStep) - Qwait)* tau)))
    
    # update values 
    if(tIdx < nTrial){
      nextReward = trialEarnings[tIdx]
      getReward = ifelse(nextReward == tokenValue, T, F)
      T = Ts[tIdx]
      returns = sapply(1 : (T-1), function(t) gamma^(T-t-1) *nextReward + gamma^(T-t) * Viti)
      if(getReward){
        Gs[1 : (T-1), tIdx] = returns[1 : (T-1)];
        deltas[1 : (T-1), tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
        Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
      }else{
        if(T > 2){
          Gs[1 : (T-2), tIdx] = returns[1 : (T-2)]
          deltas[1 : (T-2), tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phiP*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
        }
      }
      # update Viti
      if(getReward){
        Viti = Viti + phi*(gamma^(iti / stepDuration) * returns[1] - Viti)
      }else{
        Viti = Viti + phiP*(gamma^(iti / stepDuration) * returns[1] - Viti)
      }
      
      # update Qquit by counterfactual learning
      if(tIdx > 1){
        if(trialEarnings[tIdx - 1] == 0){
          if(getReward){
            Qquit = Qquit + phi*(gamma^(iti / stepDuration + 1) * returns[1] - Qquit)
          }else{
            Qquit = Qquit + phiP*(gamma^(iti / stepDuration + 1) * returns[1] - Qquit)
          }
        }
      }
      # record updated values
      Qwaits[,tIdx + 1] = Qwait
      Qquits[tIdx + 1] = Qquit
      Vitis[tIdx + 1] = Viti
    }# end of the value update section
    
  } # end of the trial loop
  
  outputs = list( 
    "logLik_" = logLik_,
    "Qwaits" = Qwaits, "Qquits" = Qquits, "Gs" = Gs, "deltas" = deltas, "Vitis" = Vitis
  )
  return(outputs)
}



