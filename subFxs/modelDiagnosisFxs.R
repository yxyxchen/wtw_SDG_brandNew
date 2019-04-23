getLogLikFun = function(modelName){
  if(modelName == "curiosityTrialR"){
    logLikFun = LL_curiosityTrialR
  }else if(modelName == "curiosityTrial"){
    logLikFun = LL_curiosityTrial
  }else{
    return("wrong model name!")
  }
  return(logLikFun)
}

LL_curiosityTrialR = function(paras, blockData){
  # parse paras
  phi = paras[1]
  tau = paras[2]
  phiR = paras[3]
  #phiR = 0.005
  
  # parse blockData
  timeWaited = blockData$timeWaited
  trialEarnings = blockData$trialEarnings
  cond = unique(blockData$condition)
  
  # coefficient of curiosity
  curSlope = 0.2
  curIntercept = 2
  
  # determine number of trials and nTimeSteps 
  nTrial = length(timeWaited)
  Ts =  round(ifelse(trialEarnings >0, ceiling(timeWaited / stepDuration) + 1,
                     floor(timeWaited / stepDuration) + 1))
  tMax = ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  nTimeStep = tMax / stepDuration
  
  # initialize actionValues
  # here we use the optimal reward rates from the normative analysis in Lempert 2018
  # it is more accurate then the one I calcualte in wtwSettings.R
  wIni = 0
  Qwait = rep(0, nTimeStep) 
  Qquit = 0
  Viti = 0
  Rrate = wIni
  
  # initialize varibles for recording action values
  Qwaits = matrix(NA, nTimeStep, nTrial);
  Qwaits[,1] = Qwait
  Qquits = vector(length = nTrial);
  Qquits[1] = Qquit
  Vitis = vector(length = nTrial);
  Vitis[1] = Viti
  Rrates = vector(length = nTrial);
  Rrates[1] = Rrate
  
  #initialize LL
  logLik = vector(length = sum(Ts) - nTrial)
  LL = 0 # sum of logLik
  stepNo = 0;
  
  for(tIdx in 1 : nTrial){
    # determine curiosity, terminal state T 
    curiosity = curIntercept * exp(-curSlope*(tIdx-1))
    T = Ts[tIdx]
    thisTrialEarn = trialEarnings[tIdx]
    thisTimeWaited = timeWaited[tIdx]
    
    # update LL
    if(thisTrialEarn == 0) {
      thisTrialLL = c(sapply(1 : (T-2), function(t) log(1 / (1 + exp(tau * (Qquit - Qwait[t] - curiosity)))) ),
                      log(1 - 1 / (1 + exp(tau * (Qquit - Qwait[T-1] - curiosity)))))
    }else{
      thisTrialLL = sapply(1 : (T-1), function(t) log(1 / (1 + exp(tau * (Qquit - Qwait[t] - curiosity)))) )  
    }
    logLik[(stepNo + 1) : (stepNo + T - 1)] = thisTrialLL
    LL = LL + sum(thisTrialLL)
    
    # update actionValues
    if(thisTrialEarn > 0){
      G = sapply(1 : (T - 1), function(t) thisTrialEarn - Rrate * (T - t) + Viti) # update target
      Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi * (G - Qwait[1 : (T-1)])
    }else{
      G = thisTrialEarn - Rrate  + Viti
      Qquit = Qquit + phi * (G - Qquit)
      if(T > 2){
        G = sapply(1 : (T - 2), function(t) thisTrialEarn - Rrate * (T - t) + Viti) # update target
        Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phi * (G - Qwait[1 : (T-2)])        
      }
    }
    
    # update conterfactual thinking
    G1 =  thisTrialEarn - Rrate * (T - 1) + Viti;
    Qquit = Qquit + phi * (G1 - Rrate * (iti /stepDuration + 1) - Qquit);
    
    # update Viti and Rrate
    delta = (G1 - Rrate * (iti /stepDuration) - Viti);
    Viti = Viti + phi * delta;
    Rrate = Rrate + phiR * delta;
    
    # record
    if(tIdx < (nTrial)){
      Qwaits[,tIdx + 1] = Qwait
      Qquits[tIdx + 1] = Qquit
      Rrates[tIdx + 1] = Rrate
      Vitis[tIdx + 1] = Viti
    }
  }
  outputs = list("Qwaits" = Qwaits,
                 "Qquits" = Qquits,
                 "Rrates" = Rrates,
                 "Vitis" = Vitis,
                 "LL" = LL,
                 "logLik" = logLik)
  return(outputs)
}

outputs = LL_curiosityTrialR(paras, thisTrialData)
plot(outputs$Rrates)


plot(outputs$Qwaits[,41])
plot(outputs$Qwaits[,20])

