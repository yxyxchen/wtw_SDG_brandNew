# select modelFun by modelName
getRepModelFun = function(modelName){
  if(modelName == "full_model"){
    repModelFun = full_model
  }else if(modelName == "curiosityTrialR"){
    repModelFun = curiosityTrialR
  }else if(modelName == "curiosityTrial"){
    repModelFun = curiosityTrial
  }else if(modelName == "heuristicRL"){
    repModelFun = heuristicRL
  }else if(modelName == "heuristicRLAve"){
    repModelFun = heuristicRLAve
  }else{
    return("wrong model name!")
  }
  return(repModelFun)
}

################ curiosityTrial model using the Rlearning  ######################
curiosityTrialR = function(paras, cond, scheduledWait){
  # parse para
  phi = paras[1]
  tau = paras[2]
  phiR = paras[3]
  
  # coefficient of curiosity
  curSlope = 0.2
  curIntercept = 2
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  nTimeStep = tMax / stepDuration
  
  # initialize actionValues
  # here we use the optimal reward rates from the normative analysis in Lempert 2018
  # it is more accurate then the one I calcualte in wtwSettings.R
  wIni = (5/6 + 0.93)/ 2 * stepDuration
  
  Qwait = c(rep(1, nTimeStep))
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
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  
  # initialize elapsed time
  elapsedTime = 0
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # determine 
    thisScheduledWait = scheduledWait[tIdx]
    #curiosity =  curIntercept * exp(-curSlope*(tIdx-1))
    curiosity =  0
    
    # loop for each timestep t and determine At
    t = 1
    while(t <= nTimeStep){
      # determine At
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t] - curiosity)* tau))
      action = ifelse(runif(1) < waitRate, 'wait', 'quit')
      # observe St+1 and Rt+1
      rewardOccur = thisScheduledWait <= (t * stepDuration) && thisScheduledWait > ((t-1) * stepDuration)
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue, 0) 
      
      # dertime whether St+1 is the terminal state
      # if the trial terminates, track terminal timestep index T, trialEarnings, timeWaited, sellTime and elapsedTime
      # otherwise, continue
      nextStateTerminal = (getReward || action == "quit")
      if(nextStateTerminal){
        T = t+1
        trialEarnings[tIdx] = ifelse(nextReward == tokenValue, tokenValue, 0);
        timeWaited[tIdx] = ifelse(getReward, thisScheduledWait, t * stepDuration)
        sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
        elapsedTime = elapsedTime + timeWaited[tIdx] + iti
        break
      }else{
        t = t + 1
      }
    }# end of the action selection section
    
    # update values 
    if(tIdx < nTrial){
      
      # update action values for each timestep t
      returns = sapply(1 : (T-1), function(t) nextReward - (T-t) * Rrate + Viti)
      # when the agent always wait and get the reward, update Qwait[1:(T-1)]
      # otherwise, update Qquit and Qwait[1 : (T-2)]      
      if(getReward){
        Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])        
      }else{
        Qquit = Qquit + phi*(returns[T-1] - Qquit)
        if(T > 2) Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phi*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
      }

      # update Qquit by counterfactual learning
      Qquit = Qquit + phi*(returns[1] - Rrate * (iti / stepDuration + 1) - Qquit)
      
      # update Viti and Rrate
      deltaIti = returns[1] - Rrate * (iti / stepDuration) - Viti
      Viti = Viti + phi*deltaIti
      Rrate = Rrate + phiR *deltaIti
      
      # record updated values
      Qwaits[,tIdx + 1] = Qwait
      Qquits[tIdx + 1] = Qquit
      Vitis[tIdx + 1] = Viti
      Rrates[tIdx + 1] = Rrate
    }# end of the value update section
    
  } # end of the trial loop
  
  outputs = list( 
    "trialNum" = 1 : nTrial,
    "trialEarnings" = trialEarnings,
    "timeWaited" = timeWaited,
    "sellTime" = sellTime, # used in wtw analysis
    "scheduledWait" = scheduledWait,
    "Qwaits" = Qwaits,
    "Qquits" = Qquits,
    "Vitis" = Vitis,
    "Rrates" = Rrates
  )
  return(outputs)
} #end of the function

################ monte ######################
curiosityTrial = function(paras, cond, scheduledWait){
  # parse para
  phi = paras[1]
  tau = paras[2]
  gamma = paras[3]
  
  # coefficient of curiosity
  curSlope = 0.2
  curIntercept = 2
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  nTimeStep = tMax / stepDuration
  
  # initialize actionValues
  # here we use the optimal reward rates from the normative analysis in Lempert 2018
  # it is more accurate then the one I calcualte in wtwSettings.R
  # in addition, I use the gamma from 0.5s stepDuration, just hope the Q is similiar to the asympototic value in this RL
  # finally, we use / (1 - gamma) instead of the gamma / (1 - gamma), it assumes the results always happen as the begging 
  # so it is a upper
  # here we use 0.9 as the discount rate for one stepDuration
  QHPApOptim = 5 / 6 * stepDuration / (1 - 0.9) 
  QLPApOptim = 0.93 * stepDuration / (1 - 0.9) 
  wIni = (QHPApOptim + QLPApOptim)/ 2
  
  Qwait = c(rep(wIni, 80), rep(0, nTimeStep - 80))
  Qquit = wIni * 0.8
  Viti = wIni * 0.8
  
  # initialize varibles for recording action values
  Qwaits = matrix(NA, nTimeStep, nTrial);
  Qwaits[,1] = Qwait
  Qquits = vector(length = nTrial);
  Qquits[1] = Qquit
  Vitis = vector(length = nTrial);
  Vitis[1] = Viti
  deltas = matrix(NA, nTimeStep, nTrial)
  Gs = matrix(NA, nTimeStep, nTrial)
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  
  # initialize elapsed time
  elapsedTime = 0
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # determine 
    thisScheduledWait = scheduledWait[tIdx]
    # curiosity =  curIntercept * exp(-curSlope*(tIdx-1))
    curiosity = 0
    # loop for each timestep t and determine At
    t = 1
    while(t <= nTimeStep){
      # determine At
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t] - curiosity)* tau))
      action = ifelse(runif(1) < waitRate, 'wait', 'quit')
      # observe St+1 and Rt+1
      rewardOccur = thisScheduledWait <= (t * stepDuration) && thisScheduledWait > ((t-1) * stepDuration)
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue, 0) 

      # dertime whether St+1 is the terminal state
      # if the trial terminates, track terminal timestep index T, trialEarnings, timeWaited, sellTime and elapsedTime
      # otherwise, continue
      nextStateTerminal = (getReward || action == "quit")
      if(nextStateTerminal){
        T = t+1
        trialEarnings[tIdx] = ifelse(nextReward == tokenValue, tokenValue, 0);
        timeWaited[tIdx] = ifelse(getReward, thisScheduledWait, t * stepDuration)
        sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
        elapsedTime = elapsedTime + timeWaited[tIdx] + iti
        break
      }else{
        t = t + 1
      }
    }# end of the action selection section
    
    # update values 
    if(tIdx < nTrial){
      
      # update action values for each timestep t
      returns = sapply(1 : (T-1), function(t) gamma^(T-t-1) *nextReward + gamma^(T-t) * Viti)
      # returns = sapply(1 : (T-1), function(t) gamma^(T-t-1) *nextReward + gamma^(T-t) * Viti)
      # when the agent always wait and get the reward, update Qwait[1:(T-1)]
      # otherwise, update Qquit and Qwait[1 : (T-2)]      
      if(getReward){
        Gs[1 : (T-1), tIdx] = returns[1 : (T-1)]
        deltas[1 : (T-1), tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
        Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
      }else{
        Qquit = Qquit + phi*(returns[T-1] - Qquit)
        if(T > 2){
          Gs[1 : (T-2), tIdx] = returns[1 : (T-2)]
          deltas[1 : (T-2), tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phi*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
        }
      }
      
      # update Viti
      Viti = Viti + phi*(gamma^(iti / stepDuration) * returns[1] - Viti)
      
      # update Qquit by counterfactual learning
      Qquit = Qquit + phi*(gamma^(iti / stepDuration + 1) * returns[1] - Qquit)
      
      # record updated values
      Qwaits[,tIdx + 1] = Qwait
      Qquits[tIdx + 1] = Qquit
      Vitis[tIdx + 1] = Viti
    }# end of the value update section
    
  } # end of the trial loop
  
  outputs = list( 
    "trialNum" = 1 : nTrial,
    "trialEarnings" = trialEarnings,
    "timeWaited" = timeWaited,
    "sellTime" = sellTime, # used in wtw analysis
    "scheduledWait" = scheduledWait,
    "Qwaits" = Qwaits,
    "Qquits" = Qquits,
    "Gs" = Gs,
    "deltas" = deltas,
    "Vitis" = Vitis
  )
  return(outputs)
}




heuristicRL =  function(paras, cond, scheduledWait){
  phi = paras[1]
  threshd = paras[2]
  ini = paras[3]
  
  nTrial = length(scheduledWait)
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  expectedDurations = rep(0, nTrial)
  
  # initialize elapsed time
  elapsedTime = 0
  
  # initialize expectedDuration
  expectedDuration = ini
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # determine thisScheduledWait
    thisScheduledWait = scheduledWait[tIdx]
    # determine timeWaited, trialEarnings, sellTime and elapsedTime 
    if(thisScheduledWait <= (expectedDuration * threshd)){
      trialEarnings[tIdx] = tokenValue
      timeWaited[tIdx] = thisScheduledWait
      sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
      elapsedTime = elapsedTime + timeWaited[tIdx] + iti
    }else{
      trialEarnings[tIdx] = 0
      timeWaited[tIdx] = expectedDuration * threshd
      sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
      elapsedTime = elapsedTime + timeWaited[tIdx] + iti
    }
    # record expectedDurations
    expectedDurations[tIdx] = expectedDuration
    # update threshold
    expectedDuration = expectedDuration + phi*(timeWaited[tIdx] - expectedDuration)
  }  
  
  outputs = list( 
    "trialNum" = 1 : nTrial,
    "trialEarnings" = trialEarnings,
    "timeWaited" = timeWaited,
    "sellTime" = sellTime, # used in wtw analysis
    "scheduledWait" = scheduledWait,
    "expectedDurations" = expectedDuration
  )
  return(outputs)
}


heuristicRLAve =  function(paras, cond, scheduledWait){
  threshd = paras[1]
  ini = paras[2]
  
  nTrial = length(scheduledWait)
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  expectedDurations = rep(0, nTrial)
  
  # initialize elapsed time
  elapsedTime = 0
  
  # initialize expectedDuration
  expectedDuration = ini
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # determine thisScheduledWait
    thisScheduledWait = scheduledWait[tIdx]
    # determine timeWaited, trialEarnings, sellTime and elapsedTime 
    if(thisScheduledWait <= (expectedDuration * threshd)){
      trialEarnings[tIdx] = tokenValue
      timeWaited[tIdx] = thisScheduledWait
      sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
      elapsedTime = elapsedTime + timeWaited[tIdx] + iti
    }else{
      trialEarnings[tIdx] = 0
      timeWaited[tIdx] = expectedDuration * threshd
      sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
      elapsedTime = elapsedTime + timeWaited[tIdx] + iti
    }
    # record expectedDurations
    expectedDurations[tIdx] = expectedDuration
    # update threshold
    expectedDuration = expectedDuration + 1 / (tIdx + 1) *(timeWaited[tIdx] - expectedDuration)
  }  
  
  outputs = list( 
    "trialNum" = 1 : nTrial,
    "trialEarnings" = trialEarnings,
    "timeWaited" = timeWaited,
    "sellTime" = sellTime, # used in wtw analysis
    "scheduledWait" = scheduledWait,
    "expectedDurations" = expectedDuration
  )
  return(outputs)
}
