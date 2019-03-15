# select modelFun by modelName
getRepModelFun = function(modelName){
  if(modelName == "full_model"){
    repModelFun = full_model
  }else if(modelName == "reduce_one_QwaitIni"){
    repModelFun = reduce_one_QwaitIni
  }else if(modelName == "reduce_two_QwaitIni"){
    repModelFun = reduce_Qwait_Ini
  }else if(modelName == "reduce_one_phi"){
    repModelFun = reduce_one_phi
  }else if(modelName == "reduce_one_gamma"){
    repModelFun = reduce_one_gamma
  }else if(modelName == "cons_arbitrary"){
    repModelFun = cons_arbitrary
  }else if(modelName == "R_learning"){
    repModelFun = R_learning
  }else if(modelName == "R_learning2"){
    repModelFun = R_learning2
  }else if(modelName == "R_learning3"){
    repModelFun = R_learning3
  }
  else{
    return("wrong model name!")
  }
  return(repModelFun)
}
################ monte ######################
reduce_one_QwaitIni = function(para, cond, scheduledWait){
  # parse para
  phi = para[1]
  tau = para[2]
  gamma = para[3]
  
  # determine number of trials 
  nTrial = length(scheduledWait)
  wIni = mean(wInisExp)
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  Qwait = rep(wIni, nTimeStep) 
  Qquit = wIni * gamma ^(iti / stepDuration)
  
  # initialize varibles for recording
  vaWaits = matrix(NA, nTimeStep, nTrial);
  vaWaits[,1] = Qwait
  vaQuits = vector(length = nTrial);
  vaQuits[1] = Qquit
  rewardDelays = rep(0, nTrial)
  
  # initialize totalSecs 
  totalSecs = 0
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # sample rewardDelay
    rewardDelay = scheduledWait[tIdx]
    
    # loop over steps
    t = 1
    while(t <= nTimeStep){
      # determine action
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
      action = ifelse(runif(1) < waitRate, 'wait', 'quit')
      # next reward 
      rewardOccur = rewardDelay <= timeTicks[t + 1] && rewardDelay > timeTicks[t] 
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue, 0) 
      
      # dertime whether to continue
      trialGoOn= (action == 'wait' && !rewardOccur && t < nTimeStep)
      if(!trialGoOn){
        # if the trial stops, track trialEarnings, timeWaited and rewardDelays
        trialEarnings[tIdx] = ifelse(nextReward == tokenValue, tokenValue, 0);
        timeWaited[tIdx] = ifelse(getReward, rewardDelay, ifelse(action == "quit", timeTicks[t], timeTicks[t+1]))
        rewardDelays[tIdx] = rewardDelay
        sellTime[tIdx] = totalSecs + ifelse(getReward, rewardDelay, timeWaited[tIdx]) 
        break
      }else{
        # otherwise, continue
        t = t + 1
      }
    }# end of the trial 
    
    # update totalSecs 
    totalSecs = totalSecs + iti+ ifelse(getReward, rewardDelay, timeWaited[tIdx]) 
    
    # update Qwait and Qquit and go to the next trail if t < nTimeStep
    if(tIdx < nTrial){
      if(nextReward == 0){
        nextWaitRateHat =  1 / sum(1  + exp((Qquit - Qwait[1])* tau))
        trialReward = nextWaitRateHat * Qwait[1] * gamma ^(iti / stepDuration) +
          (1 - nextWaitRateHat) * Qquit * gamma ^(iti / stepDuration)
      }else{
        trialReward = nextReward
      }
      
      if(action == 'wait'){
        Qwait[1 : t] = (1 - phi) * Qwait[1 : t] + phi * trialReward * gamma ^ rev((0 : (t - 1 )))
      }else{
        Qquit =  (1 - phi) * Qquit + phi *  trialReward
        if(t > 1){
          Qwait[1 : (t - 1)] = (1 - phi) * Qwait[1 : (t - 1)] +
            phi * trialReward * gamma ^ rev((1 : (t - 1 )))
        }
      }
      # track vaWaits and vaQuits 
      vaWaits[,tIdx + 1] = Qwait
      vaQuits[tIdx + 1] = Qquit
      
    }# end of the update
    
    if(Qquit < 0 || sum(Qwait < 0) > 0){
      browser()
    }
  } # end of all trials 
  
  outputs = list( 
    "trialNum" = 1 : nTrial,
    "trialEarnings" = trialEarnings,
    "timeWaited" = timeWaited,
    "sellTime" = sellTime, # used in wtw analysis
    "scheduledWait" = rewardDelays,
    "vaWaits" = vaWaits,
    "vaQuits" = vaQuits
  )
  return(outputs)
} #end of the function


################ monte ######################
reduce_one_QwaitIni = function(para, cond, scheduledWait){
  # parse para
  phi = para[1]
  tau = para[2]
  gamma = para[3]
  
  # determine number of trials 
  nTrial = length(scheduledWait)
  wIni = ifelse(cond == "HP", wInisExp[1], wInisExp[2])
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  Qwait = rep(wIni, nTimeStep) 
  Qquit = wIni * gamma ^(iti / stepDuration)
  
  # initialize varibles for recording
  vaWaits = matrix(NA, nTimeStep, nTrial);
  vaWaits[,1] = Qwait
  vaQuits = vector(length = nTrial);
  vaQuits[1] = Qquit
  rewardDelays = rep(0, nTrial)
  
  # initialize totalSecs 
  totalSecs = 0
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # sample rewardDelay
    rewardDelay = scheduledWait[tIdx]
    
    # loop over steps
    t = 1
    while(t <= nTimeStep){
      # determine action
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
      action = ifelse(runif(1) < waitRate, 'wait', 'quit')
      # next reward 
      rewardOccur = rewardDelay <= timeTicks[t + 1] && rewardDelay > timeTicks[t] 
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue, 0) 
      
      # dertime whether to continue
      trialGoOn= (action == 'wait' && !rewardOccur && t < nTimeStep)
      if(!trialGoOn){
        # if the trial stops, track trialEarnings, timeWaited and rewardDelays
        trialEarnings[tIdx] = ifelse(nextReward == tokenValue, tokenValue, 0);
        timeWaited[tIdx] = ifelse(getReward, rewardDelay, ifelse(action == "quit", timeTicks[t], timeTicks[t+1]))
        rewardDelays[tIdx] = rewardDelay
        sellTime[tIdx] = totalSecs + ifelse(getReward, rewardDelay, timeWaited[tIdx]) 
        break
      }else{
        # otherwise, continue
        t = t + 1
      }
    }# end of the trial 
    
    # update totalSecs 
    totalSecs = totalSecs + iti+ ifelse(getReward, rewardDelay, timeWaited[tIdx]) 
    
    # update Qwait and Qquit and go to the next trail if t < nTimeStep
    if(tIdx < nTrial){
      if(nextReward == 0){
        nextWaitRateHat =  1 / sum(1  + exp((Qquit - Qwait[1])* tau))
        trialReward = nextWaitRateHat * Qwait[1] * gamma ^(iti / stepDuration) +
          (1 - nextWaitRateHat) * Qquit * gamma ^(iti / stepDuration)
      }else{
        trialReward = nextReward
      }
      
      if(action == 'wait'){
        Qwait[1 : t] = (1 - phi) * Qwait[1 : t] + phi * trialReward * gamma ^ rev((0 : (t - 1 )))
      }else{
        Qquit =  (1 - phi) * Qquit + phi *  trialReward
        if(t > 1){
          Qwait[1 : (t - 1)] = (1 - phi) * Qwait[1 : (t - 1)] +
            phi * trialReward * gamma ^ rev((1 : (t - 1 )))
        }
      }
      # track vaWaits and vaQuits 
      vaWaits[,tIdx + 1] = Qwait
      vaQuits[tIdx + 1] = Qquit
      
    }# end of the update
    
    if(Qquit < 0 || sum(Qwait < 0) > 0){
      browser()
    }
  } # end of all trials 
  
  outputs = list( 
    "trialNum" = 1 : nTrial,
    "trialEarnings" = trialEarnings,
    "timeWaited" = timeWaited,
    "sellTime" = sellTime, # used in wtw analysis
    "scheduledWait" = rewardDelays,
    "vaWaits" = vaWaits,
    "vaQuits" = vaQuits
  )
  return(outputs)
} #end of the function


################ monte ######################
full_model = function(para, cond, scheduledWait){
  # parse para
  phi = para[1]
  tau = para[2]
  gamma = para[3]
  QwaitIni = para[4]
  
  # determine number of trials 
  nTrial = length(scheduledWait)
  wIni = ifelse(cond == "HP", wInisTheory[1], wInisTheory[2])
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  Qwait = rep(QwaitIni, nTimeStep) 
  Qquit = QwaitIni * gamma ^(iti / stepDuration)
  
  # initialize varibles for recording
  vaWaits = matrix(NA, nTimeStep, nTrial);
  vaWaits[,1] = Qwait
  vaQuits = vector(length = nTrial);
  vaQuits[1] = Qquit
  rewardDelays = rep(0, nTrial)
  
  # initialize totalSecs 
  totalSecs = 0
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # sample rewardDelay
    rewardDelay = scheduledWait[tIdx]
    
    # loop over steps
    t = 1
    while(t <= nTimeStep){
      # determine action
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
      action = ifelse(runif(1) < waitRate, 'wait', 'quit')
      # next reward 
      rewardOccur = rewardDelay <= timeTicks[t + 1] && rewardDelay > timeTicks[t] 
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue, 0) 
      
      # dertime whether to continue
      trialGoOn= (action == 'wait' && !rewardOccur && t < nTimeStep)
      if(!trialGoOn){
        # if the trial stops, track trialEarnings, timeWaited and rewardDelays
        trialEarnings[tIdx] = ifelse(nextReward == tokenValue, tokenValue, 0);
        timeWaited[tIdx] = ifelse(getReward, rewardDelay, ifelse(action == "quit", timeTicks[t], timeTicks[t+1]))
        rewardDelays[tIdx] = rewardDelay
        sellTime[tIdx] = totalSecs + ifelse(getReward, rewardDelay, timeWaited[tIdx]) 
        break
      }else{
        # otherwise, continue
        t = t + 1
      }
    }# end of the trial 
    
    # update totalSecs 
    totalSecs = totalSecs + iti+ ifelse(getReward, rewardDelay, timeWaited[tIdx]) 
    
    # update Qwait and Qquit and go to the next trail if t < nTimeStep
    if(tIdx < nTrial){
      if(nextReward == 0){
        nextWaitRateHat =  1 / sum(1  + exp((Qquit - Qwait[1])* tau))
        trialReward = nextWaitRateHat * Qwait[1] * gamma ^(iti / stepDuration) +
          (1 - nextWaitRateHat) * Qquit * gamma ^(iti / stepDuration)
      }else{
        trialReward = nextReward
      }
      
      if(action == 'wait'){
        Qwait[1 : t] = (1 - phi) * Qwait[1 : t] + phi * trialReward * gamma ^ rev((0 : (t - 1 )))
      }else{
        Qquit =  (1 - phi) * Qquit + phi *  trialReward
        if(t > 1){
          Qwait[1 : (t - 1)] = (1 - phi) * Qwait[1 : (t - 1)] +
            phi * trialReward * gamma ^ rev((1 : (t - 1 )))
        }
      }
      # track vaWaits and vaQuits 
      vaWaits[,tIdx + 1] = Qwait
      vaQuits[tIdx + 1] = Qquit
      
    }# end of the update
    
    if(Qquit < 0 || sum(Qwait < 0) > 0){
      browser()
    }
  } # end of all trials 
  
  outputs = list( 
    "trialNum" = 1 : nTrial,
    "trialEarnings" = trialEarnings,
    "timeWaited" = timeWaited,
    "sellTime" = sellTime, # used in wtw analysis
    "scheduledWait" = rewardDelays,
    "vaWaits" = vaWaits,
    "vaQuits" = vaQuits
  )
  return(outputs)
} #end of the function


################ monte ######################
cons_theoretic = function(para, cond, scheduledWait){
  # parse para
  phi = para[1]
  tau = para[2]
  gamma = para[3]
  
  # determine number of trials 
  nTrial = length(scheduledWait)
  wIni = ifelse(cond == "HP", wInisTheory[1], wInisTheory[2])
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  Qwait = rep(wIni, nTimeStep) 
  Qquit = wIni * gamma ^(iti / stepDuration)
  
  # initialize varibles for recording
  vaWaits = matrix(NA, nTimeStep, nTrial);
  vaWaits[,1] = Qwait
  vaQuits = vector(length = nTrial);
  vaQuits[1] = Qquit
  rewardDelays = rep(0, nTrial)
  
  # initialize totalSecs 
  totalSecs = 0
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # sample rewardDelay
    rewardDelay = scheduledWait[tIdx]
    
    # loop over steps
    t = 1
    while(t <= nTimeStep){
      # determine action
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
      action = ifelse(runif(1) < waitRate, 'wait', 'quit')
      # next reward 
      rewardOccur = rewardDelay <= timeTicks[t + 1] && rewardDelay > timeTicks[t] 
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue, 0) 
      
      # dertime whether to continue
      trialGoOn= (action == 'wait' && !rewardOccur && t < nTimeStep)
      if(!trialGoOn){
        # if the trial stops, track trialEarnings, timeWaited and rewardDelays
        trialEarnings[tIdx] = ifelse(nextReward == tokenValue, tokenValue, 0);
        timeWaited[tIdx] = ifelse(getReward, rewardDelay, ifelse(action == "quit", timeTicks[t], timeTicks[t+1]))
        rewardDelays[tIdx] = rewardDelay
        sellTime[tIdx] = totalSecs + ifelse(getReward, rewardDelay, timeWaited[tIdx]) 
        break
      }else{
        # otherwise, continue
        t = t + 1
      }
    }# end of the trial 
    
    # update totalSecs 
    totalSecs = totalSecs + iti+ ifelse(getReward, rewardDelay, timeWaited[tIdx]) 
    
    # update Qwait and Qquit and go to the next trail if t < nTimeStep
    if(tIdx < nTrial){
      if(nextReward == 0){
        nextWaitRateHat =  1 / sum(1  + exp((Qquit - Qwait[1])* tau))
        trialReward = nextWaitRateHat * Qwait[1] * gamma ^(iti / stepDuration) +
          (1 - nextWaitRateHat) * Qquit * gamma ^(iti / stepDuration)
      }else{
        trialReward = nextReward
      }
      
      if(action == 'wait'){
        Qwait[1 : t] = (1 - phi) * Qwait[1 : t] + phi * trialReward * gamma ^ rev((0 : (t - 1 )))
      }else{
        Qquit =  (1 - phi) * Qquit + phi *  trialReward
        if(t > 1){
          Qwait[1 : (t - 1)] = (1 - phi) * Qwait[1 : (t - 1)] +
            phi * trialReward * gamma ^ rev((1 : (t - 1 )))
        }
      }
      # track vaWaits and vaQuits 
      vaWaits[,tIdx + 1] = Qwait
      vaQuits[tIdx + 1] = Qquit
      
    }# end of the update
    
    if(Qquit < 0 || sum(Qwait < 0) > 0){
      browser()
    }
  } # end of all trials 
  
  outputs = list( 
    "trialNum" = 1 : nTrial,
    "trialEarnings" = trialEarnings,
    "timeWaited" = timeWaited,
    "sellTime" = sellTime, # used in wtw analysis
    "scheduledWait" = rewardDelays,
    "vaWaits" = vaWaits,
    "vaQuits" = vaQuits
  )
  return(outputs)
} #end of the function

##### cons_theoretic
cons_arbitrary= function(para, cond, scheduledWait){
  # parse para
  phi = para[1]
  tau = para[2]
  gamma = para[3]
  
  # determine number of trials 
  nTrial = length(scheduledWait)
  wIni = ifelse(cond == "HP", wInisExp[1], wInisExp[2])
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  Qwait = rep(wInis[1], nTimeStep) 
  Qquit = wInis[2]
  
  # initialize varibles for recording
  vaWaits = matrix(NA, nTimeStep, nTrial);
  vaWaits[,1] = Qwait
  vaQuits = vector(length = nTrial);
  vaQuits[1] = Qquit
  rewardDelays = rep(0, nTrial)
  
  # initialize totalSecs 
  totalSecs = 0
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # sample rewardDelay
    rewardDelay = scheduledWait[tIdx]
    
    # loop over steps
    t = 1
    while(t <= nTimeStep){
      # determine action
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
      action = ifelse(runif(1) < waitRate, 'wait', 'quit')
      # next reward 
      rewardOccur = rewardDelay <= timeTicks[t + 1] && rewardDelay > timeTicks[t] 
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue, 0) 
      
      # dertime whether to continue
      trialGoOn= (action == 'wait' && !rewardOccur && t < nTimeStep)
      if(!trialGoOn){
        # if the trial stops, track trialEarnings, timeWaited and rewardDelays
        trialEarnings[tIdx] = ifelse(nextReward == tokenValue, tokenValue, 0);
        timeWaited[tIdx] = ifelse(getReward, rewardDelay, ifelse(action == "quit", timeTicks[t], timeTicks[t+1]))
        rewardDelays[tIdx] = rewardDelay
        sellTime[tIdx] = totalSecs + ifelse(getReward, rewardDelay, timeWaited[tIdx]) 
        break
      }else{
        # otherwise, continue
        t = t + 1
      }
    }# end of the trial 
    
    # update totalSecs 
    totalSecs = totalSecs + iti+ ifelse(getReward, rewardDelay, timeWaited[tIdx]) 
    
    # update Qwait and Qquit and go to the next trail if t < nTimeStep
    if(tIdx < nTrial){
      if(nextReward == 0){
        nextWaitRateHat =  1 / sum(1  + exp((Qquit - Qwait[1])* tau))
        trialReward = nextWaitRateHat * Qwait[1] * gamma ^(iti / stepDuration) +
          (1 - nextWaitRateHat) * Qquit * gamma ^(iti / stepDuration)
      }else{
        trialReward = nextReward
      }
      
      if(action == 'wait'){
        Qwait[1 : t] = (1 - phi) * Qwait[1 : t] + phi * trialReward * gamma ^ rev((0 : (t - 1 )))
      }else{
        Qquit =  (1 - phi) * Qquit + phi *  trialReward
        if(t > 1){
          Qwait[1 : (t - 1)] = (1 - phi) * Qwait[1 : (t - 1)] +
            phi * trialReward * gamma ^ rev((1 : (t - 1 )))
        }
      }
      # track vaWaits and vaQuits 
      vaWaits[,tIdx + 1] = Qwait
      vaQuits[tIdx + 1] = Qquit
      
    }# end of the update
    
    if(Qquit < 0 || sum(Qwait < 0) > 0){
      browser()
    }
  } # end of all trials 
  
  outputs = list( 
    "trialNum" = 1 : nTrial,
    "trialEarnings" = trialEarnings,
    "timeWaited" = timeWaited,
    "sellTime" = sellTime, # used in wtw analysis
    "scheduledWait" = rewardDelays,
    "vaWaits" = vaWaits,
    "vaQuits" = vaQuits
  )
  return(outputs)
} #end of the function



################ monte ######################
R_learning = function(para, cond, scheduledWait){
  # parse para
  phi1 = para[1]
  phi2 = para[2]
  tau = para[3]
  QwaitIni = 0.1
  rewardRateIni = 0.1
  
  # determine number of trials 
  nTrial = length(scheduledWait)
  wIni = ifelse(cond == "HP", wInisTheory[1], wInisTheory[2])
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  Qwait = rep(QwaitIni, nTimeStep) 
  rewardRate = rewardRateIni
  
  # initialize varibles for recording
  vaWaits = matrix(NA, nTimeStep, nTrial);
  vaWaits[,1] = Qwait
  vaRewardRates = vector(length = nTrial);
  vaRewardRates[1] = rewardRate
  rewardDelays = rep(0, nTrial)
  
  # initialize totalSecs 
  totalSecs = 0
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  
  # initialize rewardDelay, action, t
  rewardDelay = scheduledWait[1]
  junk = cumsum(Qwait)
  sumExpJunk = sum(exp(tau * junk))
  actionProbs = unlist(lapply(1 : nTimeStep, function(x) exp(tau * junk[x]) / sumExpJunk ))
  action = sample(1:nTimeStep, size=1, replace=TRUE, prob= actionProbs) # time step to quit waiting
  t = ifelse(action * stepDuration >rewardDelay, ceiling(rewardDelay / stepDuration), action)
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # determine timeWaited, sellTime
    getReward = action * stepDuration >= rewardDelay
    nextReward = ifelse(getReward, tokenValue, 0);
    trialEarnings[tIdx] = nextReward
    timeWaited[tIdx] = ifelse(getReward, rewardDelay, t * stepDuration)
    rewardDelays[tIdx] = rewardDelay
    sellTime[tIdx] = totalSecs + ifelse(getReward, rewardDelay, timeWaited[tIdx]) 
    
    # update Qwait and rewardRate 
    if(tIdx < nTrial){
      # determine nextRewardDelay, nextAction and nextT
      nextRewardDelay = scheduledWait[tIdx + 1]
      junk = cumsum(Qwait)
      sumExpJunk = sum(exp(tau * junk))
      actionProbs = unlist(lapply(1 : nTimeStep, function(x) exp(tau * junk[x]) / sumExpJunk ))
      nextAction = sample(1:nTimeStep, size=1, replace=TRUE, prob= actionProbs) # time step to quit waiting
      nextT = ifelse(nextAction * stepDuration > nextRewardDelay, ceiling(nextRewardDelay / stepDuration), nextAction)
      
      delta = nextReward  - rewardRate * t + sum(Qwait[1 : nextT]) - sum(Qwait[1:t])
      Qwait[1:t] = Qwait[1:t] + phi1 * delta / t # divide by t to calculate the relative contributions
      rewardRate = rewardRate + phi2 * delta / t
        
      # track vaWaits and vaQuits 
      vaWaits[,tIdx + 1] = Qwait
      vaRewardRates[tIdx + 1] = rewardRate
      
      # prepare to the next state
      rewardDelay = nextRewardDelay
      action = nextAction
      t= nextT
      
    }# end of the update
  } # end of all trials 
  
  outputs = list( 
    "trialNum" = 1 : nTrial,
    "trialEarnings" = trialEarnings,
    "timeWaited" = timeWaited,
    "sellTime" = sellTime, # used in wtw analysis
    "scheduledWait" = rewardDelays,
    "vaWaits" = vaWaits,
    "vaRewardRates" = vaRewardRates
  )
  return(outputs)
} #end of the function


R_learning2 = function(para, cond, scheduledWait){
  # parse para
  phi1 = para[1]
  phi2 = para[2]
  tau = para[3]
  QwaitIni = 0
  rewardRateIni = 0
  
  # determine number of trials 
  nTrial = length(scheduledWait)
  wIni = ifelse(cond == "HP", wInisTheory[1], wInisTheory[2])
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
  bufferWindow = nTimeStep
  
  # initialize action values
  Qwait = seq(0, QwaitIni, length.out = nTimeStep) 
  rewardRate = rewardRateIni
  
  # initialize varibles for recording
  vaWaits = matrix(NA, nTimeStep, nTrial);
  vaWaits[,1] = Qwait
  vaRewardRates = vector(length = nTrial);
  vaRewardRates[1] = rewardRate
  rewardDelays = rep(0, nTrial)
  
  # initialize totalSecs 
  totalSecs = 0
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    rewardDelay = scheduledWait[tIdx]
    actionProbs = unlist(lapply(1 : nTimeStep, function(x) exp(tau * Qwait[x]) / sum(exp(Qwait * tau))))
    action = sample(1:nTimeStep, size=1, replace=TRUE, prob= actionProbs) # time step to quit waiting
    t = ifelse(action * stepDuration >rewardDelay, ceiling(rewardDelay / stepDuration), action)
    
    # determine timeWaited, sellTime
    getReward = action * stepDuration >= rewardDelay
    nextReward = ifelse(getReward, tokenValue, 0);
    trialEarnings[tIdx] = nextReward
    timeWaited[tIdx] = ifelse(getReward, rewardDelay, t * stepDuration)
    rewardDelays[tIdx] = rewardDelay
    sellTime[tIdx] = totalSecs + ifelse(getReward, rewardDelay, timeWaited[tIdx]) 
    
    # update Qwait and rewardRate 
    if(tIdx < nTrial){
      if(nextReward > 0){
        delta = (nextReward + rewardRate * (nTimeStep - t)) / nTimeStep - rewardRate - Qwait[t : nTimeStep]
        Qwait[t : nTimeStep] = Qwait[t : nTimeStep] + phi1 * delta  
        if(t > 1){
          delta = rewardRate * (nTimeStep - 1 : (t-1)) / nTimeStep - rewardRate - Qwait[1 : (t -1)]
          Qwait[1 : (t - 1)] =   Qwait[1 : (t - 1)] + phi1 * delta
        }
      }else{
        if(t > 1){
          delta = rewardRate * (nTimeStep - 1 : (t-1)) / nTimeStep - rewardRate - Qwait[1 : (t -1)]
          Qwait[1 : (t - 1)] =   Qwait[1 : (t - 1)] + phi1 * delta
        }
        delta = rewardRate * (nTimeStep - t) / nTimeStep - rewardRate - Qwait[t : nTimeStep]
        Qwait[t : nTimeStep] = Qwait[t : nTimeStep] + phi1 * delta  
      }
      rewardRate = rewardRate + (nextReward + rewardRate * (nTimeStep - t)) / nTimeStep * phi2
      
      # track vaWaits and vaQuits 
      vaWaits[,tIdx + 1] = Qwait
      vaRewardRates[tIdx + 1] = rewardRate
      
      txt = sprintf("R = %d, T = %d", nextReward, t)
      plotData = data.frame(time = 1 : nTimeStep, Qwait = Qwait)
      ggplot(plotData, aes(time, Qwait)) + geom_point() + ggtitle(txt)
    }# end of the update
  } # end of all trials 
  
  outputs = list( 
    "trialNum" = 1 : nTrial,
    "trialEarnings" = trialEarnings,
    "timeWaited" = timeWaited,
    "sellTime" = sellTime, # used in wtw analysis
    "scheduledWait" = rewardDelays,
    "vaWaits" = vaWaits,
    "vaRewardRates" = vaRewardRates
  )
  return(outputs)
} #end of the function


mv = function(para, cond, scheduledWait){
  # parse para
  phi1 = para[1]
  phi2 = para[2]
  tau = para[3]
  
  # determine number of trials 
  nTrial = length(scheduledWait)
  wIni = ifelse(cond == "HP", wInisTheory[1], wInisTheory[2])
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
  bufferWindow = nTimeStep
  
  # initialize action values
  Pt = seq(0, 1, length.out = nTimeStep) 
  tauT = log(1 : nTimeStep)
  
  # initialize varibles for recording
  rewardDelays = rep(0, nTrial)
  
  # initialize totalSecs 
  totalSecs = 0
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    rewardDelay = scheduledWait[tIdx]
    actionProbs = unlist(lapply(1 : nTimeStep, function(x) exp(tau * Qwait[x]) / sum(exp(Qwait * tau))))
    action = sample(1:nTimeStep, size=1, replace=TRUE, prob= actionProbs) # time step to quit waiting
    t = ifelse(action * stepDuration >rewardDelay, ceiling(rewardDelay / stepDuration), action)
    
    # determine timeWaited, sellTime
    getReward = action * stepDuration >= rewardDelay
    nextReward = ifelse(getReward, tokenValue, 0);
    trialEarnings[tIdx] = nextReward
    timeWaited[tIdx] = ifelse(getReward, rewardDelay, t * stepDuration)
    rewardDelays[tIdx] = rewardDelay
    sellTime[tIdx] = totalSecs + ifelse(getReward, rewardDelay, timeWaited[tIdx]) 
    
    # update Qwait and rewardRate 
    if(tIdx < nTrial){
      if(nextReward > 0){
        delta = (nextReward + rewardRate * (nTimeStep - t)) / nTimeStep - rewardRate - Qwait[t : nTimeStep]
        Qwait[t : nTimeStep] = Qwait[t : nTimeStep] + phi1 * delta  
        if(t > 1){
          delta = rewardRate * (nTimeStep - 1 : (t-1)) / nTimeStep - rewardRate - Qwait[1 : (t -1)]
          Qwait[1 : (t - 1)] =   Qwait[1 : (t - 1)] + phi1 * delta
        }
      }else{
        if(t > 1){
          delta = rewardRate * (nTimeStep - 1 : (t-1)) / nTimeStep - rewardRate - Qwait[1 : (t -1)]
          Qwait[1 : (t - 1)] =   Qwait[1 : (t - 1)] + phi1 * delta
        }
        delta = rewardRate * (nTimeStep - t) / nTimeStep - rewardRate - Qwait[t : nTimeStep]
        Qwait[t : nTimeStep] = Qwait[t : nTimeStep] + phi1 * delta  
      }
      rewardRate = rewardRate + (nextReward + rewardRate * (nTimeStep - t)) / nTimeStep * phi2
      
      # track vaWaits and vaQuits 
      vaWaits[,tIdx + 1] = Qwait
      vaRewardRates[tIdx + 1] = rewardRate
      
      txt = sprintf("R = %d, T = %d", nextReward, t)
      plotData = data.frame(time = 1 : nTimeStep, Qwait = Qwait)
      ggplot(plotData, aes(time, Qwait)) + geom_point() + ggtitle(txt)
    }# end of the update
  } # end of all trials 
  
  outputs = list( 
    "trialNum" = 1 : nTrial,
    "trialEarnings" = trialEarnings,
    "timeWaited" = timeWaited,
    "sellTime" = sellTime, # used in wtw analysis
    "scheduledWait" = rewardDelays,
    "vaWaits" = vaWaits,
    "vaRewardRates" = vaRewardRates
  )
  return(outputs)
} #end of the function


################ monte ######################
R_learning3 = function(para, cond, scheduledWait){
  # parse para
  phi1 = para[1]
  phi2 = para[2]
  tau = para[3]
  
  # determine number of trials 
  nTrial = length(scheduledWait)
  wIni = 2
  rewardRateIni = 0.5
    
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  Qwait = rep(wIni, nTimeStep) 
  rewardRate = rewardRateIni
  Qquit = wIni - rewardRateIni * iti
  
  # initialize varibles for recording
  vaWaits = matrix(NA, nTimeStep, nTrial);
  vaWaits[,1] = Qwait
  vaQuits = vector(length = nTrial);
  vaQuits[1] = Qquit
  rewardDelays = rep(0, nTrial)
  
  # initialize totalSecs 
  totalSecs = 0
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # sample rewardDelay
    rewardDelay = scheduledWait[tIdx]
    
    # loop over steps
    t = 1
    while(t <= nTimeStep){
      # determine action
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
      action = ifelse(runif(1) < waitRate, 'wait', 'quit')
      # next reward 
      rewardOccur = rewardDelay <= timeTicks[t + 1] && rewardDelay > timeTicks[t] 
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue, 0) 
      
      # dertime whether to continue
      trialGoOn= (action == 'wait' && !rewardOccur && t < nTimeStep)
      if(!trialGoOn){
        # if the trial stops, track trialEarnings, timeWaited and rewardDelays
        trialEarnings[tIdx] = ifelse(nextReward == tokenValue, tokenValue, 0);
        timeWaited[tIdx] = ifelse(getReward, rewardDelay, ifelse(action == "quit", timeTicks[t], timeTicks[t+1]))
        rewardDelays[tIdx] = rewardDelay
        sellTime[tIdx] = totalSecs + ifelse(getReward, rewardDelay, timeWaited[tIdx]) 
        break
      }else{
        # otherwise, continue
        t = t + 1
      }
    }# end of the trial 
    
    # update totalSecs 
    totalSecs = totalSecs + iti+ ifelse(getReward, rewardDelay, timeWaited[tIdx]) 
    
    # update Qwait and Qquit and go to the next trail if t < nTimeStep
    if(tIdx < nTrial){
      if(nextReward == 0){
        nextWaitRateHat =  1 / sum(1  + exp((Qquit - Qwait[1])* tau))
        trialReward = nextWaitRateHat * Qwait[1] * gamma ^(iti / stepDuration) +
          (1 - nextWaitRateHat) * Qquit * gamma ^(iti / stepDuration)
      }else{
        trialReward = nextReward
      }
      
      if(action == 'wait'){
        Qwait[1 : t] = (1 - phi) * Qwait[1 : t] + phi * trialReward * gamma ^ rev((0 : (t - 1 )))
      }else{
        Qquit =  (1 - phi) * Qquit + phi *  trialReward
        if(t > 1){
          Qwait[1 : (t - 1)] = (1 - phi) * Qwait[1 : (t - 1)] +
            phi * trialReward * gamma ^ rev((1 : (t - 1 )))
        }
      }
      # track vaWaits and vaQuits 
      vaWaits[,tIdx + 1] = Qwait
      vaQuits[tIdx + 1] = Qquit
      
    }# end of the update
    
    if(Qquit < 0 || sum(Qwait < 0) > 0){
      browser()
    }
  } # end of all trials 
  
  outputs = list( 
    "trialNum" = 1 : nTrial,
    "trialEarnings" = trialEarnings,
    "timeWaited" = timeWaited,
    "sellTime" = sellTime, # used in wtw analysis
    "scheduledWait" = rewardDelays,
    "vaWaits" = vaWaits,
    "vaQuits" = vaQuits
  )
  return(outputs)
} #end of the function