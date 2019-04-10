# select modelFun by modelName
getRepModelFun = function(modelName){
  if(modelName == "full_model"){
    repModelFun = full_model
  }else if(modelName == "monteRatio"){
    repModelFun = monteRatio
  }else if(modelName == "curiosityTrial"){
    repModelFun = curiosityTrial
  }else{
    return("wrong model name!")
  }
  return(repModelFun)
}

################ monte ######################
curiosityTrial = function(paras, cond, scheduledWait){
  # parse para
  phi = paras[1]
  tau = paras[2]
  gamma = paras[3]
  
  # coefficient of curiosity
  curSlope = 0.2
  curIntercept = 2
  
  # determine number of trials 
  nTrial = length(scheduledWait)
  wIni = mean(wInisExp)
  QHPApOptim = 3.937851
  QLPApOptim = 4.396877
  wIni = (QHPApOptim + QLPApOptim)/ 2
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  Qwait = rep(wIni, nTimeStep) 
  Qquit = wIni 
  
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
  
  curiosity = curIntercept
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # sample rewardDelay
    rewardDelay = scheduledWait[tIdx]
    
    # loop over steps
    t = 1
    while(t <= nTimeStep){
      # determine action
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t] - curiosity)* tau))
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
    
    # update curiosity
    nextCuriosity =  curIntercept * exp(-curSlope*tIdx)
    # update Qwait and Qquit and go to the next trail if t < nTimeStep
    if(tIdx < nTrial){
      # determine the update target by one-step bellman backup
      # in determing nextQ, we use the unupdated policy, yet using the new curiosity
      nextWaitRateS1 =  1 / sum(1  + exp((Qquit - Qwait[1] - nextCuriosity)* tau))
      nextQ = nextWaitRateS1 * Qwait[1] +
        (1 - nextWaitRateS1) * Qquit 
      trialReward = nextReward + nextQ * gamma ^(iti / stepDuration)
      
      # the learning rule is very naive,
      # first, it only update Qquit when the agent quit
      # more importantly, the update target of Qwait[1] should not be trialReward * gamma ^ rev((1 : t))
      if(action == 'wait'){     
        Qwait[1 : t] = (1 - phi) * Qwait[1 : t] + phi * trialReward * gamma ^ rev((1 : t))
        # counterfactual thinking
        Qquit = (1 - phi) * Qquit + phi * trialReward * gamma ^ ((iti / stepDuration) + t)
      }else{
        Qquit =  (1 - phi) * Qquit + phi *  trialReward
        if(t > 1){
          Qwait[1 : (t - 1)] = (1 - phi) * Qwait[1 : (t - 1)] +
            phi * trialReward * gamma ^ rev((1 : (t - 1 )))
        }
        # counterfactual thinking
        # Qquit = (1 - phi) * Qquit + phi * trialReward * gamma ^ ((iti / stepDuration) + t)
      }
      # track vaWaits and vaQuits 
      vaWaits[,tIdx + 1] = Qwait
      vaQuits[tIdx + 1] = Qquit
      
    }# end of the update
    
    # update curiosity
    curiosity = nextCuriosity
    
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
bellmanStep2 = function(paras, cond, scheduledWait){
  # parse para
  phi = paras[1]
  tau = paras[2]
  gamma = paras[3]
  
  # coefficient of curiosity
  curSlope = 0.2
  curIntercept = 2
  
  # determine number of trials 
  nTrial = length(scheduledWait)
  wIni = mean(wInisExp)
  QHPApOptim = 3.937851
  QLPApOptim = 4.396877
  wIni = (QHPApOptim + QLPApOptim)/ 2
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  Qwait = rep(wIni, nTimeStep) 
  Qquit = wIni 
  
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
  
  curiosity = curIntercept
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # sample rewardDelay
    rewardDelay = scheduledWait[tIdx]
    
    # loop over steps
    t = 1
    while(t <= nTimeStep){
      # determine action
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t] - curiosity)* tau))
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
    
    # update curiosity
    curiosity =  curIntercept * exp(-curSlope*(tIdx - 1))
    # update Qwait and Qquit and go to the next trail if t < nTimeStep
    if(tIdx < nTrial) {
      # determine the update target by one-step bellman backup
      # in determing nextQ, we use the unupdated policy, yet using the new curiosity
      nextWaitRateS1 =  1 / sum(1  + exp((Qquit - Qwait[1] - curiosity)* tau))
      nextQ = nextWaitRateS1 * Qwait[1] +
        (1 - nextWaitRateS1) * Qquit 
      trialReward = nextReward + nextQ * gamma ^(iti / stepDuration)
      
      # the learning rule is very naive,
      # first, it only update Qquit when the agent quit
      # more importantly, the update target of Qwait[1] should not be trialReward * gamma ^ rev((1 : t))
      if(action == 'wait'){     
        Qwait[1 : t] = (1 - phi) * Qwait[1 : t] + phi * trialReward * gamma ^ rev((1 : t))
        # counterfactual thinking
        Qquit = (1 - phi) * Qquit + phi * trialReward * gamma ^ ((iti / stepDuration) + t)
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
