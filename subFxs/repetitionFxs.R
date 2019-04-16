# select modelFun by modelName
getRepModelFun = function(modelName){
  if(modelName == "full_model"){
    repModelFun = full_model
  }else if(modelName == "curiosityTrialR"){
    repModelFun = curiosityTrialR
  }else if(modelName == "curiosityTrial"){
    repModelFun = curiosityTrial
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
  phiR = paras[3] # learning rate for Rrate
  
  # coefficient of curiosity
  curSlope = 0.2
  curIntercept = 2
  
  # determine number of trials 
  nTrial = length(scheduledWait)
  wIni = (5/6 + 0.93) / 2 * stepDuration# initial value for Rrate
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  Rrate = wIni
  Qwait = rep(0, nTimeStep) 
  Qquit = 0 
  
  # initialize varibles for recording
  vaWaits = matrix(NA, nTimeStep, nTrial);
  vaWaits[,1] = Qwait
  vaQuits = vector(length = nTrial);
  vaQuits[1] = Qquit
  vaRrates =  vector(length = nTrial);
  vaRrates[1] = Rrate
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
      # no change in this part
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
      trialReward = nextReward + nextQ - (iti / stepDuration) * Rrate
      
      # the learning rule is very naive,
      # first, it only update Qquit when the agent quit
      # more importantly, the update target of Qwait[1] should not be trialReward * gamma ^ rev((1 : t))
      if(action == 'wait'){     
        Qwait[1 : t] = (1 - phi) * Qwait[1 : t] + phi * (trialReward - Rrate * rev(1:t))
        Rrate = Rrate * (1 - phiR) + phiR * (trialReward - Rrate * (t + iti/stepDuration))
        # counterfactual thinking
        Qquit = (1 - phi) * Qquit + phi * (trialReward - Rrate * (t+ iti / stepDuration))

      }else{
        # counterfactual thinking here
        Qquit =  (1 - phi) * Qquit + phi *  trialReward
        Rrate = Rrate * (1 - phiR) + phiR * (trialReward - Rrate * (t-1 + iti/stepDuration))
        # counterfactual thinking
        Qquit = (1 - phi) * Qquit + phi * (trialReward - Rrate * (t - 1 + iti / stepDuration))
        if(t > 1){
          Qwait[1 : (t - 1)] = (1 - phi) * Qwait[1 : (t - 1)] +
            phi * (trialReward - rev(1 : (t - 1)) * Rrate)
        }
        # counterfactual thinking
        # Qquit = (1 - phi) * Qquit + phi * trialReward * gamma ^ ((iti / stepDuration) + t)
      }
      # track vaWaits and vaQuits 
      vaWaits[,tIdx + 1] = Qwait
      vaQuits[tIdx + 1] = Qquit
      vaRrates[tIdx + 1] = Rrate
      
    }# end of the update
    
    # update curiosity
    curiosity = nextCuriosity
  } # end of all trials 
  
  outputs = list( 
    "trialNum" = 1 : nTrial,
    "trialEarnings" = trialEarnings,
    "timeWaited" = timeWaited,
    "sellTime" = sellTime, # used in wtw analysis
    "scheduledWait" = rewardDelays,
    "vaWaits" = vaWaits,
    "vaQuits" = vaQuits,
    "vaRrates" = vaRrates
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
  
  # determine number of trials 
  nTrial = length(scheduledWait)
  wIni = mean(wInisExp)
  # here we use the optimal reward rates from the normative analysis in Lempert 2018
  # it is more accurate then the one I calcualte in wtwSettings.R
  # in addition, I use the gamma from 0.5s stepDuration, just hope the Q is similiar to the asympototic value in this RL
  # finally, we use / (1 - gamma) instead of the gamma / (1 - gamma), it assumes the results always happen as the begging 
  # so it is a upper
  QHPApOptim = 5 / 6 * stepDuration / (1 - 0.9) 
  QLPApOptim = 0.93 * stepDuration / (1 - 0.9) 
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
        Qwait[1 : t] = (1 - phi) * Qwait[1 : t] + phi * trialReward * gamma ^ rev((0 : (t-1)))
        # counterfactual thinking
        Qquit = (1 - phi) * Qquit + phi * trialReward * gamma ^ ((iti / stepDuration) + t-1)
      }else{
        Qquit =  (1 - phi) * Qquit + phi *  trialReward
        # counterfactual thinking
        Qquit = (1 - phi) * Qquit + phi * trialReward * gamma ^ ((iti / stepDuration) + t-1)
        if(t > 1){
          Qwait[1 : (t - 1)] = (1 - phi) * Qwait[1 : (t - 1)] +
            phi * trialReward * gamma ^ rev((1 : (t - 1 )))
        }

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