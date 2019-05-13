# select modelFun by modelName
getRepModelFun = function(modelName){
  if(modelName == "full_model"){
    repModelFun = full_model
  }else if(modelName == "curiosityTrialRSp"){
    repModelFun = curiosityTrialRSp
  }else if(modelName == "curiosityTrialSp"){
    repModelFun = curiosityTrialSp
  }else if(modelName == "functionRL"){
    repModelFun = functionRL
  }else if(modelName == "functionLinear"){
    repModelFun = functionLinear
  }else if(modelName == "functionParabolic"){
    repModelFun = functionParabolic
  }else if(modelName == "curiosityTrialRUp"){
    repModelFun = curiosityTrialRUp
  }else if(modelName == "curiosityTrialRBack"){
    repModelFun = curiosityTrialRBack
  }else if(modelName == "curiosityTrialRLog"){
    repModelFun = curiosityTrialRLog
  }else if(modelName == "fullModel"){
    repModelFun = fullModel
  }else{
    return("wrong model name!")
  }
  return(repModelFun)
}

################ monte ######################
full_model = function(paras, cond, scheduledWait){
  # parse para
  phi = paras[1]
  tau = paras[2]
  gamma = paras[3]
  QwaitIni = paras[4]
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  nTimeStep = tMax / stepDuration
  
  #Qwait = rep(wIni*0.93, nTimeStep)
  Qwait = rep(QwaitIni, nTimeStep)
  #Qwait = rep(wIni, nTimeStep)
  Qquit = QwaitIni * gamma ^ (iti / stepDuration)
  Viti =  QwaitIni * gamma ^ (iti / stepDuration)
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
    # loop for each timestep t and determine At
    t = 1
    while(t <= nTimeStep){
      # determine At
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
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

################ monte ######################
fullModel = function(paras, cond, scheduledWait){
  # parse para
  phi = paras[1]
  tau = paras[2]
  gamma = paras[3]
  QwaitIni = paras[4]
  
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
  
  # Qwait = rep(wIni, nTimeStep)
  # since the participants start the trial with , we assume max(Qwait0) = wini * 0.8
  # again, we assume it has a slope 
  Qquit = wIni * 0.9
  Viti = wIni * 0.9
  #Qwait = rep(wIni*0.93, nTimeStep)
  Qwait = rep(QwaitIni, nTimeStep)
  #Qwait = rep(wIni, nTimeStep)
  
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
    # loop for each timestep t and determine At
    t = 1
    while(t <= nTimeStep){
      # determine At
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
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

################ curiosityTrial model using the Rlearning  ######################
curiosityTrialRBack = function(paras, cond, scheduledWait){
  # parse para
  phi = paras[1]
  tau = paras[2]
  phiR = paras[3]
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  nTimeStep = tMax / stepDuration
  
  # initialize actionValues
  # here we use the optimal reward rates from the normative analysis in Lempert 2018
  # it is more accurate then the one I calcualte in wtwSettings.R
  # wIni = (5/6 + 0.93)/ 2 * stepDuration # change wIni didn't change anything
  wIni = 0.3
  # wIni = 0 # chaneg to 0 therefore to encourage explore
  
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
  targets = matrix(NA, nTimeStep, nTrial);
  deltas = matrix(NA, nTimeStep, nTrial);
  
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
    # loop for each timestep t and determine At
    t = 1
    while(t <= nTimeStep){
      # determine At
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
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
      targets[1 : (T-1),tIdx] = returns
      # when the agent always wait and get the reward, update Qwait[1:(T-1)]
      # otherwise, update Qquit and Qwait[1 : (T-2)]      
      if(getReward){
        Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])       
        deltas[1 : (T-1),tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
      }else{
        Qquit = Qquit + phi*(returns[T-1] - Qquit)
        if(T > 2){
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phi*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
          deltas[1 : (T-2),tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
        }
      }
      

      # update Viti and Rrate
      deltaIti = returns[1] - Rrate * (iti / stepDuration) - Viti
      Viti = Viti + phi*deltaIti
      Rrate = Rrate + phiR *deltaIti
      # bellman backup
      # Qwait looks better, very close to optimal but not to the behaviors
      phiB = phi
      Qwait[2 : nTimeStep] = (1 - phiB) * Qwait[2: nTimeStep] + phiB* (Qwait[1 : (nTimeStep - 1)]- Rrate)
      Qquit = Qquit + phi*(returns[1] - Rrate * (iti / stepDuration + 1) - Qquit)
      
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
    "Rrates" = Rrates,
    "deltas" = deltas,
    "targets" = targets
  )
  return(outputs)
} #end of the function

################ curiosityTrial model using the Rlearning  ######################
curiosityTrialRSp = function(paras, cond, scheduledWait){
  # parse para
  phi = paras[1]
  tau = paras[2]
  phiR = paras[3]
  zeroPoint = paras[4]

  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  nTimeStep = tMax / stepDuration
  
  # initialize actionValues
  # here we use the optimal reward rates from the normative analysis in Lempert 2018
  # it is more accurate then the one I calcualte in wtwSettings.R
  wIni = (5/6 + 0.93)/ 2 * stepDuration * 0.5# change wIni didn't change anything
  # wIni = 0 # chaneg to 0 therefore to encourage explore
  

  # Qwait = c(5 -(1:35) * 0.2, 5 -(36:40) * 2) # linear with sudden fall
  # Qwait = -(1:40 / 10)^1.5 + 5 # the tail doesn't fall enough
  # Qwait = -(1:20 / 10)^1.5 + 2 # good for HP
  Qwait = zeroPoint*0.2 - 0.2 *(1 : nTimeStep - 1)
  

  # Qwait = -6 / (rev(1:40) + 4) + 0.5
  # plot(-6 / (rev(1:40) + 100) + 0.5) 
  # Qwait = 3.5 -(1:20) * 0.2
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
  targets = matrix(NA, nTimeStep, nTrial);
  deltas = matrix(NA, nTimeStep, nTrial);
  
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
    # loop for each timestep t and determine At
    t = 1
    while(t <= nTimeStep){
      # determine At
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
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
      targets[1 : (T-1),tIdx] = returns
      # when the agent always wait and get the reward, update Qwait[1:(T-1)]
      # otherwise, update Qquit and Qwait[1 : (T-2)]      
      if(getReward){
        Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])       
        deltas[1 : (T-1),tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
      }else{
        Qquit = Qquit + phi*(returns[T-1] - Qquit)
        if(T > 2){
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phi*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
          deltas[1 : (T-2),tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
        }
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
    "Rrates" = Rrates,
    "deltas" = deltas,
    "targets" = targets
  )
  return(outputs)
} #end of the function

################ curiosityTrial model using the Rlearning  ######################
curiosityTrialRLog = function(paras, cond, scheduledWait){
  # parse para
  phi = paras[1]
  tau = paras[2]
  phiR = paras[3]
  
  # 

  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  #nTimeStep = tMax / stepDuration
  base = 1.03
  junk = 1 - tMax / base * (1-base)
  nTimeStep = ceiling(log(junk) / log(base))
  #nTimeStep = 40
  gapEndTimes = cumsum(base ^ (1:nTimeStep))
  gapStartTimes = c(0,cumsum(base ^ (1:(nTimeStep - 1))))
  
  # initialize actionValues
  # here we use the optimal reward rates from the normative analysis in Lempert 2018
  # it is more accurate then the one I calcualte in wtwSettings.R
  wIni = (5/6 + 0.93)/ 2 * stepDuration # change wIni didn't change anything
  # wIni = 0 # chaneg to 0 therefore to encourage explore
  
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
  targets = matrix(NA, nTimeStep, nTrial);
  deltas = matrix(NA, nTimeStep, nTrial);
  
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
    # loop for each timestep t and determine At
    t = 1
    while(t <= nTimeStep){
      # determine At
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
      action = ifelse(runif(1) < waitRate, 'wait', 'quit')
      # observe St+1 and Rt+1
      rewardOccur = thisScheduledWait <= gapEndTimes[t] 
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue, 0) 
      
      # dertime whether St+1 is the terminal state
      # if the trial terminates, track terminal timestep index T, trialEarnings, timeWaited, sellTime and elapsedTime
      # otherwise, continue
      nextStateTerminal = (getReward || action == "quit")
      if(nextStateTerminal){
        T = t+1
        trialEarnings[tIdx] = ifelse(nextReward == tokenValue, tokenValue, 0);
        timeWaited[tIdx] = ifelse(getReward, thisScheduledWait, gapEndTimes[t])
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
      targets[1 : (T-1),tIdx] = returns
      # when the agent always wait and get the reward, update Qwait[1:(T-1)]
      # otherwise, update Qquit and Qwait[1 : (T-2)]      
      if(getReward){
        Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])       
        deltas[1 : (T-1),tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
      }else{
        Qquit = Qquit + phi*(returns[T-1] - Qquit)
        if(T > 2){
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phi*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
          deltas[1 : (T-2),tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
        }
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
    "Rrates" = Rrates,
    "deltas" = deltas,
    "targets" = targets
  )
  return(outputs)
} #end of the function

################ curiosityTrial model using the Rlearning  ######################
curiosityTrialRUp = function(paras, cond, scheduledWait){
  # parse para
  phi = paras[1]
  tau = paras[2]
  phiR = paras[3]
  c = paras[4]
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  nTimeStep = tMax / stepDuration
  
  # initialize actionValues
  # here we use the optimal reward rates from the normative analysis in Lempert 2018
  # it is more accurate then the one I calcualte in wtwSettings.R
  # wIni = (5/6 + 0.93)/ 2 * stepDuration # change wIni didn't change anything
  wIni = 0 # chaneg to 0 therefore to encourage explore
  
  Qwait = c(rep(1, nTimeStep))
  nUpdate = c(rep(0, nTimeStep)) # count of update
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
  targets = matrix(NA, nTimeStep, nTrial);
  deltas = matrix(NA, nTimeStep, nTrial);
  nUpdates = matrix(NA, nTimeStep, nTrial);
  nUpdates[,1] = nUpdate
  uppers = matrix(NA, nTimeStep, nTrial);
  
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
    # loop for each timestep t and determine At
    t = 1
    while(t <= nTimeStep){
      # determine At
      if(nUpdate[t] == 0) waitRate = 1
      else waitRate =  1 / sum(1  + exp(Qquit - Qwait[t] - c * (sqrt(log(tIdx) / nUpdate[t]) -
                                                                        sqrt(log(tIdx) / (tIdx-1)))))
      uppers[t, tIdx] = sqrt(log(tIdx) / nUpdate[t]) - sqrt(log(tIdx) / (tIdx-1))
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
      
      # update nUpdate
      if(getReward) nUpdate[1 : (T - 1)] = nUpdate[1 : (T - 1)] + 1
      # update action values for each timestep t
      returns = sapply(1 : (T-1), function(t) nextReward - (T-t) * Rrate + Viti)
      targets[1 : (T-1),tIdx] = returns
      # when the agent always wait and get the reward, update Qwait[1:(T-1)]
      # otherwise, update Qquit and Qwait[1 : (T-2)]      
      if(getReward){
        Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])       
        deltas[1 : (T-1),tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
      }else{
        Qquit = Qquit + phi*(returns[T-1] - Qquit)
        if(T > 2){
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phi*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
          deltas[1 : (T-2),tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
        }
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
      nUpdates[,tIdx + 1] = nUpdate
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
    "Rrates" = Rrates,
    "deltas" = deltas,
    "targets" = targets,
    "nUpdates" = nUpdates,
    "uppers" = uppers
  )
  return(outputs)
} #end of the function
################ monte ######################
curiosityTrialSp = function(paras, cond, scheduledWait){
  # parse para
  phi = paras[1]
  tau = paras[2]
  gamma = paras[3]
  zeroPoint = paras[4]

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
  
  # Qwait = rep(wIni, nTimeStep)
  # since the participants start the trial with , we assume max(Qwait0) = wini * 0.8
  # again, we assume it has a slope 
  Qquit = wIni * 0.9
  Viti = wIni * 0.9
  #Qwait = rep(wIni*0.93, nTimeStep)
  Qwait = zeroPoint*0.1 - 0.1*(0 : (nTimeStep - 1)) + Qquit
  #Qwait = rep(wIni, nTimeStep)

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
    # loop for each timestep t and determine At
    t = 1
    while(t <= nTimeStep){
      # determine At
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
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

################ monte ######################
curiosityTrialSp = function(paras, cond, scheduledWait){
  # parse para
  phi = paras[1]
  tau = paras[2]
  gamma = paras[3]
  zeroPoint = paras[4]
  
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
  
  # Qwait = rep(wIni, nTimeStep)
  # since the participants start the trial with , we assume max(Qwait0) = wini * 0.8
  # again, we assume it has a slope 
  Qquit = wIni * 0.9
  Viti = wIni * 0.9
  #Qwait = rep(wIni*0.93, nTimeStep)
  Qwait = zeroPoint*0.1 - 0.1*(0 : (nTimeStep - 1)) + Qquit
  #Qwait = rep(wIni, nTimeStep)
  
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
    # loop for each timestep t and determine At
    t = 1
    while(t <= nTimeStep){
      # determine At
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
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


 
functionRL =  function(paras, cond, scheduledWait){
  # parse paras
  phi = paras[1]
  tau = paras[2]
  
  # parse the data
  nTrial = length(scheduledWait)
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  nTimeStep = tMax / stepDuration
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)

  # initialize elapsed time
  elapsedTime = 0
  
  # initialize varibales for debugs 
  rrs = matrix(NA, nTimeStep, nTrial);
  rrBars = vector(length = nTrial);
  deltas = matrix(NA, nTimeStep, nTrial)
  targets = matrix(NA, nTimeStep, nTrial)
  
  rrIni = 1.2
  rr = rep(rrIni, nTimeStep)
  rrBar = ifelse(cond == "HP", 5 / 6, 0.93)
  rrs[,1] = rr
  rrBars[1] = rrBar
  
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # determine 
    thisScheduledWait = scheduledWait[tIdx]
    # clear cached values
    target = rep(0, nTimeStep)
    delta = rep(0, nTimeStep)
    # loop for each timestep t and determine At
    t = 1
    while(t <= nTimeStep){
      # determine At
      waitRate =  1 / sum(1  + exp((rrBar - rr[t])* tau))
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
        trialEarnings[tIdx] = nextReward;
        timeWaited[tIdx] = ifelse(getReward, thisScheduledWait, t * stepDuration)
        sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
        elapsedTime = elapsedTime + timeWaited[tIdx] + iti
        break
      }else{
        t = t + 1
      }
    }# end of the action selection section
    
    # update reward rates 
    if(tIdx < nTrial){
      target[1 : (T-1)]= rev(nextReward / (1 : (T-1) * stepDuration))
      delta[1 : (T-1)] = target[1 : (T-1)] - rr[1 : (T-1)]
      rr = rr + phi*delta
      rrBar = rrBar + phi * (nextReward / ((T-1) * stepDuration + iti) - rrBar)
      
      # record updated values
      rrs[,tIdx + 1] = rr
      rrBars[tIdx + 1] = rrBar
      deltas[,tIdx] = delta
      targets[,tIdx] = target
    }# end of the value update section
  } # end of the trial loop
  
  outputs = list( 
    "trialNum" = 1 : nTrial,
    "trialEarnings" = trialEarnings,
    "timeWaited" = timeWaited,
    "sellTime" = sellTime, # used in wtw analysis
    "scheduledWait" = scheduledWait,
    "rrs" = rrs,
    "rrBars" = rrBars,
    "targets" = targets,
    "deltas" = deltas
  )
  return(outputs)
}


functionLinear =  function(paras, cond, scheduledWait){
  # parse the data
  nTrial = length(scheduledWait)
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  nTimeStep = tMax / stepDuration
  
  # parse paras
  phi = paras[1]
  tau = paras[2]
  rrIni = paras[3]
  sigma = paras[4]
  
  yIni = 0
  
  beta.prior = matrix(c(rrIni, -0.05), ncol = 1)
  sigmaSq.prior = diag(c(sigma, sigma))
  x.star = cbind(rep(1, nTimeStep), (1 : nTimeStep) * stepDuration)
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  Ts = rep(0, nTrial)
  
  # initialize elapsed time
  elapsedTime = 0
  
  # initialize varibales for debugs 
  rrs = matrix(NA, nTimeStep, nTrial);
  ws = matrix(NA, nTimeStep, nTrial);
  ys = matrix(NA, nTimeStep, nTrial);
  rrBars = vector(length = nTrial);
  betaPosts =  matrix(NA, 2, nTrial);
  betaLiks =  matrix(NA, 2, nTrial);
  sigmaSqs =  vector(length = nTrial);
  
  # initialize cached values
  rr = x.star %*% beta.prior
  rrBar = ifelse(cond == "HP", 5 / 6, 0.93)
  w = rep(0, nTimeStep)
  x.mu = (1 : nTimeStep) * stepDuration
  y.mu = rep(yIni, nTimeStep)
  rrs[,1] = rr
  rrBars[1] = rrBar

  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # determine 
    thisScheduledWait = scheduledWait[tIdx]
    # clear cached values
    target = rep(0, nTimeStep)
    delta = rep(0, nTimeStep)
    # loop for each timestep t and determine At
    t = 1
    while(t <= nTimeStep){
      # determine At
      waitRate =  1 / sum(1  + exp((rrBar - rr[t])* tau))
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
        Ts[tIdx] = T
        trialEarnings[tIdx] = nextReward;
        timeWaited[tIdx] = ifelse(getReward, thisScheduledWait, t * stepDuration)
        sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
        elapsedTime = elapsedTime + timeWaited[tIdx] + iti
        break
      }else{
        t = t + 1
      }
    }# end of the action selection section
    
    # update reward rates 
    if(tIdx < nTrial){
      # updated averaged observation
      if(T > 2){
        r = trialEarnings[tIdx] / rev(stepDuration * (1 : (T - 1)))
        w[1 : (T - 1)] = w[1 : (T - 1)]  + 1
        y.mu[1 : (T - 1)] = y.mu[1 : (T - 1)] + 1 / w[1 : (T - 1)] * (r - y.mu[1 : (T - 1)] )
        x = unlist(lapply(1 : nTimeStep, function(i) rep(x.mu[i], w[i])))
        y = unlist(lapply(1 : nTimeStep, function(i) rep(y.mu[i], w[i])))
        X = cbind(rep(1, length(x)), x)
        Y = matrix(y, ncol = 1)
        
        # Bayesian 
        betaLik = solve(t(X) %*% X) %*% t(X) %*% Y
        sigmaSq = sum((y - X %*% betaLik)^2) / (length(y) - 2) + 10 / sqrt(tIdx)
        A = solve(sigmaSq.prior) + t(X) %*% X / sigmaSq 
        betaPost = solve(A)  %*% (t(X) %*% Y/sigmaSq + solve(sigmaSq.prior) %*% beta.prior)
        yHat.mu = x.star %*% betaPost
        rr = yHat.mu
        rrBar = rrBar + phi * (nextReward / ((T-1) * stepDuration + iti) - rrBar)
      }
      # record updated values
      rrs[,tIdx + 1] = rr
      rrBars[tIdx + 1] = rrBar
      betaPosts[,tIdx] = betaPost
      betaLiks[,tIdx] = betaLik
      sigmaSqs[tIdx] = sigmaSq
      ys[,tIdx] = y.mu
      ws[,tIdx] = w
    }# end of the value update section
  } # end of the trial loop
  
  outputs = list( 
    "trialNum" = 1 : nTrial,
    "trialEarnings" = trialEarnings,
    "timeWaited" = timeWaited,
    "sellTime" = sellTime, # used in wtw analysis
    "scheduledWait" = scheduledWait,
    "rrs" = rrs,
    "rrBars" = rrBars,
    "betaPosts" = betaPosts,
    "betaLiks" = betaLiks,
    "sigmaSqs" = sigmaSqs,
    "ys" = ys,
    "ws" = ws
  )
  return(outputs)
}

