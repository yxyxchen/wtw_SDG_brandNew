# select modelFun by modelName
getSimModelFun = function(modelName){
  if(modelName == "full_model"){
    simModelFun = full_model
  }else if(modelName == "reduce_one_QwaitIni"){
    simModelFun = reduce_one_QwaitIni
  }else if(modelName == "reduce_two_QwaitIni"){
    simModelFun = reduce_Qwait_Ini
  }else if(modelName == "reduce_one_phi"){
    simModelFun = reduce_one_phi
  }else if(modelName == "reduce_one_gamma"){
    simModelFun == reduce_one_gamma
  }else{
    return("wrong model name!")
  }
  return(simModelFun)
}

################ full_model ######################
full_model = function(para, cond, scheduledWait, nBlock){
  # parse para
  phi = para[1]
  tau = para[2]
  gamma = para[3]
  QwaitIni = para[4]

  # determine simBlockSecs
  simBlockSecs = blockSecs * nBlock
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  Qwait = rep(QwaitIni, nTimeStep)
  Qquit = QwaitIni * gamma^(iti)
  
  # initialize varibles for recording
  vaWaits = matrix(NA, nTimeStep, simBlockSecs / iti + 1);
  vaWaits[,1] = Qwait
  vaQuits = vector(length = simBlockSecs / iti + 1);
  vaQuits[1] = Qquit
  rewardDelays = rep(0, simBlockSecs / iti + 1)
  
  # initialize totalSecs
  totalSecs = 0
  
  # initialize outputs 
  trialEarnings = rep(0, simBlockSecs / iti + 1)
  timeWaited = rep(0, simBlockSecs / iti + 1)
  sellTime = rep(0, simBlockSecs / iti + 1)
  
  # loop over trials
  tIdx = 1
  while(totalSecs < simBlockSecs) {
    # sample rewardDelay
    rewardDelay = drawSample(cond)
    # calculaye available time steps
    # since we use floor there maybe 0.5 sec error (less than 90 s)
    nAvaStep = min(floor((simBlockSecs - totalSecs) / stepDuration), nTimeStep)
    
    # loop over steps
    t = 1
    while(t <= nAvaStep){
      # determine action
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
      action = ifelse(runif(1) < waitRate, 'wait', 'quit')
      # next reward 
      rewardOccur = rewardDelay <= timeTicks[t + 1] && rewardDelay > timeTicks[t] 
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue, 0) 
    
      # dertime whether to continue
      trialGoOn= (action == 'wait' && !rewardOccur && t < nAvaStep)
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
    
    # update Qwait and Qquit and go to the next trail if totalSecs < simBlockSecs
    if(totalSecs < simBlockSecs){
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
      
      # update tIdx
      tIdx = tIdx + 1
    }# end of the update
    
    if(Qquit < 0 || sum(Qwait < 0) > 0){
      browser()
    }
  } # end of all trials 
  
  outputs = list( 
    "trialNum" = 1 : tIdx,
     "trialEarnings" = trialEarnings[1 : tIdx],
     "timeWaited" = timeWaited[1 : tIdx],
     "sellTime" = sellTime[1 : tIdx], # used in wtw analysis
     "scheduledWait" = rewardDelays[1 : tIdx],
     "vaWaits" = vaWaits[, 1 : tIdx],
     "vaQuits" = vaQuits[1 : tIdx]
    )
  return(outputs)
} #end of the function


simulate = function(modelName, nBlock, nRep){
  dir.create("genData/simulation")
  dir.create(sprintf("genData/simulation/%s", modelName))
  # choose modelFun
  simModelFun = getSimModelFun(modelName)
  
  # determine paraComb
  paraTable = data.frame(phi = c(0.02, 0.05, 0.08), tau = c(5, 10, 15),
                         gamma = c(0.85, 0.90, 0.95), QwaitIni = c(2, 3, 4))
  paraComb = getParaComb(paraTable)
  nComb = nrow(paraTable) ^ ncol(paraTable)
  simNo = matrix(seq(1 : nComb * nRep), nrow = nComb, ncol = nRep)
  save("paraComb", "nComb", "nRep", "simNo", file = sprintf("genData/simulation/%s/simParas.RData", modelName))
  # initialize outputs
  trialData = vector(length = nComb * nRep, mode ='list')
  # loop over conditions
  for(condIdx in 1 : 2){
    cond = conditions[condIdx];
    # loop over repetions 
    for(h in 1 : nrow(paraComb)){
      para = paraComb[h,];
      # calculate wIni
      for(j in 1 : nRep ){
        tempt = simModelFun(para, cond, nBlock)
        trialData[[simNo[h, j]]] = tempt
      }  
    }
    # save 
    if(cond == "HP"){
      trialHPData = trialData
      fileName = sprintf('genData/simulation/%s/trialHPData.RData', modelName)
      save(trialHPData,file = fileName)
    }else{
      trialLPData = trialData
      fileName =  sprintf('genData/simulation/%s/trialLPData.RData', modelName)
      save(trialLPData,file = fileName)
    }
  }
}

################ R_learning ######################
R_learning = function(para, cond, nBlock){
  # parse para
  phi1 = para[1]
  phi2 = para[2]
  tau = para[3]
  QwaitIni = para[4]
  costIni = para[5]
  
  # determine simBlockSecs
  simBlockSecs = blockSecs * nBlock
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
}






