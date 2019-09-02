# simulate behavioral data for a single participant 
# using the given parameters and the empirical reward delay sequence

# inputs:
# paras: parameter values 
# condition: HP or LP
# scheduledWait : empirical trial-wise reward delay sequences 

# outputs



modelSimSingle = function(modelName, paras, scheduledWait){
  if(modelName == "QL1") genModel  = QL1
  else if(modelName == "QL2") genModel = QL2
  else if(modelName == "RL1") genModel  = RL1
  else if(modelName == "RL2") genModel  = RL2
  else if(modelName == "BL") genModel  = BL
  else{
    return("wrong model name!")
  }
  return(genModel)
}

# this function replicates the behavioral data using inividual fitted parameters
modelRepitation = function(modelName, summaryData, trialData,  nComb){
  # determine repFun
  repFun = getRepFun(modelName)
  
  # load inividual fitted parameters
  paraNames = getParaNames(modelName)
  parentDir ="genData/expModelFitting"; dirName = sprintf("%s/%sdb",parentDir, modelName)
  expPara = loadExpPara(paraNames, dirName)
  ids = expPara$id; nSub = length(ids)
  
  # initialize outputs
  repTrialData = vector(length = nSub * nComb, mode ='list')
  repNo = matrix(1 : (nSub * nComb), nrow = nComb, ncol = nSub)
  
  # draw nComb parameter combination samples, and simualte nComb times for each participant
  set.seed(231)
  for(sIdx in 1 : nSub){
    # prepare inputs
    id = ids[[sIdx]]
    paras_ = read.table(sprintf("%s/%sdb/s%s.txt", parentDir, modelName, id),sep = ",", row.names = NULL)
    thisTrialData = trialData[[id]] # here we id instead of sIdx
    # excluded some trials
    excluedTrials1 = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[1]) &
                             thisTrialData$condition == conditions[1])
    excluedTrials2 = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[2]) &
                             thisTrialData$condition == conditions[2])
    excluedTrials = c(excluedTrials1, excluedTrials2)
    thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excluedTrials,]
    cond = unique(thisTrialData$condition)
    scheduledWait = thisTrialData$scheduledWait
    
    # simulate nComb times
    for(cbIdx in 1 : nComb){
      paras = as.double(paras_[sample(1 : nrow(paras_), 1), 1 : length(paraNames)])
      tempt = repFun(paras, cond, scheduledWait)
      repTrialData[[repNo[cbIdx, sIdx]]] = tempt
    }
  }
  outputs = list(expPara = expPara, repTrialData = repTrialData, repNo = repNo)
  return(outputs)
}


QL1 = function(paras, condition, scheduledWait){
  # extract parameter values
  phi = paras[1]; tau = paras[2]; gamma = paras[3]; prior = paras[4]
  
  # other constants 
  nTrial = length(scheduledWait) # num of trials 
  stepSec = 1 # duration of one time step 
  nStepMax = ifelse(unique(condition) == "HP", tMaxs[1] / stepSec, tMaxs[2] / stepSec) # maximal number of steps in a trial
  
  # initialize action values 
  Viti = 0.9 * mean(unlist(optimRewardRates) * stepSec / (1 - 0.85))
  Qwaits = (prior  - 1 : nStepMax) * 0.1 + Viti
  
  # initialize output variables
  Qwaits_ = matrix(NA, nTimeStep, nTrial); Qwaits_[,1] = Qwaits
  Viti_ = vector(length = nTrial); Viti_[1] = Viti
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  
  # track elpased time from the beginning of the task 
  elapsedTime = 0 
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    
    # loop over steps until a trial ends 
    t = 1
    while(t <= nStepMax){
      # determine the current action At 
      pWait =  1 / sum(1  + exp((Viti - Qwaits[t])* tau))
      action = ifelse(runif(1) < pWait, 'wait', 'quit')
      
      # if a reward occurs at t and the agent is still waiting, Rt+1 = 10 otherwise Rt+1 = 0
      rewardOccur = scheduledWait[tIdx] <= (t * stepSec) && scheduledWait[tIdx] > ((t-1) * sepSec)
      getReward = (action == 'wait' && rewardOccur)
      nextReward = ifelse(getReward, tokenValue, 0) 
      
      # a trial ends if the agent gets a reward or quits. 
      # if the trial ends, record relavant variables.Otherwise, proceed to the next step 
      isTerminal = (getReward || action == "quit")
      if(isTerminal){
        # record the terminal state T 
        T = t + 1 
        # record output variables 
        trialEarnings[tIdx] = nextReward
        timeWaited[tIdx] = ifelse(getReward, scheduledWait[tIdx], t * stepSec)
        sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
        # update the elapsed time
        elapsedTime = elapsedTime + timeWaited[tIdx] + iti
        break
      }else{
        t = t + 1
      }
    }
    
    # update action values at the end of each trial
    if(tIdx < nTrial){
      # calculate the reward signal for updating action value 
      # which equals Rt+1 + V(St+1) * gamma. Noticably, St+1 is the iti state
      rwdSignal = nextReward + Viti * gamma
      # discounted reward signals for step 1 - (T-1)
      stepRwdSignals = sapply(1 : (T-1), function(t) gamma^(T-t-1) * rwdSignal)
      # discounted reward signals for the iti state
      itiRwdSignal = rwdSignal * gamma^(T-2 + iti / stepDuration)
      
      # update Qwaits 
      if(getReward){
        Qwaits[1 : (T-1)] = Qwaits[1 : (T-1)] + phi *(discReturns[1 : (T-1)] - Qwaits[1 : (T-1)])
      }else{
        if(T > 2){
          # in non-rewarded trials, Qwait in the last step will not be updated
          # since the agent chooses to quit on that step 
          Qwaits[1 : (T-2)] = Qwaits[1 : (T-2)] + phi * (discReturns[1 : (T-2)] - Qwaits[1 : (T-2)])
        }
      }
      
      # update Viti
      itiDelta = itiRwdSignal - Viti
      Viti =  Viti + phi * itiDelta
      
      # record updated values
      Qwaits_[,tIdx + 1] = Qwaits
      Viti_[tIdx + 1] = Viti
    }# end of the value update section
  } # end of the loop over trials
  
  # return outputs
  outputs = list( 
    "trialNum" = 1 : nTrial, 
    "condition" = condition,
    "trialEarnings" = trialEarnings, "timeWaited" = timeWaited,
    "sellTime" = sellTime, "scheduledWait" = scheduledWait,
    "Qwaits_" = Qwait, 
    "Viti_" = Viti_
  )
  return(outputs)
}

QL2 = function(paras, condition, scheduledWait){
  # extract parameter values
  phi_pos = paras[1]; phi_neg = paras[2]; tau = paras[3]; gamma = paras[4]; prior = paras[5]
  
  # other constants 
  nTrial = length(scheduledWait) # num of trials 
  stepSec = 1 # duration of one time step 
  nStepMax = ifelse(unique(condition) == "HP", tMaxs[1] / stepSec, tMaxs[2] / stepSec) # maximal number of steps in a trial
  
  # initialize action values 
  Viti = 0.9 * mean(unlist(optimRewardRates) * stepSec / (1 - 0.85))
  Qwaits = (prior  - 1 : nStepMax) * 0.1 + Viti
  
  # initialize output variables
  Qwaits_ = matrix(NA, nTimeStep, nTrial); Qwaits_[,1] = Qwaits
  Viti_ = vector(length = nTrial); Viti_[1] = Viti
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  
  # track elpased time from the beginning of the task 
  elapsedTime = 0 

  # loop over trials
  for(tIdx in 1 : nTrial) {
    
    # loop over steps until a trial ends 
    t = 1
    while(t <= nStepMax){
      # determine the current action At 
      pWait =  1 / sum(1  + exp((Viti - Qwaits[t])* tau))
      action = ifelse(runif(1) < pWait, 'wait', 'quit')
    
      # if a reward occurs at t and the agent is still waiting, Rt+1 = 10 otherwise Rt+1 = 0
      rewardOccur = scheduledWait[tIdx] <= (t * stepSec) && scheduledWait[tIdx] > ((t-1) * sepSec)
      getReward = (action == 'wait' && rewardOccur)
      nextReward = ifelse(getReward, tokenValue, 0) 
      
      # a trial ends if the agent gets a reward or quits. 
      # if the trial ends, record relavant variables.Otherwise, proceed to the next step 
      isTerminal = (getReward || action == "quit")
      if(isTerminal){
        # record the terminal state T 
        T = t + 1 
        # record output variables 
        trialEarnings[tIdx] = nextReward
        timeWaited[tIdx] = ifelse(getReward, scheduledWait[tIdx], t * stepSec)
        sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
        # update the elapsed time
        elapsedTime = elapsedTime + timeWaited[tIdx] + iti
        break
      }else{
        t = t + 1
      }
    }
    
    # update action values at the end of each trial
    if(tIdx < nTrial){
      # calculate the reward signal for updating action value 
      # which equals Rt+1 + V(St+1) * gamma. Noticably, St+1 is the iti state
      rwdSignal = nextReward + Viti * gamma
      # discounted reward signals for step 1 - (T-1)
      stepRwdSignals = sapply(1 : (T-1), function(t) gamma^(T-t-1) * rwdSignal)
      # discounted reward signals for the iti state
      itiRwdSignal = rwdSignal * gamma^(T-2 + iti / stepDuration)
        
      # update Qwaits 
      if(getReward){
        Qwaits[1 : (T-1)] = Qwaits[1 : (T-1)] + phi_pos *(discReturns[1 : (T-1)] - Qwaits[1 : (T-1)])
      }else{
        if(T > 2){
          # in non-rewarded trials, Qwait in the last step will not be updated
          # since the agent chooses to quit on that step 
          Qwaits[1 : (T-2)] = Qwaits[1 : (T-2)] + phi_neg * (discReturns[1 : (T-2)] - Qwaits[1 : (T-2)])
        }
      }
      
      # update Viti
      itiDelta = itiRwdSignal - Viti
      Viti = ifelse(nextReward > 0, Viti + phi_pos * itiDelta, Viti + phi_neg * itiDelta)
      
      # record updated values
      Qwaits_[,tIdx + 1] = Qwaits
      Viti_[tIdx + 1] = Viti
    }# end of the value update section
  } # end of the loop over trials
  
  # return outputs
  outputs = list( 
    "trialNum" = 1 : nTrial, 
    "condition" = condition,
    "trialEarnings" = trialEarnings, "timeWaited" = timeWaited,
    "sellTime" = sellTime, "scheduledWait" = scheduledWait,
    "Qwaits_" = Qwait, 
    "Viti_" = Viti_
  )
  return(outputs)
}

RL1 = function(paras, condition, scheduledWait){
  # extract parameter values
  phi = paras[1]; tau = paras[2]; prior = paras[3]; beta = paras[4];
  
  # other constants 
  nTrial = length(scheduledWait) # num of trials 
  stepSec = 1 # duration of one time step 
  nStepMax = ifelse(unique(condition) == "HP", tMaxs[1] / stepSec, tMaxs[2] / stepSec) # maximal number of steps in a trial
  
  # initialize action values 
  Viti = 0.9 * mean(unlist(optimRewardRates) * stepSec / (1 - 0.85))
  Qwaits = (prior  - 1 : nStepMax) * 0.1 + Viti
  
  # initialize output variables
  Qwaits_ = matrix(NA, nTimeStep, nTrial); Qwaits_[,1] = Qwaits
  Viti_ = vector(length = nTrial); Viti_[1] = Viti
  reRate_ = vector(length = nTrial); reRate_[1] = reRate
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  
  # track elpased time from the beginning of the task 
  elapsedTime = 0 
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    
    # loop over steps until a trial ends 
    t = 1
    while(t <= nStepMax){
      # determine the current action At 
      pWait =  1 / sum(1  + exp((Viti - Qwaits[t])* tau))
      action = ifelse(runif(1) < pWait, 'wait', 'quit')
      
      # if a reward occurs at t and the agent is still waiting, Rt+1 = 10 otherwise Rt+1 = 0
      rewardOccur = scheduledWait[tIdx] <= (t * stepSec) && scheduledWait[tIdx] > ((t-1) * sepSec)
      getReward = (action == 'wait' && rewardOccur)
      nextReward = ifelse(getReward, tokenValue, 0) 
      
      # a trial ends if the agent gets a reward or quits. 
      # if the trial ends, record relavant variables.Otherwise, proceed to the next step 
      isTerminal = (getReward || action == "quit")
      if(isTerminal){
        # record the terminal state T 
        T = t + 1 
        # record output variables 
        trialEarnings[tIdx] = nextReward
        timeWaited[tIdx] = ifelse(getReward, scheduledWait[tIdx], t * stepSec)
        sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
        # update the elapsed time
        elapsedTime = elapsedTime + timeWaited[tIdx] + iti
        break
      }else{
        t = t + 1
      }
    }
    
    # update action values at the end of each trial
    if(tIdx < nTrial){
      # calculate the reward signal for updating action value 
      # which equals Rt+1 + V(St+1). Noticably, St+1 is the iti state
      rwdSignal = nextReward + Viti 
      # net reward signals for step 1 - (T-1)
      stepRwdSignals = sapply(1 : (T-1), function(t) rwdSignal - (T-t) * reRate)
      # net reward signals for the iti state
      itiRwdSignal = rwdSignal - (T - 1 - iti / stepSec) * reRate
      
      # update Qwaits 
      if(getReward){
        Qwaits[1 : (T-1)] = Qwaits[1 : (T-1)] + phi *(stepRwdSignals[1 : (T-1)] - Qwaits[1 : (T-1)])
      }else{
        if(T > 2){
          # in non-rewarded trials, Qwait in the last step will not be updated
          # since the agent chooses to quit on that step 
          Qwaits[1 : (T-2)] = Qwaits[1 : (T-2)] + phi * (stepRwdSignals[1 : (T-2)] - Qwaits[1 : (T-2)])
        }
      }
      
      # update Viti
      itiDelta = itiRwdSignal - Viti
      Viti =  Viti + phi * itiDelta
      reRate = reRate + beta * itiDelta
      
      # record updated values
      Qwaits_[,tIdx + 1] = Qwaits
      Viti_[tIdx + 1] = Viti
      reRate_[tIdx + 1] = reRate 
    }# end of the value update section
  } # end of the loop over trials
  
  # return outputs
  outputs = list( 
    "trialNum" = 1 : nTrial, 
    "condition" = condition,
    "trialEarnings" = trialEarnings, "timeWaited" = timeWaited,
    "sellTime" = sellTime, "scheduledWait" = scheduledWait,
    "Qwaits_" = Qwait, 
    "Viti_" = Viti_,
    "Rerate_" = "Rerate_"
  )
  return(outputs)
}

RL2 = function(paras, condition, scheduledWait){
  # extract parameter values
  phi_pos = paras[1]; phi_neg = paras[2]; tau = paras[3]; prior = paras[4]; beta = paras[5];
  
  # other constants 
  nTrial = length(scheduledWait) # num of trials 
  stepSec = 1 # duration of one time step 
  nStepMax = ifelse(unique(condition) == "HP", tMaxs[1] / stepSec, tMaxs[2] / stepSec) # maximal number of steps in a trial
  
  # initialize action values 
  Viti = 0.9 * mean(unlist(optimRewardRates) * stepSec / (1 - 0.85))
  Qwaits = (prior  - 1 : nStepMax) * 0.1 + Viti
  
  # initialize output variables
  Qwaits_ = matrix(NA, nTimeStep, nTrial); Qwaits_[,1] = Qwaits
  Viti_ = vector(length = nTrial); Viti_[1] = Viti
  reRate_ = vector(length = nTrial); reRate_[1] = reRate
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  
  # track elpased time from the beginning of the task 
  elapsedTime = 0 
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    
    # loop over steps until a trial ends 
    t = 1
    while(t <= nStepMax){
      # determine the current action At 
      pWait =  1 / sum(1  + exp((Viti - Qwaits[t])* tau))
      action = ifelse(runif(1) < pWait, 'wait', 'quit')
      
      # if a reward occurs at t and the agent is still waiting, Rt+1 = 10 otherwise Rt+1 = 0
      rewardOccur = scheduledWait[tIdx] <= (t * stepSec) && scheduledWait[tIdx] > ((t-1) * sepSec)
      getReward = (action == 'wait' && rewardOccur)
      nextReward = ifelse(getReward, tokenValue, 0) 
      
      # a trial ends if the agent gets a reward or quits. 
      # if the trial ends, record relavant variables.Otherwise, proceed to the next step 
      isTerminal = (getReward || action == "quit")
      if(isTerminal){
        # record the terminal state T 
        T = t + 1 
        # record output variables 
        trialEarnings[tIdx] = nextReward
        timeWaited[tIdx] = ifelse(getReward, scheduledWait[tIdx], t * stepSec)
        sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
        # update the elapsed time
        elapsedTime = elapsedTime + timeWaited[tIdx] + iti
        break
      }else{
        t = t + 1
      }
    }
    
    # update action values at the end of each trial
    if(tIdx < nTrial){
      # calculate the reward signal for updating action value 
      # which equals Rt+1 + V(St+1). Noticably, St+1 is the iti state
      rwdSignal = nextReward + Viti 
      # net reward signals for step 1 - (T-1)
      stepRwdSignals = sapply(1 : (T-1), function(t) rwdSignal - (T-t) * reRate)
      # net reward signals for the iti state
      itiRwdSignal = rwdSignal - (T - 1 - iti / stepSec) * reRate
      
      # update Qwaits 
      if(getReward){
        Qwaits[1 : (T-1)] = Qwaits[1 : (T-1)] + phi_pos *(stepRwdSignals[1 : (T-1)] - Qwaits[1 : (T-1)])
      }else{
        if(T > 2){
          # in non-rewarded trials, Qwait in the last step will not be updated
          # since the agent chooses to quit on that step 
          Qwaits[1 : (T-2)] = Qwaits[1 : (T-2)] + phi_neg * (stepRwdSignals[1 : (T-2)] - Qwaits[1 : (T-2)])
        }
      }
      
      # update Viti
      itiDelta = itiRwdSignal - Viti
      Viti =  ifelse(nextReward > 0, Viti + phi_pos * itiDelta, Viti + phi_neg * itiDelta)
      reRate = reRate + beta * itiDelta
      
      # record updated values
      Qwaits_[,tIdx + 1] = Qwaits
      Viti_[tIdx + 1] = Viti
      reRate_[tIdx + 1] = reRate 
    }# end of the value update section
  } # end of the loop over trials
  
  # return outputs
  outputs = list( 
    "trialNum" = 1 : nTrial, 
    "condition" = condition,
    "trialEarnings" = trialEarnings, "timeWaited" = timeWaited,
    "sellTime" = sellTime, "scheduledWait" = scheduledWait,
    "Qwaits_" = Qwait, 
    "Viti_" = Viti_,
    "Rerate_" = "Rerate_"
  )
  return(outputs)
}

BL = function(paras, cond, scheduledWait){
  # parse 
  pWait = paras[1]
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  nTimeStep = tMax / stepDuration
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial); timeWaited = rep(0, nTrial); sellTime = rep(0, nTrial); elapsedTime = 0
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # action selections
    t = 1
    thisScheduledWait = scheduledWait[tIdx]
    while(t <= nTimeStep){
      # determine At
      action = ifelse(runif(1) < pWait, 'wait', 'quit')
      # observe St+1 and Rt+1
      rewardOccur = thisScheduledWait <= (t * stepDuration) && thisScheduledWait > ((t-1) * stepDuration)
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue, 0) 
      
      # dertime whether St+1 is the terminal state
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
  }
  
  # return outputs
  outputs = list( 
    "trialNum" = 1 : nTrial, "trialEarnings" = trialEarnings, "timeWaited" = timeWaited,
    "sellTime" = sellTime, "scheduledWait" = scheduledWait, "condition" = cond
  )
  return(outputs)
}


