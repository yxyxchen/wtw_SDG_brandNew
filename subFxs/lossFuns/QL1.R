# log likelyhood of paras given choice data 
# QL1: Q-learning model with a single learning rate
# QL2: Q-learning model with separate learning rates for rewards and punishments
# RL1: R-learning model with a single learning rate 
# RL2: R-learning model with separate learning rates for rewards and punishments

# inputs:
# paras: parameter values 
# condition_: [nTrialx1 factor] HP or LP
# trialEarnings_ : [nTrialx1 num] trial-wise rewards
# timeWaited_ : [nTrialx1 real] trial-wise waiting durations 

# outputs
# condition : [nTrialx1 factor] from inputs: 
# LLTrial_ : [nTrialx1 real] log likelyhood of paras given choice data in each trial
# pWaits : [nStepMax x nTrial, real] # the probability of waiting at each step  and each trial


QL1 = function(paras, condition_, trialEarnings_, timeWaited_){
  # extract parameter values
  phi = paras[1]; tau = paras[2]; gamma = paras[3]; prior = paras[4]
  
  # num of trials
  nTrial = length(trialEarnings_) # num of trials 
  
  # parameters for discrete temproal states 
  stepSec = 1 # duration of one time step 
  nStepMax = ifelse(unique(condition_) == "HP", tMaxs[1] / stepSec, tMaxs[2] / stepSec) # maximal number of steps in a trial
  
  # initialize action values 
  Viti = 0.9 * mean(unlist(optimRewardRates) * stepSec / (1 - 0.85))
  Qwaits = (prior  - 1 : nStepMax) * 0.1 + Viti
  
  # initialize output variables
  pWaits_ = matrix(nrow = nStepMax, ncol = nTrial) # the probability of waiting at each step and each trial
  LLTrial_ = vector(length = nTrial); # the loglikelyhood of the input parameters in this trial
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # calculate the probability of waiting at each step 
    pWaits =  sapply(1 : nStepMax, function(i) 1 / (1  + exp((Viti - Qwaits[i])* tau)))
    pWaits_[,tIdx] = pWaits
    
    # number of waiting steps
    trialEarnings = trialEarnings_[tIdx]
    timeWaited = timeWaited_[tIdx]
    T = ceiling(timeWaited / stepSec) + 1
    getReward = trialEarnings != 0
    ## if the agent gets rewards, it waits until the last step. Otherwise
    ## it quits at the final step
    if(getReward ){
      nWait = ceiling(timeWaited / stepSec)
    }else{
      nWait = ceiling(timeWaited / stepSec) - 1
    }
    
    # the loglikelyhood of the input parameters in this trial
    if(getReward){
      junk = log(pWaits[1 : nWait])
      junk[is.na(junk)] = -1000
      LLTrial_[tIdx] = sum(junk)
    }else{
      if(nWait >= 1){
        junk = c(log(pWaits[1 : nWait]), log(1 - pWaits[nWait + 1]))
        junk[is.na(junk)] = -1000
        LLTrial_[tIdx] = sum(junk)
      }else{
        junk = log(1 - pWaits[nWait + 1])
        junk[is.na(junk)] = -1000       
        LLTrial_[tIdx] =  sum(junk)
      }
      
    }
    
    # update action values at the end of each trial
    if(tIdx < nTrial){
      # calculate the reward signal for updating action value which equals 
      # trialEarnings + discounted value of the successor state gamma. Noticably,
      # the successor state at the end of trial is always the iti state before the next trial
      rwdSignal = trialEarnings + Viti * gamma
      # discounted reward signals for step 1 - (T-1)
      stepRwdSignals = sapply(1 : (T-1), function(t) gamma^(T-t-1) * rwdSignal)
      # discounted reward signals for the iti state
      itiRwdSignal = rwdSignal * gamma^(T-2 + iti / stepSec)
      
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
      Viti = Viti + phi * itiDelta
    }# end of the value update section
  } # end of the loop over trials
  
  # return outputs
  outputs = list( 
    "condition" = condition_,
    "LLTrial_" = LLTrial_,
    "pWaits_" = pWaits_
  )
  return(outputs)
}