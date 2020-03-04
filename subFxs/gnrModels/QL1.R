# our reinfocement learning generative models simulate persistence behavior as wait-or-quit choices 
# on discrete time steps:
# QL1: Q-learning model with a single learning rate
# QL2: Q-learning model with separate learning rates for rewards and non-rewards
# RL1: R-learning model with a single learning rate 
# RL2: R-learning model with separate learning rates for rewards and non-rewards

# inputs:
# paras: parameter values 
# condition_: [nTrialx1 factor] HP or LP
# scheduledWait_ : [nTrialx1 num] trial-wise reward delays

# outputs
# trialNum : [nTrialx1 int] 1 : nTrial
# condition : [nTrialx1 factor] from inputs
# scheduledWait : [nTrialx1 num] from inputs 
# trialEarnings : [nTrialx1 int] payment for each trial, either 10 or 0
# timeWaited : [nTrialx1 num] waiting duration for each trial 
# sellTime : [nTrialx1 num]  when the agent sells then token on each trial 
# Qwaits_ : [nStepMax x nTrial num] action value of waiting for each time step (2<= t <= nStepMax)at each trial
# V0_ : [nTrialx1 num] state value of the iti stage at each trial 

QL1 = function(paras, condition_, scheduledWait_){
  # extract learning parameter 
  alpha = paras[1]; tau = paras[2]; gamma = paras[3]; prior = paras[4]
  
  # num of trials
  nTrial = length(scheduledWait_) # num of trials 
  
  # parameters for discrete temproal states 
  stepSec = 1 # duration of one time step 
  nStepMax = ifelse(unique(condition_) == "HP", (tMaxs[1] + iti) / stepSec, (tMaxs[2] + iti) / stepSec) # maximal number of steps in a trial
  
  # initialize action values 
  V0 = mean(unlist(optimRewardRates)) * stepSec / (1 - 0.85) # state value for t = 0
  Qwaits = -0.1 * (0 : (nStepMax-1)) * stepSec + prior + gamma * V0 # action values for 0 <= t < tMax
  
  # initialize output variables
  Qwaits_ = matrix(NA, (nStepMax - 2), nTrial); Qwaits_[,1] = Qwaits 
  V0_ = vector(length = nTrial); V0_[1] = V0
  trialEarnings_ = rep(0, nTrial)
  timeWaited_ = rep(0, nTrial)
  sellTime_ = rep(0, nTrial)
  
  # track elpased time from the beginning of the task 
  elapsedTime = 0 
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # current scheduledWait 
    scheduledWait = scheduledWait_[tIdx]
    
    # loop over steps until a trial ends
    step = 0
    while(step <= nStepMax){
      t = (step - 1) * stepSec # this time step [t, t + stepSec)
      # take actions after the iti
      if(t >= iti){
        # decide whether to wait or quit
        pWait =  1 / sum(1  + exp((V0 - Qwaits[step])* tau))
        action = ifelse(runif(1) < pWait, 'wait', 'quit')
        
        # if a reward occurs and the agent is still waiting, the agent gets the reward
        rewardOccur = ((scheduledWait + iti) >= t) & (scheduledWait + iti) < (t + stepSec)
        getReward = (action == 'wait' && rewardOccur)
        
        # a trial ends if the agent gets a reward or quits. if the trial ends, 
        # proceed to the next trial.Otherwise, proceed to the next step 
        isTerminal = (getReward || action == "quit")
        if(isTerminal){
          # update trial-wise variables 
          T = step # terminal step
          trialEarnings = ifelse(getReward, tokenValue, 0) # payment 
          timeWaited =  ifelse(getReward, scheduledWait, t + stepSec - iti) # waiting duration 
          sellTime = elapsedTime + timeWaited # elapsed task time when the agent sells the token
          elapsedTime = elapsedTime + timeWaited + iti  # elapsed task time before the next trial
          # record trial-wise variables
          trialEarnings_[tIdx] = trialEarnings
          timeWaited_[tIdx] = timeWaited
          sellTime_[tIdx] = sellTime
          break
        }else{
          t = t + 1
        }
      }

    }
    
    # update action values at the end of each trial
    if(tIdx < nTrial){
      # calculate the reward signal for updating action value which equals 
      # trialEarnings + discounted value of the successor state gamma. Noticably,
      # the successor state at the end of trial is always the iti state before the next trial
      rwdSignal = trialEarnings + V0 * gamma
      # discounted reward signals for step 1 - (T-1)
      stepRwdSignals = sapply(1 : (T-1), function(t) gamma^(T-t-1) * rwdSignal)
      # discounted reward signals for the iti state
      itiRwdSignal = rwdSignal * gamma^(T-2 + iti / stepSec)
      
      # update Qwaits
      Qwaits[1 : T] = 
      
      # update Qwaits 
      if(getReward){
        Qwaits[1 : (T-1)] = Qwaits[1 : (T-1)] + alpha *(stepRwdSignals[1 : (T-1)] - Qwaits[1 : (T-1)])
      }else{
        if(T > 2){
          # in non-rewarded trials, Qwait in the last step will not be updated
          # since the agent chooses to quit on that step 
          Qwaits[1 : (T-2)] = Qwaits[1 : (T-2)] + alpha * (stepRwdSignals[1 : (T-2)] - Qwaits[1 : (T-2)])
        }
      }
      
      # update V0
      itiDelta = itiRwdSignal - V0
      V0 = V0 + alpha * itiDelta
      
      # record updated values
      Qwaits_[,tIdx + 1] = Qwaits
      V0_[tIdx + 1] = V0
    }# end of the value update section
  } # end of the loop over trials
  
  # return outputs
  outputs = list( 
    "trialNum" = 1 : nTrial, 
    "condition" = condition_,
    "trialEarnings" = trialEarnings_, 
    "timeWaited" = timeWaited_,
    "sellTime" = sellTime_,
    "scheduledWait" = scheduledWait_,
    "Qwaits_" = Qwaits_, 
    "V0_" = V0_
  )
  return(outputs)
}