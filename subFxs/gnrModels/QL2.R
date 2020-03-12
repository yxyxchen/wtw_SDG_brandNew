# our reinfocement learning generative models simulate adapative persistence behavior as wait-or-quit choices 
# QL1: Q-learning model with a single learning rate
# QL2: Q-learning model with separate learning rates for rewards and non-rewards
# RL1: R-learning model with a single learning rate 
# RL2: R-learning model with separate learning rates for rewards and non-rewards

# inputs:
# paras: parameter values 
# condition_: [nTrialx1 factor] HP or LP
# scheduledWait_ : [nTrialx1 num] delay for each trial

# outputs
# trialNum : [nTrialx1 int] 1 : nTrial
# condition : [nTrialx1 factor] from inputs
# scheduledWait : [nTrialx1 num] from inputs 
# trialEarnings : [nTrialx1 int] payment for each trial, either 10 or 0
# timeWaited : [nTrialx1 num] how long the agent waits after the iti in each trial 
# sellTime : [nTrialx1 num]  when the agent sells the token on each trial 
# Qwaits_ : [20/40 x nTrial num] value of waiting at each second in each trial
# V0_ : [nTrialx1 num] value of entering a pre-trial iti, namely t = 0

QL2 = function(paras, condition_, scheduledWait_){
  # extract learning parameter 
  alphaR = paras[1]; alphaU = paras[2]; tau = paras[3]; gamma = paras[4]; prior = paras[5]
  
  # num of trials
  nTrial = length(scheduledWait_) # num of trials 
  
  # parameters for discrete temproal states 
  stepSec = 1 # duration of a sampling interval 
  tMax = ifelse(unique(condition_) == "HP", tMaxs[1], tMaxs[2])
  
  # initialize action values 
  V0 = mean(unlist(optimRewardRates)) * stepSec / (1 - 0.85) # state value for t = 0
  tWaits = seq(2, tMax + iti - stepSec, by = stepSec) # decision points 
  Qwaits = -0.1 * (tWaits) + prior + gamma * V0 # value of waiting at each decision points
  
  # initialize output variables
  Qwaits_ = matrix(NA, length(tWaits), nTrial); Qwaits_[,1] = Qwaits 
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
    
    # sample at a temporal resolution of 1 sec until a trial ends
    t = 0
    while(t < (tMax + iti)){
      # take actions after the iti
      if(t >= iti){
        # decide whether to wait or quit
        pWait =  1 / sum(1  + exp((V0 - Qwaits[tWaits == t])* tau))
        action = ifelse(runif(1) < pWait, 'wait', 'quit')
        
        # if a reward occurs and the agent is still waiting, the agent gets the reward
        rewardOccur = ((scheduledWait + iti) >= t) & (scheduledWait + iti) < (t + stepSec)
        getReward = (action == 'wait' && rewardOccur)
        
        # a trial ends if the agent gets a reward or quits. if the trial ends, 
        # return to t = 0. Otherwise, proceed to t + 1.
        isTerminal = (getReward || action == "quit")
        if(isTerminal){
          # update trial-wise variables 
          T = t
          Tidx = which(tWaits == T) # where tWaits == T
          trialEarnings = ifelse(getReward, tokenValue, 0) 
          timeWaited =  ifelse(getReward, scheduledWait, t + stepSec - iti) 
          sellTime = elapsedTime + timeWaited # elapsed task time when the agent sells the token
          elapsedTime = elapsedTime + timeWaited + iti  # elapsed task time before the next token appears
          # record trial-wise variables
          trialEarnings_[tIdx] = trialEarnings
          timeWaited_[tIdx] = timeWaited
          sellTime_[tIdx] = sellTime
          break
        }
      }
      t = t + stepSec
    }
    
    # update value functions at the end of each trial
    if(tIdx < nTrial){
      # determine the learning rate depending on trialEarnings
      if(trialEarnings > 0){
        alpha = alphaR
      }else{
        alpha = alphaU
      }
      
      # calculate expected returns for 2 <= t <= T
      Gts = gamma ^ (T - tWaits[1 : Tidx]) * trialEarnings + gamma ^ (T + 1 - tWaits[1 : Tidx]) * V0
      # update Qwaits
      Qwaits[1 : Tidx] = Qwaits[1 : Tidx] + alpha * (Gts - Qwaits[1 : Tidx])
      
      # calculate expected returns for t == 0 and update V0
      Gt = gamma ^ (T - 0) * trialEarnings + gamma ^ (T + 1 - 0) * V0
      V0 = V0 + alpha * (Gt - V0)
      
      # record updated values
      Qwaits_[,tIdx + 1] = Qwaits
      V0_[tIdx + 1] = V0
    }# end of the loop within a trial 
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