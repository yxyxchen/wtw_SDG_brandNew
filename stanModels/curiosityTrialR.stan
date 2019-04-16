data {
  // depending on the condition
  real wIni;
  int tMax;
  int nTimeSteps; // nTimeSteps = tMax / stepDuration
  
  // depending on each subject
  int N; // number of trials
  vector[N] timeWaited;
  vector[N] trialEarnings;
  int nActionsPerTrial[N]; // total actions participants made for each trial
}
transformed data {
  // constant
  real stepDuration = 0.5;
  real iti = 2;
  real tokenValue = 10;
  int totalSteps = sum(nActionsPerTrial);
  }
parameters {
  real<lower = 0, upper = 0.3> phi;
  real<lower = 2, upper = 50> tau;
  real<lower = 0, upper = 0.3> phiR;
}
transformed parameters{
  // initialize action values 
  vector[nTimeSteps] Qwait = rep_vector(0, nTimeSteps);
  real Qquit = 0;
  real Rrate = wIni;
  
  
  // initialize recordings of action values 
  matrix[nTimeSteps, N] Qwaits = rep_matrix(0, nTimeSteps, N);
  vector[N] Qquits = rep_vector(0, N);
  vector[N] Rrates = rep_vector(0, N);
  
  // initialize trialReward and nextWaitRateHat
  real trialReward;
  real nextWaitRateHat;
  
  // define toRewardTimeList
  // we use the 
  vector[nTimeSteps] toRewardTimeList;
  for(i in 1 : nTimeSteps){
    toRewardTimeList[i] = nTimeSteps - i + 1;
  }
  // fill the first trial of Qwaits and Quits
  Qwaits[,1] = Qwait;
  Qquits[1] = Qquit;
  Rrates[1] = Rrate;
  
  // initialize opportunity time length
  
  //loop over trial
  for(tIdx in 1 : (N -1)){
    // determine thisNActions
    int thisNActions = nActionsPerTrial[tIdx];
    // determine curiosity
    // tIdx + 1 - 1 = tIdx
    real nextCuriosity =  2 * exp(-0.2 * tIdx);
    // update and track action values
    nextWaitRateHat =  1 / (1  + exp((Qquit - Qwait[1] - nextCuriosity)* tau));
    if(trialEarnings[tIdx] > 0){ 
      trialReward = tokenValue + nextWaitRateHat * Qwait[1]  + (1 - nextWaitRateHat) * Qquit - Rrate * (iti / stepDuration);
      // given gammaList[nTimeSteps] = 1
      Qwait[1 : thisNActions] = (1 - phi) * Qwait[1 : thisNActions] + phi * (trialReward - Rrate * toRewardTimeList[(nTimeSteps - thisNActions + 1): nTimeSteps]);  
      Rrate = (1 - phiR) * Rrate + phiR * (trialReward - Rrate * (toRewardTimeList[(nTimeSteps - thisNActions + 1)] + iti/stepDuration));
      # counterfactual thinking 
      Qquit = (1 - phi) * Qquit + phi * (trialReward - Rrate * (toRewardTimeList[(nTimeSteps - thisNActions + 1)] + iti/stepDuration));
    }else{
      trialReward = nextWaitRateHat * Qwait[1] + (1 - nextWaitRateHat) * Qquit - Rrate * (iti / stepDuration);
      Qquit =  (1 - phi) * Qquit + phi *  trialReward;
      Rrate = (1 - phiR) * Rrate + phiR * (trialReward - Rrate * (toRewardTimeList[(nTimeSteps - thisNActions + 2)] + iti/stepDuration));
      if(thisNActions > 1){
        Qwait[1 : (thisNActions - 1)] = (1 - phi) * Qwait[1 : (thisNActions - 1)] + phi * (trialReward - Rrate * toRewardTimeList[(nTimeSteps - thisNActions + 2): nTimeSteps]);
      }
      # counterfactual thinking 
      if(thisNActions > 1){
        Qquit = (1 - phi) * Qquit + phi * (trialReward - Rrate * (toRewardTimeList[nTimeSteps - thisNActions + 2] + iti/stepDuration));
      }else{
        Qquit = (1 - phi) * Qquit + phi * (trialReward - Rrate * (iti/stepDuration));
      }
    }
    Qwaits[,tIdx+1] = Qwait;
    Qquits[tIdx+1] = Qquit;
    Rrates[tIdx+1] = Rrate;
  }// end of the loop
}
model {
  phi ~ uniform(0, 0.3);
  tau ~ uniform(2, 50);
  phiR ~ uniform(0, 0.3);
  
  // calculate the likelihood 
  for(tIdx in 1 : N){
    int action;
    vector[2] values;
    real curiosity = 2 * exp(-0.2 * (tIdx - 1));
    for(i in 1 : nActionsPerTrial[tIdx]){
    if(trialEarnings[tIdx] == 0 && i == nActionsPerTrial[tIdx]){
      action = 2; // quit
    }else{
      action = 1; // wait
    }
      values[1] = (Qwaits[i, tIdx] + curiosity) * tau;
      values[2] = Qquits[tIdx] * tau;
      //action ~ categorical_logit(values);
      target += categorical_logit_lpmf(action | values);
    } 
  }
}
generated quantities {
// initialize log_lik
  vector[totalSteps] log_lik = rep_vector(0, totalSteps);
  vector[2] values;
  real LL_all;
  int no = 1;
  // loop over trials
  for(tIdx in 1 : N){
    int action;
    real curiosity = 2 * exp(-0.2 * (tIdx - 1));
    for(i in 1 : nActionsPerTrial[tIdx]){
      if(trialEarnings[tIdx] == 0 && i == nActionsPerTrial[tIdx]){
        action = 2; // quit
      }else{
        action = 1; // wait
      }
      values[1] = (Qwaits[i, tIdx] + curiosity) * tau;
      values[2] = Qquits[tIdx] * tau;
      log_lik[no] =categorical_logit_lpmf(action | values);
      no = no + 1;
    }
  }// end of the loop
  LL_all =sum(log_lik);
}



