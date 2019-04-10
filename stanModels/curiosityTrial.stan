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
  real<lower = 2, upper = 22> tau;
  real<lower = 0.7, upper = 1> gamma;
}
transformed parameters{
  // initialize action values 
  vector[nTimeSteps] Qwait = rep_vector(wIni, nTimeSteps);
  real Qquit = wIni;
  
  // initialize recordings of action values 
  matrix[nTimeSteps, N] Qwaits = rep_matrix(0, nTimeSteps, N);
  vector[N] Qquits = rep_vector(0, N);
  
  
  // initialize trialReward and nextWaitRateHat
  real trialReward;
  real nextWaitRateHat;

  // define gamma List
  // we assume rewards occur at the end of the last step
  // and the action value is for the beginning of the step
  // e.g, wait at the first step and get the reward at the end of the step
  // Qwait = tokenValue * gamma * 1
  // for Qquit, gamma ^ iti /stepDuration
  vector[nTimeSteps] gammaList;
  for(i in 1 : nTimeSteps){
    gammaList[i] = gamma ^ (nTimeSteps - i + 1);
  }
  
  // fill the first trial of Qwaits and Quits
  Qwaits[,1] = Qwait;
  Qquits[1] = Qquit;

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
      trialReward = tokenValue + nextWaitRateHat * Qwait[1] * gamma ^(iti / stepDuration) + (1 - nextWaitRateHat) * Qquit * gamma ^(iti / stepDuration);
      Qwait[1 : thisNActions] = (1 - phi) * Qwait[1 : thisNActions] + phi * trialReward * gammaList[(nTimeSteps - thisNActions + 1):nTimeSteps];
      # counterfactual thinking 
      Qquit = (1 - phi) * Qquit + phi * trialReward * gammaList[nTimeSteps - thisNActions + 1] * gamma ^(iti / stepDuration);
    }else{
      trialReward = nextWaitRateHat * Qwait[1] * gamma ^(iti / stepDuration) + (1 - nextWaitRateHat) * Qquit * gamma ^(iti / stepDuration);
      Qquit =  (1 - phi) * Qquit + phi *  trialReward;
      if(thisNActions > 1){
        Qwait[1 : (thisNActions - 1)] = (1 - phi) * Qwait[1 : (thisNActions - 1)] + phi * trialReward * gammaList[(nTimeSteps - thisNActions+2): nTimeSteps];
      }
    }
    Qwaits[,tIdx+1] = Qwait;
    Qquits[tIdx+1] = Qquit;
  }// end of the loop
}
model {
  phi ~ uniform(0, 0.3);
  tau ~ uniform(2, 50);
  gamma ~ uniform(0.5, 1);
  
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



