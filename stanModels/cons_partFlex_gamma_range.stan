data {
  // depending on the condition
  real wInis[2];
  real wIni;
  int tMax;
  int nTimeStep; // since round returns real here, so nTimeStep != tMax / stepDuration
  
  // depending on each subject
  int N; // number of trials
  vector[N] timeWaited;
  vector[N] trialEarnings;
  int nTimePoints[N]; // list of available time points 
  int nScheduledWaitPoints[N];
}
transformed data {
  // constant
  real stepDuration = 0.5;
  real iti = 2;
  real rangeSlope = 0; // how range increase with time, both side included
  int rangeConstant = 2; // in the unit of step, that is 0.5s
  real tokenValue = 10;
  int totalSteps = sum(nTimePoints);
  }
  parameters {
  real<lower = 0, upper = 0.3> phi;
  real<lower = 2, upper = 22> tau;
  real<lower = 0.7, upper = 1> gamma;
  real<lower = 0.5, upper = 9.5> QwaitIni;
}
transformed parameters{
  // initialize action values 
  vector[nTimeStep] Qwait = rep_vector(QwaitIni, nTimeStep);
  real Qquit = QwaitIni * gamma ^ (iti / stepDuration);
  
  // initialize recordings of action values 
  matrix[nTimeStep, N] Qwaits = rep_matrix(0, nTimeStep, N);
  vector[N] Qquits = rep_vector(0, N);
  
  
  // initialize trialReward and nextWaitRateHat
  real trialReward;
  real nextWaitRateHat = 1.0;

  // define gamma List
  vector[nTimeStep] gammaList;
  for(i in 1 : nTimeStep){
    gammaList[i] = gamma ^ (nTimeStep - i);
  }
  
  // fill the first trial of Qwaits and Quits
  Qwaits[,1] = Qwait;
  Qquits[1] = Qquit;

  //loop over trial
  for(tIdx in 1 : (N -1)){
    // determine nTimePoint
    int nTimePoint = nTimePoints[tIdx]; 
    // update and track action values
    if(trialEarnings[tIdx] > 0){
      trialReward = tokenValue;
      Qwait[1 : nTimePoint] = (1 - phi) * Qwait[1 : nTimePoint] + phi * trialReward * gammaList[(nTimeStep - nTimePoint + 1):nTimeStep];
    }else{
      nextWaitRateHat =  1 / (1  + exp((Qquit - Qwait[1])* tau));
      trialReward = nextWaitRateHat * Qwait[1] * gamma ^(iti / stepDuration) + (1 - nextWaitRateHat) * Qquit * gamma ^(iti / stepDuration);
      Qquit =  (1 - phi) * Qquit + phi *  trialReward;
      if(nTimePoint > 1){
        Qwait[1 : (nTimePoint - 1)] = (1 - phi) * Qwait[1 : (nTimePoint - 1)] + phi * trialReward * gammaList[(nTimeStep - nTimePoint + 1):(nTimeStep - 1)];
      }
    }
    Qwaits[,tIdx+1] = Qwait;
    Qquits[tIdx+1] = Qquit;
  }// end of the loop
}
model {
  phi ~ uniform(0, 0.3);
  tau ~ uniform(2, 22);
  gamma ~ uniform(0.7, 1);
  QwaitIni ~ uniform(0.5, 9.5);
  // calculate the likelihood 
  for(tIdx in 1 : N){
    vector[2] values;
    real largerEquelLower;
    real largerUpper;
    // have rewards
    if(trialEarnings[tIdx] > 0){
      // if the range doen't include 0
      if(nTimePoints[tIdx] > rangeConstant){
        for(i in 1 : (nTimePoints[tIdx] - rangeConstant)){
          values[1] = Qwaits[i, tIdx] * tau;
          values[2] = Qquits[tIdx] * tau;
          target += categorical_logit_lpmf(1 | values);
        }
      }else{
        target += 0;
      }
    // have no rewards
    }else{
       // prob of waiting >= the lower bound, lower bound = nTimePoints[tIdx] - rangeConstant - 1
       // e.g. quit at the 4 step, then wait >= 2 step if constant = 1
       largerEquelLower = 0;
        if(nTimePoints[tIdx] > (rangeConstant+1)){
          for(i in 1 : (nTimePoints[tIdx] - rangeConstant - 1)){
            values[1] = Qwaits[i, tIdx] * tau;
            values[2] = Qquits[tIdx] * tau;
            largerEquelLower = largerEquelLower + categorical_logit_lpmf(1 | values);
          }
        }
       //prob of waiting > the upper bound
       // e.g. quit at 4 step(namely wait 3), then wait <= 3+ 1 step, that is < 3 + 2 step
      largerUpper = 0;//??
      if(nScheduledWaitPoints[tIdx] >= (rangeConstant + nTimePoints[tIdx])){
        for(i in 1 : (nTimePoints[tIdx] + rangeConstant)){
          values[1] = Qwaits[i, tIdx] * tau;
          values[2] = Qquits[tIdx] * tau;
          largerUpper = largerUpper + categorical_logit_lpmf(1 | values);
        }
        target += log_diff_exp(largerEquelLower, largerUpper);
      }else{
          target += largerEquelLower;
      }
    }// end of having no rewards
  }// end of trial-wise loop
}
generated quantities {
// initialize log_lik
  vector[N] log_lik = rep_vector(0, N); // it is trial-wise
  real LL_all;
  // calculate the likelihood 
  for(tIdx in 1 : N){
    real largerEquelLower;
    real largerUpper;
    vector[2] values;
    // have rewards
    if(trialEarnings[tIdx] > 0){
      // if the range doen't include 0
      if(nTimePoints[tIdx] > rangeConstant){
        for(i in 1 : (nTimePoints[tIdx] - rangeConstant)){
          values[1] = Qwaits[i, tIdx] * tau;
          values[2] = Qquits[tIdx] * tau;
          log_lik[tIdx] = log_lik[tIdx]  + categorical_logit_lpmf(1 | values);
        }
      }
    // have no rewards
    }else{
       // prob of waiting >= the lower bound, lower bound = nTimePoints[tIdx] - rangeConstant - 1
       // e.g. quit at the 4 step, then wait >= 2 step if constant = 1
       largerEquelLower = 0;
        if(nTimePoints[tIdx] > (rangeConstant+1)){
          for(i in 1 : (nTimePoints[tIdx] - rangeConstant - 1)){
            values[1] = Qwaits[i, tIdx] * tau;
            values[2] = Qquits[tIdx] * tau;
            largerEquelLower = largerEquelLower + categorical_logit_lpmf(1 | values);
          }
        }
       //prob of waiting > the upper bound
       // e.g. quit at 4 step(namely wait 3), then wait <= 3+ 1 step, that is < 3 + 2 step
      largerUpper = 0;
      if(nScheduledWaitPoints[tIdx] >= (rangeConstant+nTimePoints[tIdx])){
        for(i in 1 : (nTimePoints[tIdx] + rangeConstant)){
          values[1] = Qwaits[i, tIdx] * tau;
          values[2] = Qquits[tIdx] * tau;
          largerUpper = largerUpper + categorical_logit_lpmf(1 | values);
        }
        log_lik[tIdx] = log_lik[tIdx] + log_diff_exp(largerEquelLower, largerUpper);
      }else{
          log_lik[tIdx] = largerEquelLower;
      }
    }// end of having no rewards
  }// end of trial-wise loop
  LL_all = sum(log_lik);
}



