# this script uses a different loss function
# the distance to the original durations.
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
  int nScheduledWaitPoints[N];
  int nTimePoints[N]; // list of available time points 
}
transformed data {
  // constant
  real sigma = 10; // sigma for the normal distribution, actuallt doesn't matter
  real stepDuration = 0.5;
  real iti = 2;
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
  matrix[nTimeStep, N] waitLogProb = rep_matrix(0, nTimeStep, N);
  vector[2] values;
  int I;
  phi ~ uniform(0, 0.3);
  tau ~ uniform(2, 22);
  gamma ~ uniform(0.7, 1);
  QwaitIni ~ uniform(0.5, 9.5);
  // calculate the waitLogProb
  for(tIdx in 1 : N){
    values[2] = Qquits[tIdx] * tau;
    I = 1;
    while(I <= nScheduledWaitPoints[tIdx]){
      values[1] = Qwaits[I, tIdx] * tau;
      waitLogProb[I, tIdx] = categorical_logit_lpmf(1 | values);
      I = I + 1;
    }
  }
  
  //calculate the target log likelyhood 
  for(tIdx in 1 : N){
    target += log(1 - exp(waitLogProb[1,tIdx])) + normal_lpdf(0 | timeWaited[tIdx], sigma);
    I = 1;
    while(I < nScheduledWaitPoints[tIdx]){
      target += log(1 - exp(waitLogProb[I+1,tIdx])) + sum(waitLogProb[1:I, tIdx]) + normal_lpdf(I * stepDuration | timeWaited[tIdx], sigma);
      I = I + 1;
    }
    if(nScheduledWaitPoints[tIdx] > 1){
      target += sum(waitLogProb[1:nScheduledWaitPoints[tIdx], tIdx]) + normal_lpdf(nScheduledWaitPoints[tIdx] * stepDuration | timeWaited[tIdx], sigma);
    }
  }
  
}
generated quantities {
  matrix[nTimeStep + 1, N] timeWaitedLogProb = rep_matrix(0, nTimeStep + 1, N);
  matrix[nTimeStep, N] waitLogProb = rep_matrix(0, nTimeStep, N);
  vector[N] log_lik = rep_vector(0, N);
  real LL_all;
  int I; 
  vector[N] dist = rep_vector(0, N);
  vector[2] values;
  
  // calculate waitLogProb
  for(tIdx in 1 : N){
    values[2] = Qquits[tIdx] * tau;
    I = 1;
    while(I <= nScheduledWaitPoints[tIdx]){
      values[1] = Qwaits[I, tIdx] * tau;
      waitLogProb[I, tIdx] = categorical_logit_lpmf(1 | values);
      I = I + 1;
    }
  }
  
  //calculate timeWaitedLogProb and log_lik
  for(tIdx in 1 : N){
    I = 1;
    timeWaitedLogProb[1, tIdx] =  log(1 - exp(waitLogProb[1,tIdx]));
    log_lik[tIdx] = log_lik[tIdx] + log(1 - exp(waitLogProb[1,tIdx])) + normal_lpdf(0 | timeWaited[tIdx], sigma);
    while(I < nScheduledWaitPoints[tIdx]){
      timeWaitedLogProb[I+1, tIdx] = log(1 - exp(waitLogProb[I+1,tIdx])) + sum(waitLogProb[1:I, tIdx]);
      log_lik[tIdx] = log_lik[tIdx] + log(1 - exp(waitLogProb[I+1,tIdx])) + sum(waitLogProb[1:I, tIdx]) + normal_lpdf(I * stepDuration | timeWaited[tIdx], sigma);
      I = I + 1;
    }
    if(nScheduledWaitPoints[tIdx] > 1){
    timeWaitedLogProb[nScheduledWaitPoints[tIdx]+1, tIdx] = sum(waitLogProb[1:nScheduledWaitPoints[tIdx], tIdx]);
    log_lik[tIdx] =  log_lik[tIdx] + sum(waitLogProb[1:nScheduledWaitPoints[tIdx], tIdx]) + normal_lpdf(nScheduledWaitPoints[tIdx] * stepDuration | timeWaited[tIdx], sigma);
    }
  }
  
  // calculate distance and log_lik
  for(tIdx in 1 : N){
    for(i in 1 : (nScheduledWaitPoints[tIdx] + 1)){
      dist[tIdx] = dist[tIdx] + exp(timeWaitedLogProb[i,tIdx]) * ( (i - 1) * stepDuration - timeWaited[tIdx]) * ( (i - 1) * stepDuration - timeWaited[tIdx]);
    }
  }
  
  //
  LL_all = sum(log_lik);
}
