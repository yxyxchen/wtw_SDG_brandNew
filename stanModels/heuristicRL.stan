data {
  // depending on the condition
  real wIni;
  int tMax;
  int nTimeSteps; // nTimeSteps = tMax / stepDuration
  
  // depending on each subject
  int N; // number of trials
  vector[N] trialEarnings;
  real timeWaited[N];
  int Ts[N]; 
  real ini;
}
parameters {
  real<lower = 0, upper = 1> phi;
  real<lower = 0.01, upper = 15> sigma;
  //real<lower = 0, upper = 0.3> phiR;
}
transformed parameters{
  // initialize updated values
  real expectedWait = ini;
  
  // initialize variables to record updated values 
  vector[N] expectedWaits = rep_vector(0, N);
  expectedWaits[1] = expectedWait;
  
  //loop over trial
  for(tIdx in 1 : (N -1)){
    expectedWait = expectedWait + phi * (timeWaited[tIdx] - expectedWait);
    expectedWaits[tIdx + 1] = expectedWait;
  }
}
model {
  phi ~ uniform(0, 1);
  sigma ~ uniform(0.01, 15);
  
  // calculate the likelihood 
  for(tIdx in 1 : N){
    if(trialEarnings[tIdx] < 0){
      target += normal_lpdf(timeWaited[tIdx] | expectedWaits[tIdx], sigma);
    }else{
      target += normal_lcdf(2*expectedWaits[tIdx] - timeWaited[tIdx] | expectedWaits[tIdx], sigma);
    }
      
  }
}
generated quantities {
// initialize log_lik
  vector[N] log_lik = rep_vector(0, N);
  real LL_all;
  // loop over trials
  for(tIdx in 1 : N){
    if(trialEarnings[tIdx] < 0){
      log_lik[tIdx] = normal_lpdf(timeWaited[tIdx] | expectedWaits[tIdx], sigma);
    }else{
      log_lik[tIdx] = normal_lcdf(2*expectedWaits[tIdx] - timeWaited[tIdx] | expectedWaits[tIdx], sigma);
    }
  }// end of the loop
  LL_all =sum(log_lik);
}


