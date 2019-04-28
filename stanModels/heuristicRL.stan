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
  real<lower = 2, upper = 50> tau;
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
  tau ~ uniform(2, 50);
  
  // calculate the likelihood 
  for(tIdx in 1 : N){
    if(trialEarnings[tIdx] < 0){
      target += normal_lpdf(timeWaited[tIdx] | expectedWaits[tIdx], expectedWaits[tIdx] / tau);
      print(normal_lpdf(timeWaited[tIdx] | expectedWaits[tIdx], expectedWaits[tIdx] / tau));
    }else{
      target += normal_lcdf(2*expectedWaits[tIdx] - timeWaited[tIdx] | expectedWaits[tIdx], expectedWaits[tIdx] / tau);
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
      log_lik[tIdx] = normal_lpdf(timeWaited[tIdx] | expectedWaits[tIdx], expectedWaits[tIdx] / tau);
    }else{
      log_lik[tIdx] = 1 - normal_lcdf(timeWaited[tIdx] | expectedWaits[tIdx], expectedWaits[tIdx] / tau);
    }
  }// end of the loop
  LL_all =sum(log_lik);
}


