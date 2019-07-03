data {
  // depending on the condition
  real wIni;
  int tMax;
  int nTimeSteps; // nTimeSteps = tMax / stepDuration
  
  // depending on each subject
  int N; // number of trials
  vector[N] timeWaited;
  vector[N] trialEarnings;
  int Ts[N]; // terminal time step index 
  real stepDuration;
  real iti;
  real tokenValue;
}
transformed data {
  int totalSteps = sum(Ts) - N;
}
parameters {
  real<lower = 0, upper = 0.3> phi;
  real<lower = 0, upper = 0.3> phiP; 
  real<lower = 0.5, upper = 22> tau;
  real<lower = 0, upper = 2> reRateIni; 
  real<lower = 0, upper = 5> slope;
}
transformed parameters{
  // initialize action values 
  // especially for this version use 0.9, original 1
  real reRate = reRateIni;
  vector[N] reRates = rep_vector(0, N);
  reRates[1] = reRate;

  // initialize caching variables
  //loop over trial
  for(tIdx in 1 : (N -1)){
    real RT = trialEarnings[tIdx];
    
    // update action values 
    if(RT > 0){
      real reRateHat = tokenValue / (timeWaited[tIdx] + iti);
      reRate = reRate + phi * (reRateHat - reRate);
    }else{
      real reRateHat = tokenValue / (timeWaited[tIdx] + iti);
      reRate = reRate + phi * (reRateHat - reRate);
    }

    // save action values
    reRates[tIdx+1] = reRate;
  }// end of the loop
}
model {
  phi ~ uniform(0, 0.3);
  phiP ~ uniform(0, 0.3);
  tau ~ uniform(0.5, 22);
  reRateIni ~ uniform(0, 2);
  slope ~ uniform(0, 5);
  
  // calculate the likelihood 
  for(tIdx in 1 : N){
    int action;
    vector[2] values;
    int T = Ts[tIdx];
    for(i in 1 : (T - 1)){
      if(trialEarnings[tIdx] == 0 && i == (T-1)){
        action = 2; // quit
      }else{
        action = 1; // wait
      }
      values[1] = reRates[tIdx] * tau;
      values[2] = tokenValue / (i * slope * stepDuration + iti) * tau;
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
    int T = Ts[tIdx];
    for(i in 1 : (T - 1)){
      if(trialEarnings[tIdx] == 0 && i == (T-1)){
        action = 2; // quit
      }else{
        action = 1; // wait
      }
      values[1] = reRates[tIdx] * tau;
      values[2] = tokenValue / (i * stepDuration * slope + iti) * tau;
      log_lik[no] =categorical_logit_lpmf(action | values);
      no = no + 1;
    }
  }// end of the loop
  LL_all =sum(log_lik);
}
