data {
  // depending on the task
  real stepDuration;
  real iti;
  
  // depending on the condition
  real wIni;
  int nTimeSteps; // nTimeSteps = tMax / stepDuration
  
  // depending on each subject
  int N; // number of trials
  vector[N] trialEarnings;
  int Ts[N]; // terminal time step index 
}
transformed data {
  int totalSteps = sum(Ts) - N;
}
parameters {
  real<lower = 0, upper = 1> raw_phi;
  real<lower = 0, upper = 1> raw_beta;
  real<lower = 0, upper = 1> raw_nega;
  real<lower = 0, upper = 1> raw_tau;
  real<lower = 0, upper = 1> raw_prior;
}
transformed parameters{
  // transfer paras
  real phi = raw_phi * 0.3;
  real beta = raw_beta * 0.3;
  real nega = min([raw_nega * 5, 1 / beta, 1 / phi]');
  real tau = raw_tau * 21.9 + 0.1;
  real prior = raw_prior * 65;
  // initialize action values 
  real Viti = 0;
  real reRate = wIni;// reward rate for each step
  vector[nTimeSteps] Qwait;
  
  // initialize variables to record action values 
  matrix[nTimeSteps, N] Qwaits = rep_matrix(0, nTimeSteps, N);
  vector[N] Vitis = rep_vector(0, N);
  vector[N] reRates = rep_vector(0, N);
  
  // initialize caching variables
  real G1;
  real delta;
  
  // enter the initial values
  for(i in 1 : nTimeSteps){
    Qwait[i] = prior*0.1 - 0.1*(i - 1) + Viti;
  }
  Qwaits[,1] = Qwait;
  Vitis[1] = Viti;
  reRates[1] = reRate;
  
  //loop over trials
  for(tIdx in 1 : (N -1)){
    // determine the termial timestep T 
    int T = Ts[tIdx];
    real RT = trialEarnings[tIdx];
    
    // update action values for rewarded trials
    if(RT > 0){
      for(t in 1 : (T - 1)){
        real G = RT - reRate * (T - t) + Viti;
        Qwait[t] = Qwait[t] + phi * (G - Qwait[t]);
      }
    }else{
      if(T > 2){
        for(t in 1 : (T-2)){
          real G =  RT  - reRate * (T - t) + Viti;
          Qwait[t] = Qwait[t] + phi * nega * (G - Qwait[t]);    
        }
      }
    }
    
    // update Viti
    G1 =  RT  - reRate*(T - 1) + Viti;
    delta = (G1 - reRate * (iti / stepDuration) - Viti);
    if(RT > 0){
       Viti = Viti + phi * delta;
    }else{
       Viti = Viti + phi * nega * delta;
    }
   
    // update reRate 
    if(RT > 0){
      reRate = reRate +  beta * delta;
    }else{
      reRate = reRate + beta * nega * delta;
    }
    
    // save action values
    Qwaits[,tIdx+1] = Qwait;
    Vitis[tIdx+1] = Viti;
    reRates[tIdx + 1] = reRate;
    
  }// end of the loop
}
model {
  raw_phi ~ uniform(0, 1);
  raw_beta ~ uniform(0, 1);
  raw_nega ~ uniform(0, 1);
  raw_tau ~ uniform(0, 1);
  raw_prior ~ uniform(0, 1);
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
      values[1] = Qwaits[i, tIdx] * tau;
      values[2] = Vitis[tIdx] * tau;
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
      values[1] = Qwaits[i, tIdx] * tau;
      values[2] = Vitis[tIdx] * tau;
      log_lik[no] =categorical_logit_lpmf(action | values);
      no = no + 1;
    }
  }// end of the loop
  LL_all =sum(log_lik);
}
