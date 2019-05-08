data {
  // depending on the condition
  real wIni;
  int tMax;
  int nTimeSteps; // nTimeSteps = tMax / stepDuration
  
  // depending on each subject
  int N; // number of trials
  real timeWaited[N];
  vector[N] trialEarnings;
  int Ts[N]; // terminal time step index 
}
transformed data {
  // constant
  real stepDuration = 1;
  real iti = 2;
  real tokenValue = 10;
  int totalSteps = sum(Ts) - N;
  real QwaitIni = 1;
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
  real Viti = 0;
  real Rrate = wIni;
  
  // initialize variables to record action values 
  matrix[nTimeSteps, N] Qwaits = rep_matrix(QwaitIni, nTimeSteps, N);
  vector[N] Rrates = rep_vector(0, N);
  vector[N] Qquits = rep_vector(0, N);
  vector[N] Vitis = rep_vector(0, N);
  
  // initialize caching variables
  real G1;
  real delta;
  
  // fill the first element of Qwaits, Quits and Vitis 
  Qwaits[,1] = Qwait;
  Qquits[1] = Qquit;
  Vitis[1] = Viti;
 
  //loop over trial
  for(tIdx in 1 : (N -1)){
    // determine the termial timestep T 
    int T = Ts[tIdx];
    real RT = trialEarnings[tIdx];
    
    // update action values for rewarded trials
    if(trialEarnings[tIdx] > 0){
      for(t in 1 : (T - 1)){
        real G = RT - Rrate * (T - t) + Viti;
        Qwait[t] = Qwait[t] + phi * (G - Qwait[t]);
      }
    }else{
      real G =  RT - Rrate + Viti;
      Qquit = Qquit + phi * (G - Qquit);
      if(T > 2){
        for(t in 1 : (T-2)){
          G =  RT - Rrate * (T - t) + Viti;
          Qwait[t] = Qwait[t] + phi * (G - Qwait[t]);          
        }
      }
    }
    // update Qquit by counterfactual thiking
    G1 =  RT - Rrate * (T - 1) + Viti;
    Qquit = Qquit + phi * (G1 - Rrate * (iti /stepDuration + 1) - Qquit);

    // update Viti and Rrate
    delta = (G1 - Rrate * (iti /stepDuration) - Viti);
    Viti = Viti + phi * delta;
    Rrate = Rrate + phiR * delta;
    
    // save action values
    Qwaits[,tIdx+1] = Qwait;
    Qquits[tIdx+1] = Qquit;
    Vitis[tIdx + 1] = Viti;
    Rrates[tIdx + 1] = Rrate;
    
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
    int T = Ts[tIdx];
    for(i in 1 : (T - 1)){
      if(trialEarnings[tIdx] == 0 && i == (T-1)){
        action = 2; // quit
      }else{
        action = 1; // wait
      }
      values[1] = (Qwaits[i, tIdx]) * tau;
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
    int T = Ts[tIdx];
    for(i in 1 : (T - 1)){
      if(trialEarnings[tIdx] == 0 && i == (T-1)){
        action = 2; // quit
      }else{
        action = 1; // wait
      }
      values[1] = (Qwaits[i, tIdx]) * tau;
      values[2] = Qquits[tIdx] * tau;
      log_lik[no] =categorical_logit_lpmf(action | values);
      no = no + 1;
    }
  }// end of the loop
  LL_all =sum(log_lik);
}

