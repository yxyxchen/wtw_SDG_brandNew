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
  real<lower = 2, upper = 22> tau;
  real<lower = 0.4, upper = 5> k;
  real<lower = 0, upper = nTimeSteps> zeroPoint; 
}
transformed parameters{
  // initialize action values 
  // especially for this version use 0.9, original 1
  real Qquit = wIni;
  real Viti = wIni;
  vector[nTimeSteps] Qwait;
    // initialize variables to record action values 
  matrix[nTimeSteps, N] Qwaits = rep_matrix(0, nTimeSteps, N);
  vector[N] Qquits = rep_vector(0, N);
  vector[N] Vitis = rep_vector(0, N);

  // initialize caching variables
  real G1;
  // fill values
  for(i in 1 : nTimeSteps){
    Qwait[i] = zeroPoint*0.1 - 0.1*(i - 1) + Qquit;
  }
  
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
    if(RT > 0){
      for(t in 1 : (T - 1)){
        real G = RT / (1 + k * (T- t - 1)) + Viti / (1 + k * (T - t));
        Qwait[t] = Qwait[t] + phi * (G - Qwait[t]);
      }
    }else{
      real G =  RT  + Viti / (1 + k);
      Qquit = Qquit + phiP * (G - Qquit) * (G > Qquit) + phiP *  (G - Qquit) * (G <= Qquit) ;
      if(T > 2){
        for(t in 1 : (T-2)){
          G =  RT  / (1 + k *(T - t -1)) + Viti / (1 + k *(T - t));
          Qwait[t] = Qwait[t] + phiP * (G - Qwait[t]);    
        }
      }
    }
    // update Qquit by counterfactual thiking
    G1 =  RT  / ( 1 + k *(T - 2)) + Viti / (1 + k * (T - 1));
    if(RT > 0){
      Qquit = Qquit + phi * (RT  / ( 1 + k *(T - 1 + iti / stepDuration)) + Viti / (1 + k * (T + iti / stepDuration)) - Qquit);
    }else{
      Qquit = Qquit + phiP * (RT  / ( 1 + k *(T - 1 + iti / stepDuration)) + Viti / (1 + k * (T + iti / stepDuration)) - Qquit);
    }
    
    // update Viti
    if(RT > 0){
       Viti = Viti + phi * (RT  / ( 1 + k *(T - 2 + iti / stepDuration)) + Viti / (1 + k * (T - 1 + iti / stepDuration)) - Viti);
    }else{
       Viti = Viti + phiP * (RT  / ( 1 + k *(T - 2 + iti / stepDuration)) + Viti / (1 + k * (T - 1 + iti / stepDuration)) - Viti);
    }
   
    
    // save action values
    Qwaits[,tIdx+1] = Qwait;
    Qquits[tIdx+1] = Qquit;
    Vitis[tIdx + 1] = Viti;
  }// end of the loop
}
model {
  phi ~ uniform(0, 0.3);
  phiP ~ uniform(0, 0.3);
  tau ~ uniform(2, 22);
  k ~ uniform(0.4, 5);
  zeroPoint ~ uniform(0, nTimeSteps);
  
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
      values[1] = Qwaits[i, tIdx] * tau;
      values[2] = Qquits[tIdx] * tau;
      log_lik[no] =categorical_logit_lpmf(action | values);
      no = no + 1;
    }
  }// end of the loop
  LL_all =sum(log_lik);
}
