data {
  // experiment parameters
  real stepSec;// duration of one step
  int nStepMax; // max num of steps (actions) in one trial
  real iti;// iti duration, unit = second
  
  // initialization parameters
  real VitiIni; 
  
  // experiment data
  int N; // number of trials
  int Rs[N]; // reward in each trial
  int Ts[N]; // terminal state in each trial
}
transformed data {
  // total number of steps in all trials
  // remmember terminal state - 1 gives the number of steps in a trial
  int nStepTotal = sum(Ts) - N; 
}
parameters {
  // parameters:
  // phi_pos : learning rate for postive outcomes
  // phi_neg: learning rate for negative outcomes or reward omssions 
  // tau : action consistency, namely soft-max temperature paraeter
  // gamma: discount factor
  // prior: prior belief parameter
  
  // for computational efficiency,we sample raw parameters from unif(-0.5, 0.5)
  // which are later transformed into actual parameters
  // Additionally, since phi_pos and phi_neg are correlated
  // we reparameterize phi_neg as its ratio relative to phi_pos
  real<lower = -0.5, upper = 0.5> raw_phi_pos;
  real<lower = -0.5, upper = 0.5> raw_phi_neg_ratio; 
  real<lower = -0.5, upper = 0.5> raw_tau;
  real<lower = -0.5, upper = 0.5> raw_gamma;
  real<lower = -0.5, upper = 0.5> raw_prior;
}
transformed parameters{
  // scale raw parameters 
  real phi_pos = (raw_phi_pos + 0.5) * 0.3; // phi_pos ~ unif(0, 0.3)
  real phi_neg_ratio = (raw_phi_neg_ratio + 0.5) * 5; // phi_neg_ratio ~ unif(0, 5)
  real phi_neg = min([phi_pos * phi_neg_ratio, 1]'); // convert phi_neg_ratio to phi_neg
  real tau = (raw_tau + 0.5) * 21.9 + 0.1; // tau ~ unif(0.1, 22)
  real gamma = (raw_gamma + 0.5) * 0.3 + 0.4; // gamma ~ unif(0.7, 1)
  real prior = (raw_prior + 0.5) * 65; // prior ~ unif(0, 65)
  
  // declare action values 
  // // state value of the ITI state
  real Viti; 
  // // action value of waiting in each step after ITI
  vector[nStepMax] Qwaits; 
  
  // declare variables to record action values 
  matrix[nStepMax, N] Qwaits_ = rep_matrix(0, nStepMax, N);
  vector[N] Viti_ = rep_vector(0, N);
  
  // initialize action values 
  //// the initial value of the ITI state 
  Viti = VitiIni; 
  // the initial waiting value delines with elapsed time 
  // and the prior parameter determines at which step it falls below Viti
  for(i in 1 : nStepMax){
    Qwaits[i] = (prior - i) * 0.1 + Viti;
  }
  
  // record initial action values
  Qwaits_[,1] = Qwaits;
  Viti_[1] = Viti;
 
  //loop over trials
  for(tIdx in 1 : (N - 1)){
    int T = Ts[tIdx]; // current terminal state
    int R = Rs[tIdx]; // current reward
    
    // update action values 
    if(R > 0){
      for(t in 1 : (T - 1)){
        // calculate the expected return G 
        real G = R * gamma^(T - t -1) + Viti * gamma^(T - t);
        // update Qwait towards G
        Qwaits[t] = Qwaits[t] + phi_pos * (G - Qwaits[t]);
      }
    }else{
      // if R <=0, Qwait in the last step will not be updated
      // since the agent proceed to the next ITI on that step
      if(T > 2){
        for(t in 1 : (T-2)){
          real G =  R  * gamma^(T - t -1) + Viti * gamma^(T - t);
          Qwaits[t] = Qwaits[t] + phi_neg * (G - Qwaits[t]);    
        }
      }
    }
    
    // update Viti
    if(R > 0){
      Viti = Viti + phi_pos * ((R  * gamma^(T - 2 + iti / stepSec) + Viti * gamma^(T - 1 + iti / stepSec)) - Viti);
    }else{
      Viti = Viti + phi_neg * ((R  * gamma^(T - 2 + iti / stepSec) + Viti * gamma^(T - 1 + iti / stepSec)) - Viti);
    }
    
    
    // save action values
    Qwaits_[,tIdx+1] = Qwaits;
    Viti_[tIdx+1] = Viti;
  }
}
model {
  // delcare variables 
  int action; 
  vector[2] actionValues; 
  
  // distributions for raw parameters
  raw_phi_pos ~ uniform(-0.5, 0.5);
  raw_phi_neg_ratio ~ uniform(-0.5, 0.5);
  raw_tau ~ uniform(-0.5, 0.5);
  raw_gamma ~ uniform(-0.5, 0.5);
  raw_prior ~ uniform(-0.5, 0.5);
  
  // loop over trials
  for(tIdx in 1 : N){
    int T = Ts[tIdx]; // current terminal state
    int R = Rs[tIdx]; // current reward
    // loop over steps
    for(t in 1 : (T - 1)){
      // determine the action
      // the agent wait in every steps in rewarded trials
      // and wait except for the last step in non-rewarded trials
      if(R == 0 && t == (T-1)){
        action = 2; // proceed to the next ITI
      }else{
        action = 1; // wait
      }
      // calculate the likelihood using the soft-max function
      actionValues[1] = Qwaits_[t, tIdx] * tau;
      actionValues[2] = Viti_[tIdx] * tau;
      target += categorical_logit_lpmf(action | actionValues);
    } 
  }
}
generated quantities {
 // generate action-wise log likelihood and total log likelihood
 
 // initialize variables
  vector[2] actionValues;
  int action;
  vector[nStepTotal] log_lik = rep_vector(0, nStepTotal);
  real LL_all; // total log likelihood
  int no = 1; // action index
  
  // loop over trials
  for(tIdx in 1 : N){
    int T = Ts[tIdx]; // current terminal state
    int R = Rs[tIdx]; // current reward
    // loop over steps
    for(t in 1 : (T - 1)){
      if(R == 0 && t == (T-1)){
        action = 2; // proceed to the next ITI
      }else{
        action = 1; // wait
      }
      // calculate the likelihood using the soft-max function
      actionValues[1] = Qwaits_[t, tIdx] * tau;
      actionValues[2] = Viti_[tIdx] * tau;
      log_lik[no] =categorical_logit_lpmf(action | actionValues);
      no = no + 1;
    }
  }
  // calculate total log likelihood
  LL_all =sum(log_lik);
}
