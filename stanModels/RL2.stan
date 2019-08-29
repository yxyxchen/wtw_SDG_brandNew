data {
  // experiment parameters
  real stepSec;// duration of one step
  int nStepMax; // max num of steps in one trial
  real iti;// iti duration, unit = second
  
  // initialization parameters
  real reRateIni; 
  
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
  // phi_pos : action values learning rate for postive outcomes
  // phi_neg : action values leanring rate for negative outcomes
  // beta : reward rate learning rate
  // tau : action consistency, namely soft-max temperature 
  // prior : prior belief parameter
  
  // Additionally, since phi_pos and phi_neg are correlated
  // we reparameterize phi_neg as its ratio relative to phi_pos
  // since phi_pos, phi_neg, and beta are correlated
  // we reparameterize phi_nega and beta as ratios relative to phi_pos
  real<lower = -0.5, upper = 0.5> raw_phi_pos; 
  real<lower = -0.5, upper = 0.5> raw_phi_neg_ratio; 
  real<lower = -0.5, upper = 0.5> raw_beta_ratio; 
  real<lower = -0.5, upper = 0.5> raw_tau;
  real<lower = -0.5, upper = 0.5> raw_prior;
}
transformed parameters{
  // scale raw parameters 
  real phi_pos = (raw_phi_pos + 0.5) * 0.3; // phi_pos ~ unif(0, 0.3)
  real phi_neg_ratio = (raw_phi_neg_ratio + 0.5) * 5; // np_ratio ~ unif(0, 5)
  real beta_ratio =  (raw_beta_ratio + 0.5) * 1; // bp_ratio ~ unif(0, 1)
  real phi_neg = min([phi_pos * phi_neg_ratio, 1]');// convert phi_neg_ratio to phi_neg
  real beta = phi_pos * beta_ratio; // convert beta_ratio to beta
  real tau = (raw_tau + 0.5) * 21.9 + 0.1; // tau ~ unif(0.1, 22)
  real prior = (raw_prior + 0.5) * 65; // prior ~ unif(0, 65)
  
  // declare variables
  real Viti; // value of the ITI state
  real reRate; // reward rate
  vector[nStepMax] Qwaits; // action value of waiting in each step after the ITI state
  vector[N] Viti_ = rep_vector(0, N); // recording Viti
  vector[N] reRate_ = rep_vector(0, N);// recording reRate
  matrix[nStepMax, N] Qwaits_ = rep_matrix(0, nStepMax, N); // recording Qwaits
  real delta; // prediction error to update Viti and reRate
  // initialize variables
  Viti = 0;
  reRate = reRateIni;
  // the initial waiting value delines with elapsed time 
  // and the prior parameter determines at which step it falls below Viti
  for(i in 1 : nStepMax){
    Qwaits[i] = (prior - i) * 0.1 + Viti;
  }
  
  // record initial action values
  Viti_[1] = Viti;
  reRate_[1] = reRate;
  Qwaits_[,1] = Qwaits;
  
  //loop over trials
  for(tIdx in 1 : (N -1)){
    int T = Ts[tIdx]; // current terminal state
    int R = Rs[tIdx]; // current reward
    
    // update action values for rewarded trials
    if(R > 0){
      for(t in 1 : (T - 1)){
        // calculate the expected return G 
        real G = R - reRate * (T - t) + Viti;
        Qwaits[t] = Qwaits[t] + phi_pos * (G - Qwaits[t]);
      }
    }else{
      // if R <=0, Qwait in the last step will not be updated
      // since the agent proceed to the next ITI on that step
      if(T > 2){
        for(t in 1 : (T-2)){
          real G =  R  - reRate * (T - t) + Viti;
          Qwaits[t] = Qwaits[t] + phi_neg * (G - Qwaits[t]);    
        }
      }
    }

    // update Viti and reRate
    delta =  (R - reRate * (T - 1 + iti / stepSec) + Viti - Viti);
    if(R > 0){
      Viti = Viti + phi_pos * delta;
    }else{
      Viti = Viti + phi_neg * delta;
    }
    reRate = reRate + beta * delta;

    
    // save action values
    Qwaits_[,tIdx+1] = Qwaits;
    Viti_[tIdx+1] = Viti;
    reRate_[tIdx + 1] = reRate;
  }
}
model {
  // delcare variables 
  int action; 
  vector[2] actionValues; 
  
  // distributions for raw parameters
  raw_phi_pos ~ uniform(-0.5, 0.5);
  raw_phi_neg_ratio ~ uniform(-0.5, 0.5);
  raw_beta_ratio ~ uniform(-0.5, 0.5);
  raw_tau ~ uniform(-0.5, 0.5);
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
