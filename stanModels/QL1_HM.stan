data {
  // experiment parameters
  real stepSec;// duration of one step
  int nStepMax; // max num of steps in one trial
  real iti;// iti duration
  
  // initialization parameters
  real VitiIni; 
  
  // experiment data
  int nSub;
  int nTrialMax;
  int Ns[nSub]; // number of trials
  int R_[nTrialMax, nSub]; // reward in each trial
  int T_[nTrialMax, nSub]; // terminal state in each trial
  int nWait_s_[nTrialMax, nSub]; // the number of waiting steps in each trial
  int nStepTotal;
  
}
transformed data {
  // median of the individually fitted parameters 
  real phiMu = 0.0220;
  real tauMu = 3.26;
  real gammaMu = 0.885;
  real priorMu = 25.8;
  
  real phiSe = 0.00684;
  real tauSe = 0.478;
  real gammaSe = 0.0124;
  real priorSe = 3.36;
}
parameters {
  // parameters:
  // phi : learning rate
  // tau : action consistency, namely the soft-max temperature parameter
  // gamma: discount factor
  // prior: prior belief parameter
  
  // for computational efficiency,we sample raw parameters from unif(-0.5, 0.5)
  // which are later transformed into actual parameters
  real raw_group_phi;
  real raw_group_tau;
  real raw_group_gamma;
  real raw_group_prior;
  
  real<lower = 0, upper = 1> phis[nSub];
  real<lower = 0, upper = 22> taus[nSub];
  real<lower = 0, upper = 1> gammas[nSub];
  real<lower = 0, upper = 65> priors[nSub];
}
transformed parameters{
  // declare variables 
  // // state value of the ITI state
  real Viti; 
  // // the reward signal for updating action values at the end of each trial 
  real rwdSignal;  
  // // action value of waiting in each step after ITI
  real Qwaits[nStepMax]; 
  // // variables to record action values 
  real Qwaits_[nStepMax, nTrialMax] = rep_array(0.0, nStepMax, nTrialMax);
  real Viti_[nTrialMax] = rep_array(0.0, nTrialMax);
  real QwaitArray[nStepMax, nTrialMax, nSub] = rep_array(0.0, nStepMax, nTrialMax, nSub);
  real VitiArray[nTrialMax, nSub] = rep_array(0.0, nTrialMax, nSub);
  
  for(sIdx in 1 : nSub){
    real phi = phis[sIdx];
    real tau = taus[sIdx];
    real gamma = gammas[sIdx];
    real prior = priors[sIdx];
    int N = Ns[sIdx];
    int nWait_s[N] = nWait_s_[1:N,sIdx];
    int Ts[N] = T_[1 : N,sIdx];
    int Rs[N] = R_[1 : N,sIdx];
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
      //calculate the reward signal for updating action value 
      // which equals R plus the discounted value of the successor state. Noticably, 
      // the successor state at the end of trial is always the iti state before the next trial
      rwdSignal = R + gamma * Viti;
      
      // update Qwaits and Viti towards the discounted reward signals 
      for(t in 1 : nWait_s[tIdx]) {
        real discRwdsignal = rwdSignal * gamma^(T - t);// the discounted reward signal 
        Qwaits[t] = Qwaits[t] + phi * (discRwdsignal - Qwaits[t]);
      }
      Viti = Viti + phi * (gamma ^ (T - 1 + iti / stepSec) * rwdSignal - Viti);
      
      // save action values
      Qwaits_[,tIdx+1] = Qwaits;
      Viti_[tIdx+1] = Viti;
    }
    # save action values
    QwaitArray[:, :, sIdx] = Qwaits_;
    VitiArray[:, sIdx] = Viti_;
  }
}
model {
  // delcare variables 
  int action; 
  vector[2] actionValues; 
  real group_phi;
  real group_tau;
  real group_gamma;
  real group_prior;
  
  // distributions for raw parameters
  // scale raw parameters into real parameters
  raw_group_phi ~ std_normal();
  raw_group_tau ~ std_normal();
  raw_group_gamma ~ std_normal();
  raw_group_prior ~ std_normal();
  
  group_phi = raw_group_phi * phiSe * 2 + phiMu;
  group_tau = raw_group_tau * tauSe * 2 + tauMu; 
  group_gamma = raw_group_gamma * gammaSe * 2 + gammaMu;
  group_prior = raw_group_prior * priorSe * 2 + priorMu;
  
  for(sIdx in 1 : nSub){
    phis[sIdx] ~ normal(group_phi, phiSe * 5) T[0,1];
    taus[sIdx] ~ normal(group_tau, tauSe * 5) T[0,22];
    gammas[sIdx] ~ normal(group_gamma, gammaSe * 5) T[0,1];
    priors[sIdx] ~ normal(group_prior, priorSe * 5) T[0,65];
  }
  
  // loop over participants 
  for(sIdx in 1 : nSub){
    real tau = taus[sIdx];
    real gamma = gammas[sIdx];
    int N = Ns[sIdx];
    int Ts[N] = T_[1 : N,sIdx];
    int Rs[N] = R_[1 : N,sIdx];
    // loop over trials
    for(tIdx in 1 : N){
      int T = Ts[tIdx]; // current terminal state
      int R = Rs[tIdx]; // current reward
      // loop over steps
      action = 1;
      for(t in 1 : (T - 2)){
        // calculate the likelihood using the soft-max function
        actionValues[1] = Qwaits_[t, tIdx] * tau;
        actionValues[2] = Viti_[tIdx] * gamma * tau;
        target += categorical_logit_lpmf(action | actionValues);
      }
      // determine the action
      // the agent wait in every steps in rewarded trials
      // and wait except for the last step in non-rewarded trials
      action = 2 - (R != 0);
      actionValues[1] = Qwaits_[(T-1), tIdx] * tau;
      actionValues[2] = Viti_[tIdx] * gamma * tau;
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
  
  for(sIdx in 1 : nSub){
    real tau = taus[sIdx];
    real gamma = gammas[sIdx];
    int N = Ns[sIdx];
    int Ts[N] = T_[1 : N,sIdx];
    int Rs[N] = R_[1 : N,sIdx];
    // loop over trials
    for(tIdx in 1 : N){
      int T = Ts[tIdx]; // current terminal state
      int R = Rs[tIdx]; // current reward
      // loop over steps
      action = 1;
      for(t in 1 : (T - 2)){
        // calculate the likelihood using the soft-max function
        actionValues[1] = Qwaits_[(T-1), tIdx] * tau;
        actionValues[2] = Viti_[tIdx] * gamma * tau;
        log_lik[no] =categorical_logit_lpmf(action | actionValues);
        no = no + 1;
      }
      action = 2 - (R != 0);
      actionValues[1] = Qwaits_[(T-1), tIdx] * tau;
      actionValues[2] = Viti_[tIdx] * gamma * tau;
      log_lik[no] =categorical_logit_lpmf(action | actionValues);
    }
  }
  // calculate total log likelihood
  LL_all =sum(log_lik);
}
