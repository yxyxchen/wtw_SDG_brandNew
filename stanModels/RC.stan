data {
  // experiment parameters
  real stepSec;// duration between two decision points
  real iti;// iti duration
  int  nDecPoint; // number of possible decision points
  real tWaits[nDecPoint]; // time for each decision point 
  
  // initial value for V0
  real V0_ini; 
  
  // empirical data
  int N; // number of trials
  int Rs[N]; // payoff in each trial
  real Ts[N]; // a trial ends at t == T
  int lastDecPoints[N];// at which decision point each trial ends
}
transformed data {
  // total number of decision points in all trials
  int nDecPointTotal = sum(lastDecPoints);
}
parameters {
  // parameters:
  // kappa: the prob of waiting 
  
  // for computational efficiency,we sample raw parameters from unif(-0.5, 0.5)
  // which are later transformed into actual parameters
  real<lower = -0.5, upper = 0.5> raw_kappa;

}
transformed parameters{
  // transfer paras
  real kappa = (raw_kappa + 0.5) ; // kappa ~ unif(0, 1)
}
model {
  // delcare variables 
  int action; 
  
  // sample
  raw_kappa ~ uniform(-0.5, 0.5);
  
  // loop over trials
  for(tIdx in 1 : N){
    real T = Ts[tIdx]; // this trial ends on t = T
    int R = Rs[tIdx]; // payoff in this trial
    int lastDecPoint = lastDecPoints[tIdx]; // last decision point in this trial
    
    // loop over decision points
    for(i in 1 : lastDecPoint){
      // the agent wait in every decision point in rewarded trials
      // and wait except for the last decision point in non-rewarded trials
      if(R == 0 && i == lastDecPoint){
        action = 2; // quit
      }else{
        action = 1; // wait
      }
      // calculate the likelihood 
      target +=  bernoulli_lpmf(action | kappa);
    } 
  }
}
generated quantities {
 // initialize variables
  int action;
  vector[nDecPointTotal] log_lik = rep_vector(0, nDecPointTotal); // trial-wise log likelihood 
  real totalLL; // total log likelihood
  int no = 1; // action index
  
  // loop over trials
  for(tIdx in 1 : N){
    real T = Ts[tIdx]; // this trial ends on t = T
    int R = Rs[tIdx]; // payoff in this trial
    int lastDecPoint = lastDecPoints[tIdx]; // last decision point in this trial
    
    // loop over decision points
    for(i in 1 : lastDecPoint){
      if(R == 0 && i == lastDecPoint){
        action = 2; // quit
      }else{
        action = 1; // wait
      }
      // calculate the likelihood 
      log_lik[no] = bernoulli_lpmf(action | kappa);
      no = no + 1;
    }
  }
  // calculate total log likelihood
  totalLL =sum(log_lik);
}



