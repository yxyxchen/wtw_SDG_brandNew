functions {
    matrix Qwait_rng(vector paras, vector trialEarnings, real wIni, int nTimeSteps, int N, int[] Ts) {
      real phi = paras[1];
      real tau = paras[2];
      real gamma = paras[3];
      real zeroPoint = paras[4];
      
      int stepDuration = 1;
      real Qquit = wIni * 0.9;
      real Viti = wIni * 0.9;
      vector[nTimeSteps] Qwait;
      real iti = 2;
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
        if(trialEarnings[tIdx] > 0){
          for(t in 1 : (T - 1)){
            real G = RT * gamma^(T - t -1) + Viti * gamma^(T - t);
            Qwait[t] = Qwait[t] + phi * (G - Qwait[t]);
            if(is_nan(Qwait[t])) print(T-t-1);
          }
        }else{
          real G =  RT  + Viti * gamma;
          Qquit = Qquit + phi * (G - Qquit);
          if(T > 2){
            for(t in 1 : (T-2)){
              G =  RT  * gamma^(T - t -1) + Viti * gamma^(T - t);
              Qwait[t] = Qwait[t] + phi * (G - Qwait[t]);      
            }
          }
        }
            // update Qquit by counterfactual thiking
    G1 =  RT  * gamma^(T - 2) + Viti * gamma^(T - 1);
    Qquit = Qquit + phi * (G1 * gamma^(iti / stepDuration + 1) - Qquit);
        // update Viti
        Viti = Viti + phi * (G1 * gamma^(iti / stepDuration) - Viti);
        
        // save action values
        Qwaits[,tIdx+1] = Qwait;
        Qquits[tIdx+1] = Qquit;
        Vitis[tIdx + 1] = Viti;
    }
    return Qwaits;
  }
}