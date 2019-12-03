# fit a model for multiple participants in Rstan 

# inputs:
# modelName : the model name
# trialData : a nSubx1 list, each element containing behavioral data for one participant
# config : a list containing the Rstab configuration 
# outputDir: the directory to save parameter estimations

# the config variable contains:
# nChain : number of chains 
# nIter : number of interations on each chain 
# adapt_delta: real number from 0 to 1, and increaseing it forces stan to explore the target distribution in a smaller step
# max_treedepth: maximal depth of the trees that stan evaluates during each iteration
# warningFile : file for saving warnings generated Rstan

modelFitHM = function(modelName){
  # determine parameters 
  paraNameStems = getParaNames(modelName)
    
  # load sub-functions and packages
  library("dplyr"); library("tidyr")
  source("subFxs/loadFxs.R")
  source("subFxs/helpFxs.R")
  library('rstan');library("loo");library("coda") 
  load("expParas.RData")
  stepSec = 1 # duration of one time step (namely one temporal state) 
  nTrialMax =  350 # use the max possible value 
  
  # generate output directories
  dir.create("genData")
  dir.create("genData/expModelFit")
  dir.create("stanWarnings")

  # load expData
  allData = loadAllData()
  hdrData = allData$hdrData
  trialData = allData$trialData
  ids = names(trialData)
  # nSub = length(ids)   
  nSub = 2
  paraNames = c(
    paste0("raw_group_", paraNameStems),
    unlist(lapply(1 : nSub, function(sIdx) paste0(paraNameStems, "s[", sIdx, "]")))
  )
  
  # output and stan config
  outputDir = sprintf("genData/expModelFit/%s", modelName)
  dir.create(outputDir)
  config = list(
    nChain = 4,
    nIter = 100,
    adapt_delta = 0.99,
    max_treedepth = 11,
    warningFile = sprintf("stanWarnings/exp_%s.txt", modelName)
  )
   # create the file for Rstan warnings and erros
  writeLines("", config[['warningFile']])

  # compile the Rstan model 
  options(warn= 1) 
  Sys.setenv(USE_CXX14=1)
  rstan_options(auto_write = TRUE) 
  model = stan_model(file = sprintf("stanModels/%s.stan", modelName))
  
        
  
  # parallel compuation settings
  nCore = as.numeric(Sys.getenv("NSLOTS")) # needed for cluster
  if(is.na(nCore)) nCore = 1 # needed for cluster
  # nCore = parallel::detectCores() -1 
  nChain = config[['nChain']] # number of MCMC chains
  nIter = config[['nIter']] # number of iterations on each chain
  controlList = list(adapt_delta = config[['adapt_delta']],
                     max_treedepth = config[['max_treedepth']] )
  warningFile = config[['warningFile']] # output file for stan warnings and errors
  

  Ns = vector(length = nSub)
  nStepTotals = vector(length = nSub)
  R_ = matrix(0, nTrialMax, nSub)
  T_ = matrix(0, nTrialMax, nSub)
  nWait_s_ = matrix(0, nTrialMax, nSub)
  for(i in 1 : nSub){
    id = ids[[i]]
    thisTrialData = trialData[[id]]
    excludedTrials = which(thisTrialData$trialStartTime > (blockSec - max(tMaxs)))
    thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excludedTrials,]
    thisTrialData = within(thisTrialData, {timeWaited[trialEarnings!= 0] = scheduledWait[trialEarnings!= 0]})
    ## terminal state in each trial
    ## Noticeably, ceiling(timeWaited / stepSec) gives the num of steps
    ## and adding 1 upon it gives the terminal state, namely the state following the final action. 
    N = length(thisTrialData$trialEarnings)
    Ns[i] = N
    Ts = with(thisTrialData, {round(ceiling(timeWaited / stepSec) + 1)}) 
    T_[1 : N, i] = Ts
    Rs = thisTrialData$trialEarnings
    R_[1 : N, i] = Rs
    nWait_s = vector(length = N)
    for(j in 1 : N){
      if(Rs[j] != 0){
        nWait_s[j] = Ts[j] - 1;
      }else{
        nWait_s[j] = Ts[j] - 2;
      }
    }
    nWait_s_[1 : N, i] = nWait_s
  }
  
  inputs <- list(
    iti = iti,
    stepSec = stepSec,
    nStepMax = nStepMax,
    nTrialMax = nTrialMax,
    nWait_s_ = nWait_s_,
    nStepTotal = sum(T_) - sum(Ns),
    nSub = nSub,
    Ns = Ns, # number of trials
    R_ = R_, # rewards on each trial
    T_ = T_)
  
  if(modelName %in% c("QL1", "QL2", "QL1_prime", "QL2_prime", "QL1_HM", "QL2_HM")){
    ## in Q-learning, the initial value of the iti state is proportional to 
    ## the discounted total rewards averaged across two conditions
    ## discount factor for one step is 0.85
    VitiIni = 0.9 * mean(unlist(optimRewardRates) * stepSec / (1 - 0.85))
    inputs$VitiIni = VitiIni
  }else{
    ## in R-learning, the initial reward rate is proportional to
    ## the optimal reward rates averaged across two conditions
    reRateIni = 0.9 * mean(unlist(optimRewardRates)) * stepSec;
    inputs$reRateIni = reRateIni     
  }
  
  # fit the model
  withCallingHandlers({
    fit = sampling(object = model, data = inputs, cores = 1, chains = nChain,
                   iter = nIter, control = controlList) 
    print("Finish!")
    write("Finish", warningFile, append = T, sep = "\n")
  }, warning = function(w){
    warnText = paste(modelName,  w)
    write(warnText, warningFile, append = T, sep = "\n")
  })
  
  samples = fit %>%
    rstan::extract(permuted = F, pars = c("raw_group_phi", "phis[1]"))
  
  # save posterior samples
  samples = samples %>% adply(2, function(x) x) %>% dplyr::select(-chains) 
  fitSummary <- summary(fit, pars = c(paraNames, "LL_all"), use_cache = F)$summary
  
}
