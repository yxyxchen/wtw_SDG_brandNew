stepSec = 1

## in Q-learning, the initial value of the iti state is proportional to 
## the discounted total rewards averaged across two conditions
## discount factor for one step is 0.85
VitiIni_QL = 0.9 * mean(unlist(optimRewardRates) * stepSec / (1 - 0.85))

## in R-learning, the initial reward rate is proportional to
## the optimal reward rates averaged across two conditions
## and the initial value of the iti state is 0 
reRateIni = 0.9 * mean(unlist(optimRewardRates)) * stepSec
VitiIni_RL = 0
