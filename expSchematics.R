# this script plots delay distributions and reward rates in two environments

# load experiment parameters 
load("expParas.RData")

# for display purposes, all variables on the continous time scale
# are discretized into 0.1 second time bins
bin = 0.1 # width of a time bin
time = list(
  HP = seq(bin, tMaxs[1], by = bin),
  LP = seq(bin, tMaxs[2], by = bin)
) 

### delays in HP follow unif(0, 20)
### delays in LP follow pareto(mu = 0, k =4, sigma = 2), truncated at 40
# delay CDFS
HP = 1 / length(time[['HP']]) * (1 : length(time[['HP']]))
mu = 0; k = 4; sigma = 2; # parameters for the pareto distribution
LP = 1 - (1 + k * (time[['LP']]- mu) / sigma) ^ (-1 / k)
LP[length(LP)] = 1 # truncated at 40
rewardDelayCDFs = list(
  HP = HP,
  LP = LP
)
# delay PDFs
HP = diff(c(0, rewardDelayCDFs$HP))
LP = diff(c(0, rewardDelayCDFs$LP))
rewardDelayPDFs = list(
  "HP" = HP,
  "LP" = LP
)

# average waiting durations given different policies
# Here we assume rewards occur at the middle of each time bin
HP = cumsum((time[['HP']] - 0.5 * bin) * rewardDelayPDFs$HP) / cumsum(rewardDelayPDFs$HP)
LP = cumsum((time[['LP']] - 0.5 * bin) * rewardDelayPDFs$LP) / cumsum(rewardDelayPDFs$LP)
meanRewardDelays = list('HP' = HP, 'LP' = LP)

# rewardRates given different policies
HP = tokenValue * rewardDelayCDFs$HP /
  ((meanRewardDelays$HP * rewardDelayCDFs$HP + time[['HP']] * (1 - rewardDelayCDFs$HP)) + iti)
LP = tokenValue * rewardDelayCDFs$LP /
  ((meanRewardDelays$LP * rewardDelayCDFs$LP + time[['LP']] * (1 - rewardDelayCDFs$LP)) + iti)
rewardRates = list('HP' = HP, 'LP' = LP)

# optimal raward rates and optimal policies
optimWaitThresholds = list()
optimWaitThresholds$HP = time$HP[which.max(HP)]
optimWaitThresholds$LP = time$LP[which.max(LP)]
optimRewardRates = list()
optimRewardRates$HP = max(HP)
optimRewardRates$LP = max(LP)


# plot CDFs 
library('ggplot2')
source('subFxs/plotThemes.R')
library("tidyr"); library('dplyr')
dir.create('figures/expSchematics')
data.frame(CDF = c(0,rewardDelayCDFs$HP, 0, rewardDelayCDFs$LP),
           time = c(0, time$HP, 0, time$LP),
           condition =  rep(c('HP', 'LP'), c(length(time$HP) + 1, length(time$LP) + 1))) %>%
  ggplot(aes(time, CDF)) + geom_line(size = 3) + facet_grid(~condition) +
  ylim(c(0,1)) + 
  myTheme + xlab('Delay duration (s)') + ylab('CDF')
ggsave('figures/expSchematics/CDF.png', width =6, height = 3)

# plot reward rates
optimData = data.frame(condition = c("HP", "LP"), waitThreshold = as.double(optimWaitThresholds))
data.frame(rewardRate = c(0, rewardRates[[1]], 0, rewardRates[[2]]),
           time = c(0, time[[1]], 0, time[[2]]),
           condition = rep(c("HP", "LP"), c(length(time$HP) + 1, length(time$LP) + 1))) %>%
  ggplot(aes(time, rewardRate)) +
  geom_line(size = 3)  + myTheme + 
  ylab(expression(bold("Reward rate (cent s"^"-1"*")"))) + xlab("Waiting policy (s)")  +
  geom_vline(data = optimData, aes(xintercept = waitThreshold),
             linetype = "dashed", size = 1.5) + facet_grid(~condition)
ggsave("figures/expSchematics/reward_rate.png", width = 6, height = 3)
