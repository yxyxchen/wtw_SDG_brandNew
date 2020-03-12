# this script plots delay distributions and reward rates in two environments

# load libararies
library("ggplot2")
library("tidyr")
library("dplyr")

# load experiment parameters 
load("expParas.RData")

# for display purposes, all variables on the continous time scale
# are discretized into 0.1 second time bins
bin = 0.1 # width of a time bin
time = list(
  HP = seq(bin, tMaxs[1], by = bin),
  LP = seq(bin, tMaxs[2], by = bin)
) 

# delay CDFS
### delays in HP follow unif(0, 20)
### delays in LP follow pareto(mu = 0, k =4, sigma = 2), truncated at 40
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
## might be different from the values used in expParas.R, 
## which are calcuated with a higher temporal resoluation
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
## here we extend the HP CDF to 32s for display purposes
data.frame(CDF = c(0,c(rewardDelayCDFs$HP, rep(1, length(time$LP) - length(time$HP))), 0, rewardDelayCDFs$LP),
           time = c(0, time$LP, 0, time$LP),
           condition =  rep(c('HP', 'LP'), c(length(time$LP) + 1, length(time$LP) + 1))) %>%
  ggplot(aes(time, CDF)) + geom_line(size = 3, color = themeColor) + facet_grid(~condition) +
  ylim(c(0,1)) + scale_y_continuous(breaks = c(0,0.5,1)) + 
  scale_x_continuous(breaks = c(0, max(tMaxs)/ 2, max(tMaxs)), labels = c("0", max(tMaxs)/2, max(tMaxs)),
                     limits = c(0, max(tMaxs) * 1.1)) + 
  myTheme + xlab('Delay duration (s)') + ylab('CDF') + ggtitle(expName) + 
  theme(plot.title = element_text(hjust = 0.5, color = themeColor)) + 
  annotate("text", x = max(tMaxs)/ 2, y = 0.4, label = "+10¢", size = 6)
ggsave('figures/expSchematics/CDF.eps', width =4, height = 3)
ggsave('figures/expSchematics/CDF.png', width =4, height = 3)

# plot reward rates
optimData = data.frame(condition = c("HP", "LP"), waitThreshold = as.double(optimWaitThresholds))
data.frame(rewardRate = c(0, rewardRates[[1]], 0, rewardRates[[2]]),
           time = c(0, time[[1]], 0, time[[2]]),
           condition = rep(c("HP", "LP"), c(length(time$HP) + 1, length(time$LP) + 1))) %>%
  ggplot(aes(time, rewardRate)) +
  geom_line(size = 3, color = themeColor)  + myTheme + 
  ylab(expression(bold("Reward rate (¢ s"^"-1"*")"))) + xlab("Waiting policy (s)") +
  facet_grid(~condition) + ggtitle(expName) + 
  theme(plot.title = element_text(hjust = 0.5, color = themeColor)) + 
  scale_y_continuous(breaks = c(0, 0.4, 0.8, 1.2), limits = c(0, 1.3)) +
  scale_x_continuous(breaks = c(0, max(tMaxs)/ 2, max(tMaxs)),
                    limits = c(0, max(tMaxs) * 1.1)) 
ggsave("figures/expSchematics/reward_rate.eps", width = 4, height = 3)
ggsave("figures/expSchematics/reward_rate.png", width = 4, height = 3)

# plot prior belief
eta = 20
gamma = 0.6
V0 = (optimRewardRates$HP + optimRewardRates$LP) * 0.5 / (1/6) 
Qquit = V0 * gamma
ts = seq(0, max(tMaxs), by = 0.5)
Qwaits = (eta - ts) * 0.1 + Qquit
data = data.frame(x = ts, y = Qwaits, Qquit = rep(Qquit, length = length(ts))) 
data$z = ifelse(data$y > Qquit, data$y, Qquit)
  
ggplot(data, aes(x, y)) +
  geom_line(size = 2, color = "#2166ac") + geom_line(aes(x, Qquit), size = 2, color = "#b2182b") +
  geom_segment(aes(x = eta, xend = eta, y = 0, yend = Qquit), linetype = "dashed") +
  scale_x_continuous(breaks = c(20), labels = c(expression(0.1~eta)), limits = c(0, max(tMaxs))) +
  scale_y_continuous(breaks = NULL) + xlab("t") + ylab("Initial Value") +
  myTheme +
  theme(text = element_text(size=20))
ggsave('figures/expSchematics/prior.eps', width = 3, height = 3)
ggsave('figures/expSchematics/prior.png', width = 3, height = 3)

  


