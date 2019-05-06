library("ggplot2")
load("wtwSettings.RData")
source("subFxs/plotThemes.R")
dir.create("figures")
dir.create("figures/plotFigures")


tMax = 40
k = 4
mu = 0
sigma = 2
dx = 0.1
# prob density of delays 
x = seq(0.1, 40, by = dx)
y1 = c(rep(1 / 20, length(x) / 2), rep(NA, length(x) / 2))
tempt = 1 - (1 + k * (x - mu) / sigma) ^ (-1 / k)
y2 = diff(c(0, tempt)) / dx # be careful when change delta cdf to density, you need to divide it by ...
tempt[length(x)] = 1
y2[length(x)] = y2[length(x)] + 1 - sum(y2 * dx)
plotData = data.frame(time = c(x, x), probability = c(y1, y2), condition = rep(c("HP", "LP"), each = length(x)))
## prob density for HP
ggplot(plotData[plotData$condition == "HP",], aes(time, probability))+ geom_point(color = conditionColors[1])+
  saveTheme + ylab("Probability density") + xlab("Total delay / s") + ggtitle("High Persistence") + ylim(c(0, 0.5)) + xlim(c(0,40))
ggsave(filename = "figures/plotFigures/totalDelay_HP.pdf", width = 6, height = 4)
## prob density for LP
ggplot(plotData[plotData$condition == "LP",], aes(time, probability))+ geom_point(color = conditionColors[2])+
  saveTheme + ylab("Probability density") + xlab("Total delay / s") + ggtitle("Low Persistence") + ylim(c(0, 0.5)) + xlim(c(0,40))
ggsave(filename = "figures/plotFigures/totalDelay_LP.pdf", width = 6, height = 4)

# prob density of remaining delays after 5 seconds
plotData1 = plotData[x > 5,]
plotData1$time = plotData1$time - 5
ggplot(plotData1[plotData1$condition == "HP",], aes(time, probability))+ geom_point(color = conditionColors[1])+
  saveTheme + ylab("Probability density") + xlab("Remaining delay after 5 secs / s") + ggtitle("High Persistence") + ylim(c(0, 0.5)) + xlim(c(0,40))
ggsave(filename = "figures/plotFigures/remianDelay5_HP.pdf", width = 6, height = 4)
ggplot(plotData1[plotData1$condition == "LP",], aes(time, probability))+ geom_point(color = conditionColors[2])+
  saveTheme + ylab("Probability density") + xlab("Remaining delay after 5 secs / s") + ggtitle("Low Persistence") + ylim(c(0, 0.5)) + xlim(c(0,40))
ggsave(filename = "figures/plotFigures/remianDelay5_LP.pdf", width = 6, height = 4)


# prob density of remaining delays after 15 seconds
plotData1 = plotData[x > 15,]
plotData1$time = plotData1$time - 15
ggplot(plotData1[plotData1$condition == "HP",], aes(time, probability))+ geom_point(color = conditionColors[1])+
  saveTheme + ylab("Probability density") + xlab("Remaining delay after 15 secs / s") + ggtitle("High Persistence") + ylim(c(0, 0.4)) + xlim(c(0,40))
ggsave(filename = "figures/plotFigures/remianDelay15_HP.pdf", width = 6, height = 4)
ggplot(plotData1[plotData1$condition == "LP",], aes(time, probability))+ geom_point(color = conditionColors[2])+
  saveTheme + ylab("Probability density") + xlab("Remaining delay after 15 secs / s") + ggtitle("Low Persistence") + ylim(c(0, 0.4)) + xlim(c(0,40))
ggsave(filename = "figures/plotFigures/remianDelay15_LP.pdf", width = 6, height = 4)

# plot expected remain delay as a function of 
x = seq(0.1, 40, by = 0.1);
dx = 0.1
expDelay1 = c((x[1:(length(x) / 2)] + 20) / 2 - x[1:(length(x) / 2)],
              rep(NA, length(x)/2))
expDelay2 = vector(length = length(x))
for(i in 1 : length(expDelay2)){
  expDelay2[i] = sum(y2[i : length(x)] * (x[i:length(x)]-0.5)) / sum(y2[i : length(x)])
}
plotData3 = data.frame(time = rep(x, 2), expDelay = c(expDelay1, expDelay2),
                      condition = rep(c("HP", "LP"), each = length(x)))
ggplot(plotData3[plotData3$condition == "HP",], aes(time, expDelay)) +geom_point(color = conditionColors[1]) +
  xlab("Elapsed time / s") + ylab("Expected remaining delay / s") + ggtitle("High Persistence") + saveTheme
ggsave(filename = "figures/plotFigures/eDelay_HP.pdf", width = 6, height = 4)

ggplot(plotData3[plotData3$condition == "LP",], aes(time, expDelay)) +geom_point(color = conditionColors[2]) +
  xlab("Elapsed time / s") + ylab("Expected remaining delay / s") + ggtitle("Low Persistence") + saveTheme
ggsave(filename = "figures/plotFigures/eDelay_LP.pdf", width = 6, height = 4)



# plot softMax functions
x = log(seq(0.2, 0.8, by= 0.001) / (1 - seq(0.1, 0.9, by= 0.01)))
tauList = c(1, 3, 9)
y = unlist(lapply(tauList, function(tau) 1 / (1 + exp(-x * tau))))
plotData = data.frame(x = rep(x, length(tauList)), y = y, tau = rep(tauList, each = length(x)))
plotData$tau = as.factor(plotData$tau)
ggplot(plotData, aes(x, y, color = tau)) + geom_point() + scale_color_manual(name = bquote(tau), values = c("#74a9cf",
  "#2b8cbe", "#045a8d")) + saveTheme + ylab(bquote(P[t] (wait))) + xlab(bquote(Q[t] (wait, s) - Q(quit)))
ggsave(filename = "figures/plotFigures/softMax.pdf", width = 6, height = 4) 

## plot samples
set.seed(123)
len = 10
seq_ = matrix(NA, len, 3)
for(c in 1:2){
  cond = conditions[c]
  conditionColor = conditionColors[c]
  for(i in 1:3){
    for(j in 1 : len)
      seq_[j,i] = drawSample(cond)
  }
  plotData = data.frame(delay = as.vector(seq_), seq = factor(rep(1:3, len)),
                        trial = rep(1 : len, 3))
  # ggplot(plotData, aes(trial, delay, color = seq)) + geom_point() +
  #   xlab("Trial num") + ylab("Delay / sec")+ saveTheme + geom_line() +
  #   scale_x_continuous(breaks = seq(1, 10, by = 1))
  ggplot(plotData, aes( delay)) + geom_histogram(bins = 5, fill = conditionColor) +
    xlab("Delay / sec") + ylab("Count")+ saveTheme + facet_grid(~seq)
  ggsave(sprintf("figures/sample_%s.pdf", cond), width = 9, height = 3)
}

x = 1 : 30
curSlope = 0.2
curIntercept = 2
y =  curIntercept * exp(-curSlope*(x-1))
plotData = data.frame(trialNum = x, curiosity = y)
ggplot(plotData, aes(trialNum, curiosity)) + geom_point() + saveTheme
ggsave("curiosity.png", width =6, height = 4)

###### plot rewardRate for the new model
conditions = c("HP", "LP")
conditionColors = c("#008837", "#7b3294")

######## timing variables ########
tMaxs = c(20, 40) # trial durations
blockMins = 7 # block duration in mins
blockSecs = blockMins * 60 # block duration in secs
iti = 2 # iti duration in secs
tGrid = seq(0, blockSecs, 0.1)

######### reward variable ########
tokenValue = 10 #value of the token
stepDuration = 1
########## supporting vairbales ########
# time ticks within a trial for timeEarnings or wtw analysis
trialTicks = list(
  'HP' = round(seq(0, tMaxs[1], by = stepDuration), 1),
  'LP' = round(seq(0, tMaxs[2], by = stepDuration), 1)
)

trialGapValues = list(
  "HP" = round(seq(stepDuration, tMaxs[1], by = stepDuration), 2),
  "LP" = round(seq(stepDuration, tMaxs[2], by = stepDuration), 2)
)

trialGapIdxs = list(
  "HP" = 1 : length(trialGapValues$HP),
  "LP" = 1 : length(trialGapValues$LP)
)


########## additional  variables for optimal analysis ########
# CDF of reward delays: p(t_reward <= T)
k = 4 # shape
mu = 0 # location
sigma = 2 # scale
pareto = list()
pareto[['k']] = k # shape
pareto[['mu']] = mu # location
pareto[['sigma']] = sigma # scale
HP = 1 / length(trialGapIdxs$HP) * trialGapIdxs$HP
LP = 1 - (1 + k * (trialGapValues$LP - mu) / sigma) ^ (-1 / k)
LP[length(trialGapValues$LP)] = 1 
rewardDelayCDF = list(
  HP = HP,
  LP = LP
)

#  PDF of reward delays: p(t_reward = T)
#  there is no gap before the first tick, therefore the first element is NaN
#  remmember to use the middle value to the gap when calculating the mean
HP = diff(c(0, rewardDelayCDF$HP))
LP = diff(c(0, rewardDelayCDF$LP))
rewardDelayPDF = list(
  "HP" = HP,
  "LP" = LP
)

# E(delay | elapsed time) for 0 to 19.9 or 39.9, namly at the begining of the gap
# close to the normal analysis so we assume the reward occurs at the middle of the gap
nTimeStep = tMaxs[1] / stepDuration
HP = sapply(1 : nTimeStep, function(i) sum((trialGapValues$HP[i : nTimeStep] - 0.5 * stepDuration - trialTicks$HP[i]) *
                                             rewardDelayPDF$HP[i : nTimeStep])/
              sum(rewardDelayPDF$HP[i : nTimeStep]))

nTimeStep = tMaxs[2] / stepDuration
LP = sapply(1 : nTimeStep, function(i) sum((trialGapValues$LP[i : nTimeStep] - 0.5 * stepDuration - trialTicks$LP[i]) *
                                             rewardDelayPDF$LP[i : nTimeStep])/
              sum(rewardDelayPDF$LP[i : nTimeStep]))
meanRewardDelay = list('HP' = HP, 'LP' = LP)
HP = 1 / (HP + 2)
LP = 1 / (LP + 2)
rewardRate = list('HP' = HP, 'LP' = LP)
plot(rewardRate$HP)
plot(rewardRate$LP)
plot(meanRewardDelay$LP)




# I want to compare the different policy effects
nTimeStep = tMaxs[2] / stepDuration
for(quitTime in 1 : nTimeStep){
  
}


# plot MVT 
ini = 10
gamma = 0.9
lowTh = 6
highTh = 8
t = seq(0, 10, by = 0.1)
r = ini * gamma^(t)
library(ggplot2)
lowColor = "black"
highColor = "#e31a1c"
lowE = t[which.min(abs(r - lowTh))]
highE = t[which.min(abs(r - highTh))]
plotData = data.frame(r = r, t = t)
ggplot(plotData, aes(t, r)) + geom_line(size = 2) + geom_hline(yintercept = lowTh, color = lowColor, size = 2) +
  geom_hline(yintercept = highTh, color = highColor, size = 2) + saveTheme + xlab("Time foraging in patch")+
  ylab("Reward rate") +
  geom_segment(x = lowE, xend = lowE, y = 0, yend = lowTh, color = lowColor, linetype = 2, size = 2) +
  geom_segment(x = highE, xend = highE, y = 0, yend = highTh, color = highColor, linetype = 2, size = 2)
ggsave("MTV.png", width = 5, height = 5)


ggplot(plotData, aes(t, r)) + geom_line(size = 2) + geom_hline(yintercept = lowTh, color = lowColor, size = 2) + 
  saveTheme + xlab("Time foraging in patch")+
  ylab("Reward rate") +
  geom_segment(x = lowE, xend = lowE, y = 0, yend = lowTh, color = lowColor, linetype = 2, size = 2) 
ggsave("MTV.png", width = 5, height = 5)
