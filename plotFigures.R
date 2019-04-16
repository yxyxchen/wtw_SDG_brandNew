library("ggplot2")
load("wtwSettings.RData")
source("subFxs/plotThemes.R")
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

