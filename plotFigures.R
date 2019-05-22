library('ggplot2')
library('plyr')
library('dplyr')
library('tidyr')
load("wtwSettings.RData")
source('subFxs/repetitionFxs.R') # called by simulate 
source("subFxs/helpFxs.R") # getParas
source("subFxs/plotThemes.R")
source("subFxs/loadFxs.R") # load scheduledWait from empirical data
source("subFxs/analysisFxs.R") 
source("subFxs/taskFxs.R") # drawSample
source("subFxs/simulationFxs.R")

# plot the example simulation to demonstrate RL
paras = c(0.03, 3.5, 0.90, 35)
modelName = "curiosityTrialSp"
repModelFun = getRepModelFun(modelName)
cond = "LP"
nTimeStep = tMaxs[2] / stepDuration
tMax = tMaxs[2]
lenSeq = 100
set.seed(123)
scheduledWait = unlist(lapply(1:lenSeq, function(x) drawSample(cond)))
tempt = repModelFun(paras, cond, scheduledWait)
trialPlots(tempt, cond)

# 
trials = c(1,2,100)
data.frame(Wait = as.vector(tempt$Qwaits[,trials]),
           Quit = rep(tempt$Qquits[trials], each = nTimeStep),
           trial = rep(trials, each = nTimeStep),
           timeStep = rep(1 : nTimeStep, length(trials))) %>% 
  gather(-c("trial", "timeStep"), key = "action", value = "value") %>%
ggplot(aes(timeStep, value)) + geom_line(aes(color = action)) + facet_grid(trial~.) +
  scale_color_manual(values = c("blue", "red"))

# how do I calculate the optimal values


# 
gamma = 0.80
xCon = seq(1,5,by = 0.1)
yCon = 18.5 * gamma ^ (rev(xCon)) 
x = xCon
y = yCon
x[! (xCon %in% 1:5)] = NA
y[! (xCon %in% c(1, 3:5))] = NA
ggplot(data.frame(x = x, y = y), aes(x, y)) + geom_point(size = 3) + saveTheme + xlab("") + ylab("") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(filename = "discount.png", width =6, height = 2)
