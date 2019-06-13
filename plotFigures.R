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


# plot one ini
# so it is a upper
# here we use 0.9 as the discount rate for one stepDuration
QHPApOptim = 5 / 6 * stepDuration / (1 - 0.9) 
QLPApOptim = 0.93 * stepDuration / (1 - 0.9) 
wIni = (QHPApOptim + QLPApOptim)/ 2

zeroPoints = c(15, 30)
Qquit = wIni * 0.9
Viti = wIni * 0.9
Qwait_ = lapply(1:2, function(i){nTimeStep = tMaxs[i] / stepDuration
                zeroPoints[i]*0.05 - 0.05*(0 : (nTimeStep - 1)) + Qquit})
data.frame(Wait = Qwait_[[1]], Quit = rep(Qquit, tMaxs[1]), time = 1 :  tMaxs[1]) %>%
  gather(-c("time"), key = "action", value = "value") %>%
  ggplot(aes(time, value, color = factor(action, labels = c("Quit / Prep", "Wait")))) + geom_line(size = 3) +
  scale_color_manual(values = c("grey", "black")) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title=element_text(size= 30, face = "bold"),
        legend.position=c(0.7,0.8),
        legend.title=element_blank(),
        legend.key = element_rect(fill = NA),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 25, face = "bold")
        )  + ylab("Value") + xlab("Time") 
ggsave("prior.png", width = 5, height = 3)
# plot a softmax function
x = seq(-3, 3, by = 0.1)
tau = 2
y = 1 / (1 + exp(-tau * x))
library("latex2exp")
data.frame(x = x, y = y) %>% ggplot(aes(x, y)) + geom_line(size = 3, color = "#373737" )+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title=element_text(size= 25),
        panel.margin = unit(c(0, 0, 0, 0), "cm")) + ylab(TeX(' '))+ xlab(TeX(''))
ggsave("soft-max.png", height = 2.5, width = 3.5) 

# plot update target 
R = 10
gamma = 0.87
x = seq(0, 10, by = 0.1)
y = gamma^rev(x)* 1
data.frame(x = x, y = y) %>% ggplot(aes(x, y)) + geom_line(size = 2.5) + 
  theme(axis.text = element_text(size = 30, face = "bold"),
        axis.ticks = element_blank(),
        axis.line = element_line(color="#373737", size = 2),
        panel.background = element_rect(fill = "white"),
        axis.title=element_text(size= 30))  + ylab("") + xlab("") +
  scale_x_continuous(labels = c(" ", "T"), breaks = c(0,10),
                     limits = c(-0.5, 10.5)) +
  scale_y_continuous(labels = c("0", "1"), limits = c(0, 1), breaks = c(0,1))
ggsave("discount.png", height = 3, width = 3.5)

# discount2
data.frame(x = x, y = y) %>% ggplot(aes(x, y)) + geom_line(size = 3.5) + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title=element_text(size= 30))  + ylab("") + xlab("") 
ggsave("discount2.png", height = 3, width = 3.5)


# correlation between AUC and triats 
# load personality data
library("Hmisc")
personality = read.csv("data/SDGdataset.csv")
personality$id = personality$SubjectID
traits = c("Delay.of.Gratification", "Barratt.Impulsiveness","Intolerance.of.Uncertainty", "Trait.Anxiety..STAIT.")
traitNames = c("DG", "IMP", "UC", "AX")
nTrait = length(traits)
traitAUCCorr = list()
# plot separately for two conditions
for(i in 1 : nTrait){
  trait = traits[i];
  traitName = traitNames[i]
  input = data.frame(personality[sessionData$stress =="no stress",trait],
                     sessionData$AUC[sessionData$stress == "no stress"],
                     sessionData$condition[sessionData$stress == "no stress"])
  traitAUCCorr[[i]]= getCorrelation(input)
  p = plotCorrelation(input, isRank = T) 
  p + xlab(paste(capitalize(traitName), "(rank)")) + ylab("AUC (rank)") + myTheme
  fileName = sprintf("%s/AUC_%s_%s.png", "figures/expDataAnalysis", traitName, dataType)
  ggsave(fileName, width = 6, height = 3)
}
rhoTable = lapply(1:2, function(j) sapply(1: (nTrait), function(i) traitAUCCorr[[i]]$rhos[j]))
pTable = lapply(1:2, function(j) sapply(1: (nTrait), function(i) traitAUCCorr[[i]]$ps[j]))









