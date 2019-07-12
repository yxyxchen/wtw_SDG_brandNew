gamma = 0.9
totalT = 10
# R-learning treats prereward and postreward delay equally
# while Q-learning focus more on prereward delay.
R = 2
trialT = 2
preRdelays = seq(0, 2, by = 0.4)  # pre-reward delays 
postRdelays = trialT - preRdelays # post-reward delays

sumRs = sapply(1 : length(preRdelays),
               function(i) R * (1 - (gamma^trialT)^(totalT / trialT)) / (1 - (gamma^trialT)) *
                 gamma ^ (preRdelays[i]))
plot(preRdelays, sumRs, "tyl"= 2, ylab = "Action value", 
     xlab = "Pre-reward delay (a.u.)")

# a numerical example, where Qlearning and Rlearning make different decisions 
# because Qlearning focus more on prereward delay
R1 = 2
trialT1 = 2
preRdelay1 = 0

R2 = 2.2
trialT2 = 2
preRdelay2 = 0



smallR = 1
smallT = 1
smallGa = gamma ^ smallT

largeR = 2
largeT = 2
largeGa = gamma ^ largeT

preRtimes = seq(0, 1, by = 0.2)

# no preReward time                                                                                                                                                                                           
smallSum = smallR * (1 - smallGa ^ 10) / (1 - smallGa)
largeSum = largeR * (1 - largeGa ^ 5) / (1 - largeGa)

smallSums = sapply(1 : length(preRtimes),
                   function(i) smallSum * gamma ^ (preRtimes[i]))
largeSums = sapply(1 : length(preRtimes),
                   function(i) largeSum * gamma ^ (preRtimes[i]))

library("ggplot2")
source("subFxs/plotThemes.R")
data.frame(smallSum = smallSums, largeSum = largeSums) %>%
  ggplot(aes(smallSum, largeSum)) + geom_point() +
  geom_abline(slope = 1, xintercept = 0) + ylim(c(5.5, 7)) +
  xlim(c(5.5, 7))

