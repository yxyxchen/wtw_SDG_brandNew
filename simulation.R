
modelName = "QL1"

# create output directories
dir.create("figures/simulation/")
dir.create(sprintf("figures/simulation/%s",modelName))

# load experiment parameters
load('expParas.RData')

# load packages and sub functions 
library("ggplot2")
library("tidyr")
library("dplyr")
source("subFxs/plotThemes.R")
source("subFxs/helpFxs.R") 
source("subFxs/loadFxs.R") # 
source("subFxs/analysisFxs.R") 
source("subFxs/taskFxs.R")

# get the generative model 
source(sprintf("subFxs/gnrModels/%s.R", modelName))
gnrModel = get(modelName)
paraNames = getParaNames(modelName)
nPara = length(paraNames)

# num of repetitions 
nRep = 10

######################### simulate one typical example #####################
# specify the learning parameters
alpha = 0.03
tau = 3
gamma = 0.9
prior = 4
paras = c(alpha, tau, gamma, prior)
# initialize outputs
simHPExample = list(length = nRep)
simLPExample = list(length = nRep)
# number of trials
nTrial = 100
# simulate 
for(rIdx in 1 : nRep){
  simHPExample[[rIdx]] = gnrModel(paras, "HP", replicate(nTrial, drawSample("HP")))
  simLPExample[[rIdx]] = gnrModel(paras, "LP", replicate(nTrial, drawSample("LP")))
}
# summarise simulation results
exampleAUCs = vector(length = nRep * length(conditions))
exampleSCs = matrix(NA, length(kmGrid), nRep * length(conditions))
for(rIdx in 1 : nRep){
  kmscResults = kmsc(simHPExample[[rIdx]], min(delayMaxs), F, kmGrid)
  exampleAUCs[rIdx] = kmscResults$auc  
  exampleSCs[, rIdx] = kmscResults$kmOnGrid
  
  kmscResults = kmsc(simLPExample[[rIdx]], min(delayMaxs), F, kmGrid)
  exampleAUCs[rIdx + nRep] = kmscResults$auc  
  exampleSCs[, rIdx + nRep] = kmscResults$kmOnGrid
}
optimDf = data.frame(
  condition = rep(conditions, each = length(kmGrid)),
  t = rep(kmGrid, 2),
  pSurvival = 1
)
optimDf$pSurvival[optimDf$condition == "LP" & optimDf$t > optimWaitThresholds$LP] = 0
data.frame(
  pSurvival = as.vector(exampleSCs),
  t = rep(kmGrid, nRep * length(conditions)),
  condition = rep(rep(conditions, each = nRep), each = length(kmGrid))
) %>% group_by(condition, t) %>%
  summarise(
    mu = mean(pSurvival),
    se = sd(pSurvival) / sqrt(length(pSurvival) - 1),
    ymin = mu - se,
    ymax = mu + se
  ) %>%
  ggplot(aes(t, mu)) + geom_line(aes(color = condition), size = 1) +
  geom_ribbon(aes(x = t, ymin = ymin, ymax = ymax, fill = condition), alpha = 0.5) +
  geom_line(aes(t, pSurvival, color = condition), data = optimDf, inherit.aes = F,
            size = 1, linetype = "dashed") +
  myTheme + scale_color_manual(values = conditionColors) +
  scale_fill_manual(values = conditionColors) +
  xlab("Elapsed time (s)") + ylab("Survival rate") +
  theme(legend.title = element_blank())
ggsave(sprintf("figures/simulation/%s/survival.png",modelName), width = 4, height = 3)
  
######################### simulate for different parameters #####################


