load("expParas.RData")
library("ggplot2"); library("Hmisc"); library("coin")
library("dplyr"); library("tidyr")
source("subFxs/plotThemes.R")
source("subFxs/loadFxs.R") # load blockData and expPara
source("subFxs/helpFxs.R") # getparaNames
source("subFxs/analysisFxs.R") # plotCorrelation and getCorrelation
source('MFAnalysis.R')
library(latex2exp)

# load empirical data 
allData = loadAllData()
hdrData = allData$hdrData 
trialData = allData$trialData       
ids = hdrData$id[hdrData$stress == "no_stress"]
nSub = length(ids)

# model Name
modelName = "optim_noise_bias"
paraNames = getParaNames(modelName)
nPara = length(paraNames)

# output directories
dir.create("figures/expParaAnalysis")
saveDir = sprintf("figures/expParaAnalysis/%s", modelName)
dir.create(saveDir)

# 
MFResults = MFAnalysis(isTrct = T)
sumStats = MFResults[['sumStats']]

# load expPara
paraNames = getParaNames(modelName)
nPara = length(paraNames)
parentDir = "genData/expModelFit"
dirName = sprintf("%s/%s",parentDir, modelName)
tempt = loadExpPara(paraNames, dirName)
passCheck = checkFit(paraNames, tempt)
expPara = merge(x=tempt[,c(paraNames, "id")],y=sumStats, by="id",all.x=TRUE)

# plot hist 
paraLabels = c(expression(alpha_R), expression(alpha_U), expression(tau), expression(gamma), expression("eta"))
expPara$condition = sumStats$condition[sumStats$id %in% expPara$id]
expPara %>% filter(passCheck ) %>% select(c(paraNames, "condition")) %>%
  gather(-c("condition"), key = "para", value = "value") %>%
  mutate(para = factor(para, levels = paraNames, labels = paraLabels ))%>%
  ggplot(aes(value)) + geom_histogram(bins = 12) +
  facet_grid( ~ para, scales = "free", labeller = label_parsed) + 
  myTheme + xlab(" ") + ylab(" ")
fileName = sprintf("%s/%s/hist.pdf", "figures/expParaAnalysis", modelName)
ggsave(fileName, width = 8, height = 4)

# optimism bias
nQuits = sapply(1 : nSub, function(i) sum(trialData[[ids[i]]]$trialEarnings == 0))
wilcoxResults = wilcox.test(expPara$alphaR[passCheck & enoughQuits] / expPara$alphaU[passCheck & enoughQuits] - 1)

logOdds = log(expPara$alphaR / expPara$alphaU)
wilcoxResults = wilcox.test(logOdds[expPara$condition == "HP"])
wilcoxResults = wilcox.test(expPara$alphaR[expPara$condition == "HP" & expPara], expPara$alphaU[expPara$condition == "HP"])
wilcoxResults = wilcox.test(logOdds[expPara$condition == "LP"])

# optimism bias
expPara %>% filter(passCheck) %>%
  ggplot(aes(log(alphaR/alphaU))) +
  geom_histogram(bins = 8) +
  myTheme + facet_grid(~condition) + 
  xlab(TeX('$log(\\alpha_r/\\alpha_u)$')) +
  ylab("Count") +
  geom_vline(aes(xintercept = 0), color = "red", linetype = 2)
ggsave("figures/expParaAnalysis/optimism.eps", width = 6, height = 6)
ggsave("figures/expParaAnalysis/optimism.png", width = 6, height = 6)

# temproal discounting
wilcoxResults = wilcox.test(expPara$gamma[expPara$condition == "HP"] - 1)
wilcoxResults = wilcox.test(expPara$gamma[expPara$condition == "LP"] - 1)
expPara %>% filter(passCheck) %>% ggplot(aes(gamma)) +
  geom_histogram(bins = 8) +
  myTheme + facet_grid(~condition) + 
  xlab(TeX('$\\gamma$')) +
  ylab("Count") + xlim(c(0.65, 1.05))
ggsave("figures/expParaAnalysis/discounting.eps", width = 4, height = 4)
ggsave("figures/expParaAnalysis/discounting.png", width = 4, height = 4)

##### load and merge trait data
personality = read.csv("data/hdrData.csv")
personality = personality[personality$id %in% expPara$id,]
personality$id = expPara$id
traits = c("BDI","IUS","DoG","BIS.11") #BIS_11 trait anxiety 
nTrait = length(traits)
expPara = merge(x= expPara,y = personality[,c(traits, "id")], by="id",all.x=TRUE)
expPara$optimism = expPara$phi_pos /  expPara$phi_neg

# 
plot(log(expPara$tau), expPara$stdWTW)
# calculate trait-para correlations
rhoHP_ = matrix(NA, nrow = 1 + nPara, ncol = nTrait)
rhoLP_ = matrix(NA, nrow = 1 + nPara, ncol = nTrait)
pHP_ = matrix(NA, nrow = 1 + nPara, ncol = nTrait)
pLP_ = matrix(NA, nrow = 1 + nPara, ncol = nTrait)
for(pIdx in 1 : (nPara+1)){
    if(pIdx <= nPara){
      para = paraNames[pIdx]
    }else{
      para = "optimism"
    }
    
  for(trIdx in 1 : length(traits)){
    trait = traits[trIdx]
    # for HP 
    corResults = cor.test(expPara[expPara$condition == "HP" & passCheck,trait],
                          expPara[expPara$condition == "HP" & passCheck,para],
                          method = "kendall")
    rhoHP_[pIdx, trIdx] = as.double(corResults$estimate)
    pHP_[pIdx, trIdx] = as.double(corResults$p.value)
    # for LP 
    corResults = cor.test(expPara[expPara$condition == "LP" & passCheck,trait],
                          expPara[expPara$condition == "LP" & passCheck,para],
                          method = "kendall")
    rhoLP_[pIdx, trIdx] = as.double(corResults$estimate)
    pLP_[pIdx, trIdx] = as.double(corResults$p.value)
  }
}

pHP_
pLP_
