load("expParas.RData")
library("ggplot2"); library("Hmisc"); library("coin")
library("dplyr"); library("tidyr")
source("subFxs/plotThemes.R")
source("subFxs/loadFxs.R") # load blockData and expPara
source("subFxs/helpFxs.R") # getparaNames
source("subFxs/analysisFxs.R") # plotCorrelation and getCorrelation
source('MFAnalysis.R')

# model Name
modelName = "QL1"
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
# paraNames = c("LR", "LP", expression(tau), expression(gamma), "P")
# paraNames = c("LR", "LP", expression(tau), "P")
paraNames = paraNames
expPara$condition = sumStats$condition[sumStats$id %in% expPara$id]
expPara %>% filter(passCheck ) %>% select(c(paraNames, "condition")) %>%
  gather(-c("condition"), key = "para", value = "value") %>%
  mutate(para = factor(para, levels = paraNames, labels = paraNames ))%>%
  ggplot(aes(value)) + geom_histogram(bins = 8) +
  facet_grid(condition ~ para, scales = "free", labeller = label_parsed) + 
  myTheme + xlab(" ") + ylab(" ")
fileName = sprintf("%s/%s/hist.pdf", "figures/expParaAnalysis", modelName)
ggsave(fileName, width = 8, height = 4)

# summary stats for expPara
expParaInfo = expPara %>% filter(passCheck) %>% select(c(paraNames)) %>%
  gather(key = "para", value = "value") %>%
  mutate(para = factor(para, levels = paraNames, labels = paraNames ))%>% 
  group_by(para) %>% summarise(mu = mean(value), median = median(value),
                               se = sd(value) / sqrt(length(value)))

dir.create("genData/expParaAnalysis")
dir.create(sprintf("genData/expParaAnalysis/%s", modelName))
save(expParaInfo, file = sprintf("genData/expParaAnalysis/%s/expParaInfo.RData", modelName))

shape_ = expPara %>% filter(passCheck & condition == "HP") %>% select(c(paraNames)) %>%
  mutate(phi = phi / max(phi),
         tau = tau / max(tau),
         gamma = (gamma - min(gamma)) / (max(gamma) - min(gamma)),
         prior = (prior - min(prior)) / (max(prior) - min(prior))) %>%
  gather(key = "para", value = "value") %>%
  mutate(para = factor(para, levels = paraNames, labels = paraNames )) %>%
  group_by(para) %>% 
  summarise(mu = mean(value),
            var = var(value),
            alpha = ((1 - mu) / var - 1 / mu) * mu ^ 2,
            beta = alpha * (1 / mu - 1)) 

scale_ = expPara %>% filter(passCheck & condition == "HP") %>% select(c(paraNames)) %>%
  gather(key = "para", value = "value") %>%
  mutate(para = factor(para, levels = paraNames, labels = paraNames )) %>%
  group_by(para) %>% summarise(a = min(value), c = max(value))
i = 4
paraName = paraNames[i]
values = expPara[passCheck & expPara$condition == "HP", paraName]
x = seq(0, 1, length.out = 10)
cdf = pbeta(x, shape_$alpha[1], shape_$beta[i])
pdf = c(0, diff(cdf))
length(pdf)
xPrime = x*(scale_$c[i] - scale_$a[i]) + scale_$a[i]
hist(values, breaks = xPrime)
lines(xPrime, pdf * length(values))

# optimism bias
wilcoxResults = wilcox.test(expPara$phi_pos[passCheck] - expPara$phi_neg[passCheck])
# we use the 1.5 x the IQR for the larger of the two values 
# (phi_pos, phi_neg) as the criteria 
# for exclude outliers
junk = expPara %>% filter(passCheck) %>%
  select(c("phi_pos", "phi_neg")) %>%
  gather("paraName", "paraValue") %>%
  group_by(paraName) %>% 
  summarise(qLower = quantile(paraValue)[1],
            qUpper = quantile(paraValue)[3],
            IQR = qUpper - qLower,
            limLower = qLower - IQR * 1.5,
            limUpper = qUpper + IQR * 1.5)
limLower = max(0, min(junk$limLower))
limUpper = max(junk$limUpper)
# optimism bias
expPara %>% filter(passCheck) %>% 
  ggplot(aes(phi_neg, phi_pos)) +
  geom_point(color = themeColor,  fill = "#fdd49e", shape= 21, stroke = 1, size = 5) + 
  xlim(c(limLower, limUpper)) + ylim(c(limLower, limUpper)) + 
  geom_abline(slope = 1, intercept = 0) + 
  annotate("text", x = 0.015, y = 0.015, label = sprintf("p = %.3f", wilcoxResults$p.value)) +
  myTheme
ggsave("figures/expParaAnalysis/optimism.eps", width = 6, height = 6)


# load and merge trait data
personality = read.csv("data/hdrData.csv")
personality = personality[personality$id %in% expPara$id,]
personality$id = expPara$id
traits = c("BDI","IUS","DoG","BIS.11") #BIS_11 trait anxiety 
nTrait = length(traits)
expPara = merge(x= expPara,y = personality[,c(traits, "id")], by="id",all.x=TRUE)
expPara$optimism = expPara$phi_pos /  expPara$phi_neg
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
