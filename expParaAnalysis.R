load("expParas.RData")
library("ggplot2"); library("Hmisc"); library("coin")
library("dplyr"); library("tidyr")
source("subFxs/plotThemes.R")
source("subFxs/loadFxs.R") # load blockData and expPara
source("subFxs/helpFxs.R") # getparaNames
source("subFxs/analysisFxs.R") # plotCorrelation and getCorrelation
source('MFAnalysis.R')
library("EnvStats")

# model Name
modelName = "QL1"
paraNames = getParaNames(modelName)
nPara = length(paraNames)

# output directories
dir.create("figures/expParaAnalysis")
saveDir = sprintf("figures/expParaAnalysis/%s", modelName)
dir.create(saveDir)

# read data 
post = read.csv("genData/expModelFit/QL2/s1.txt", header = F)
mus = apply(post[,1 : nPara], MARGIN = 2, FUN = mean)
muData = data.frame(para = paraNames, mu = mus)
plotData = data.frame(post) 
names(plotData) =  c("LR", "LP", "tau", "gamma", "P", "LL")
plotData %>% gather(key = "para", value = "value", -LL) %>%
  mutate(para = factor(para, levels =  c("LR", "LP", "tau", "gamma", "P")),
         mu = rep(mus, each = 16000)) %>%
  ggplot(aes(value)) + geom_histogram() + 
  facet_wrap(~para, scales = "free", labeller = label_parsed) +
  myTheme + xlab("") + ylab("") + geom_vline(aes(xintercept = mu), color = "red")
ggsave("figures/post.png", width = 10, height = 5)
  

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
  ggplot(aes(value)) + geom_histogram(bins = 12) +
  facet_grid( ~ para, scales = "free", labeller = label_parsed) + 
  myTheme + xlab(" ") + ylab(" ")
fileName = sprintf("%s/%s/hist.pdf", "figures/expParaAnalysis", modelName)
ggsave(fileName, width = 8, height = 4)

# try model fit tau
x = expPara$gamma
up = 0.7
low = 1
breaks = seq(up, low, length.out = 10)
sampleQts = quantile(log(x), seq(0.001, 0.999, length.out = 10))
theoryQts = qnorm(seq(0.001, 0.999, length.out = 10), mean(log(x)), sd(log(x)))
plot(theoryQts, sampleQts)
abline(coef = c(0,1))

breaks = seq(up, low, length.out = 20)
hts = hist(x, breaks = breaks)
cumprobs = plnorm(breaks, meanlog = mean(log(x)), sdlog = sd(log(x)))
lines(hts$mids, diff(cumprobs) * nrow(expPara), col = "green")
cumprobs = pnorm(breaks, mean(x), sd(x))
lines(hts$mids, diff(cumprobs) * nrow(expPara), col = "red")
# pareto distribution 
n = nrow(expPara)
m = min(x)
alpha = n / sum(log(x) - log(m))
cumprobs = ppareto(breaks, m, alpha)
lines(hts$mids, diff(cumprobs) * nrow(expPara), col = "blue")
# gamma distribution
theta = var(x) / mean(x)
k = mean(x) / theta
cumprobs = pgamma(breaks, shape =  k, scale = theta)
lines(hts$mids, diff(cumprobs) * nrow(expPara), col = "pink")
# try model fit 

# summary stats for expPara
expParaInfo = expPara %>% filter(passCheck) %>% select(c(paraNames)) %>%
  gather(key = "para", value = "value") %>%
  mutate(para = factor(para, levels = paraNames, labels = paraNames ))%>% 
  group_by(para) %>% dplyr::summarise(mu = mean(value), median = median(value),
                               se = sd(value) / sqrt(length(value)),
                               log.mu = mean(log(value)), log.sd = sd(log(value)))


dir.create("genData/expParaAnalysis")
dir.create(sprintf("genData/expParaAnalysis/%s", modelName))
save(expParaInfo, file = sprintf("genData/expParaAnalysis/%s/expParaInfo.RData", modelName))



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
