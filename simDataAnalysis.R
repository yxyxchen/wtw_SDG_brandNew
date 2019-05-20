library("ggplot2")
library("dplyr")
library("Hmisc")
library("coin")
source("subFxs/plotThemes.R")
source("subFxs/loadFxs.R") # load blockData and expPara
source("subFxs/helpFxs.R") # getParas
source("subFxs/analysisFxs.R") # plotCorrelation and getCorrelation
load("wtwSettings.RData")

# input 
modelName = "curiosityTrialSp"
nTrail = 50

# create output directories
dir.create("figures/simParaAnalysis")
tempt = sprintf("figures/simParaAnalysis/%s", modelName)
dir.create(tempt)
saveDir = sprintf("figures/simParaAnalysis/%s/%dTrial", modelName, nTrial)
dir.create(tempt)

# load real parameters
load(sprintf("genData/simulation/%s/%dTrialPara.RData", modelName, nTrial))

# load simPara
paras = getParas(modelName)
nPara = length(paras)
dirName = sprintf("genData/simModelFitting/%s/%dTrial", modelName, nTrial)
tempt = loadExpPara(paras, dirName)
useID = getUseID(tempt, paras)
## only keep the estimations 
tempt = tempt[,1:nPara]
colnames(paraComb) = paste0(paras, "Real")
summaryData = cbind(rbind(paraComb, paraComb), condition = rep(c("HP", "LP"),each = nrow(paraComb)))
simPara = cbind(tempt, summaryData)

simPara$phiReal = as.factor(simPara$phiReal)
ggplot(simPara[simPara$condition == "HP", ], aes(gammaReal, gamma)) +
  geom_boxplot()

simParaLong = data.frame(realValue = unlist(simPara[,paras]),
                         estValue = unlist(simPara[,paste0(paras, "Real")]),
                         conditon = rep(conditions, each = nrow(simPara) * nPara / 2),
                         para = rep(rep(paras, each = nrow(simPara) / 2),2)
                   )
ggplot(simParaLong[simParaLong$conditon == "HP", ],
       aes(realValue, estValue)) + geom_boxplot() + facet_grid(~para)
  
# compare real parameter and estimated parameters 
ggplot(simPara, aes(realPhi, ))


## look at corrletion
modelName = "full_model"
paras = getParas(modelName )

load(sprintf("genData/simulation/%s/simParas.RData", modelName))
for(c in 1:2){
  cond = conditions[c]
  condColor = conditionColors[c]
  simPara_ = loadSimPara_(modelName, paras, cond)
  simPara = matrix(simPara_, 810, 4*5)
  if(c == 1){
    simParaHP = simPara
  }else{
    simParaLP = simPara
  }
}
simPara = rbind(simParaHP, simParaLP)
simPara = as.data.frame(simPara)
colnames(simPara) = c(paras,  paste0(paras, "SD"), paste0(paras, "Effe"), paste0(paras, "Rhat"),
                      paste0(paras, "real"))
simPara$condition = rep(c("HP", "LP"), each = nrow(simParaHP))

# choose only useful data
RhatCols = which(str_detect(colnames(simPara), "hat"))[1 : length(paras)]
EffeCols = which(str_detect(colnames(simPara), "Effe"))[1 : length(paras)]
simParaUse = simPara[apply(simPara[,RhatCols] < 1.1, MARGIN = 1, sum) == length(paras) & 
                 apply(simPara[,EffeCols] >100, MARGIN = 1, sum) == length(paras), ]

plotParaPara = function(expPara, paraX, paraY){
  # calculate correlations
  corHP = cor.test(expPara[expPara$condition == "HP", paraX], expPara[expPara$condition == "HP", paraY], method = "spearman")
  corLP = cor.test(expPara[expPara$condition == "LP", paraX], expPara[expPara$condition == "LP", paraY], method = "spearman")
  rhoHP = round(corHP$estimate, 3)
  rhoLP= round(corLP$estimate, 3)
  pHP = round(corHP$p.value, 3)
  pLP = round(corLP$p.value, 3)
  textColors = c(ifelse(pHP < 0.05, "red", "blue"), ifelse(pLP < 0.05, "red", "blue"))
  textData = data.frame(label = c(paste(rhoHP, "(p =", pHP, ")"), paste(rhoLP, "(p =", pLP, ")")),
                        condition = c("HP", "LP"))
  plotData = expPara
  p = ggplot(plotData, aes_string(paraX, paraY)) + geom_point() +
    facet_grid(~condition)  +
    geom_text(data = textData,aes(x = -Inf,y = -Inf, label = label),
              hjust   = -0.2,vjust = -1,color = "blue",size = 5, fontface = 2) +
    saveTheme + ylab(capitalize(paraY)) + xlab(capitalize(paraX))
  print(p)
}

for(i in 1 : (nPara-1)){
  paraX = paras[i]
  for(j in (i+1): nPara){
    paraY = paras[j]
    plotParaPara(simParaUse, paraX, paraY)
    fileName = sprintf("figures/simDataAnalysis/cor_%s_%s.pdf", paraX, paraY)
    ggsave(fileName, width = 6, height = 4)
  }
}

library("corrplot")
tempt = simParaHP[,1:4]
colnames(tempt) = paras
M = cor(tempt, method = "spearman")
pdf('figures/simDataAnalysis/corrplotHP.pdf')
corrplot(M, method = "color", type = "upper", outline = T, tl.col = "black",
          tl.cex = 2, cl.pos = "b", cl.cex = 1.2) 
dev.off()

tempt = simParaLP[,1:4]
colnames(tempt) = paras
M = cor(tempt, method = "spearman")
pdf('figures/simDataAnalysis/corrplotLP.pdf')
corrplot(M, method = "color", type = "upper", outline = T, tl.col = "black",
         tl.cex = 2, cl.pos = "b", cl.cex = 1.2) 
dev.off()

c= 1
cond = conditions[c]
condColor = conditionColors[c]
tempt = simParaUse[simParaUse$condition == "LP", ]
##### recaluclate the recover accuracy
for(i in 1 : nPara){
  para = paras[i]
  realName = sprintf("%sreal", para)
  tempt1 = tempt %>% group_by_at(vars(realName)) %>% summarise(mu = mean(get(para)))
  tempt2 = tempt %>% group_by_at(vars(realName)) %>% summarise(std = sd(get(para)))
  plotData = data.frame(tempt1, std = tempt2$std) 
  plotData$ymin = plotData$mu - plotData$std
  plotData$ymax = plotData$mu + plotData$std
  if(para == "phi" || para == "gamma"){
    wid = 0.01
  }else if(para == "tau" ){
    wid = 2
  }else wid = 0.5
  
  ggplot(plotData, aes_string(realName, "mu")) +
    geom_bar(stat = "identity", fill = condColor,color = condColor) +
    saveTheme + xlab(capitalize(para)) + ylab("Phi estimation")  +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = wid) +
    geom_point(aes_string(realName, realName), shape = 17, size = 5, fill = "black")
  fileName = sprintf("figures/simDataAnalysis/recovery_%s_%s.pdf", cond, para)
  ggsave(fileName, width = 4, height = 4)
}
