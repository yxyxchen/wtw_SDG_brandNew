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
modelName = "PR"
nTrial = 50

# create output directories
dir.create("figures/simParaRecovery")
tempt = sprintf("figures/simParaRecovery/%s", modelName)
dir.create(tempt)
tempt = sprintf("figures/simParaRecovery/%s/%dTrial", modelName, nTrial)
dir.create(tempt)

# load real parameters
load(sprintf("genData/simulation/%s/%dTrialPara.RData", modelName, nTrial))
# load simPara
paras = getParas(modelName)
nPara = length(paras)
dirName = sprintf("genData/simModelFitting/%s/%dTrial", modelName, nTrial)
simPara = loadSimPara(paras, dirName)

paraNames = c("LR", "LP", expression(tau), expression(gamma), "P")
sumDataPool = data.frame(matrix(NA, ncol = 4+3, nrow = nPara * 4))
colnames(sumDataPool) = c("para", "realValue", "condition", "mean", "ste", "min", "max")
for(cIdx in 1 : 2){
  cond = conditions[cIdx]
  sumData = data.frame(matrix(NA, ncol = 4+3, nrow = nPara * 2))
  colnames(sumData) = c("para", "realValue", "condition", "mean", "ste", "min", "max")
  sumData[,1] = factor(rep(paras, each = 2), levels = paras)
  sumData[,2] = as.double(unlist(paraTable[[cIdx]]))
  sumData[,3] = rep(cond, nPara * 2)
  for(pIdx in 1 : nPara){
    para = paras[pIdx]
    input = data.frame(data = thisSimParaAve[,c(para, paste0(para, "real"))])
    names(input) = c("est", "truth")
    sumData[(pIdx * 2 - 1) : (pIdx * 2),4:7] = input %>% group_by(truth) %>%
      dplyr::summarise(mean = mean(est), ste = sd(est) / sqrt(length(est)),
                       min = mean - ste, max = mean + ste) %>% select(c("mean", "ste", "min", "max"))
  }
  if(cIdx == 1){
    sumDataPool[1 : (2 * nPara), ] = sumData
  }else{
    sumDataPool[(2 * nPara + 1) : (4 * nPara), ] = sumData
  }
}
sumDataPool$realLevel = factor(rep(c(1,2), nPara * 2))
sumDataPool$para = factor(rep(rep(paraNames, each = 2), 2), levels = paraNames)
sumDataPool %>% ggplot(aes(realLevel, mean)) + geom_bar(stat = "identity", fill = "#808080")+
  facet_grid(para ~ condition, scales = "free_y", labeller = label_parsed) +
  myTheme + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  ylab(" ") + geom_point(aes(realLevel, realValue), shape = 17, size = 3, fill = "black") + 
  scale_y_continuous(breaks = scales::pretty_breaks(2))
fileName = sprintf("figures/simParaRecovery/%s/%dTrial/paraRecovery.png", modelName, nTrial)
ggsave(fileName, width = 4, height = 6)

for(cIdx in 1 : 2){
  cond = conditions[cIdx]
  paraComb = getParaComb(paraTable[[cIdx]])
  colnames(paraComb) = paste0(colnames(paraComb), "real")
  thisSimPara = simPara[[cond]][1:nPara, , ]
  thisSimParaAve = cbind(t(apply(thisSimPara, MARGIN = c(1,3), mean)), paraComb)
  for(pIdx in 1 : nPara){
    para = paras[pIdx]
    paraColor = paraColors[pIdx]
    # prepare input 
    input = data.frame(data = thisSimParaAve[,c(para, paste0(para, "real"))])
    names(input) = c("est", "truth")
    # get truth values so we can plot the real values 
    truths = unique(input$truth)
    # 
    wid = abs(diff(range(truths)) / 5)
    thisSize = ifelse(para %in% c("tau", "gamma"), 25, 15)
    # plot
    # maybe delete the error bar?
    input %>% group_by(truth) %>%
      dplyr::summarise(muEst = mean(est),stdEst = sd(est),minEst = mean(est) - sd(est),maxEst = mean(est) + sd(est)) %>%
      ggplot(aes(truth, muEst)) + geom_bar(stat = "identity", fill= paraColor) +
      scale_x_continuous(breaks = truths) + xlab(paraNames[pIdx]) + ylab("Estimation") + 
      geom_point(aes(truth, truth), shape = 17, size = 5, fill = "black") + myTheme +
       theme(axis.title.x = element_text(size = thisSize, face = "bold"))
    fileName = sprintf("figures/simParaRecovery/%s/%dTrial/%s_%s.png", modelName, nTrial, cond, para)
    ggsave(fileName, width = 4, height = 4)
  }
}

# correlation between parameters
library("corrplot")
col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061")))
paraNames = c("LR", "LP", "Tau", "Gamma", "P")
for(cIdx in 1 : 2){
  cond = conditions[cIdx]
  thisSimPara = simPara[[cond]][1:nPara, , ]
  rhoTable = matrix(NA, nrow = nPara, ncol = nPara, dimnames = list(paraNames, paraNames))
  pTable = matrix(NA, nrow = nPara, ncol = nPara)
  
  for(i in 1 : nPara){
    for(j in i : nPara){
      # junk = cor.test(as.vector(thisSimPara[i,,]),
      #                 as.vector(thisSimPara[j,,]),
      #                 method = "kendall") 
      junk1 = cor.test(as.vector(thisSimPara[i,,]),
                      as.vector(thisSimPara[j,,]),
                      method = "spearman")
      junk2 = spearman_test(as.vector(thisSimPara[i,,]) ~ as.vector(thisSimPara[j,,]))
      # pTable[i, j] = junk$p.value
      # rhoTable[i, j] = junk$estimate
      pTable[i, j] = junk1$p.value
      rhoTable[i, j] = pvalue(junk2)
    }
  }
  fileName = sprintf("figures/simParaRecovery/%s/%dTrial/sp_cor_%s.png", modelName, nTrial, cond)
  png(fileName)
  corrplot(rhoTable, 
           p.mat = pTable, 
           is.corr = T, 
           method = "color",
           insig = "label_sig",
           tl.col = "black", tl.srt = 15, tl.cex = 1.5,
           type = "upper",
           col = col2(50)) 
  dev.off()
  
}


