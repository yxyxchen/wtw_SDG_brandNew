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
    # plot
    # maybe delete the error bar?
    input %>% group_by(truth) %>%
      dplyr::summarise(muEst = mean(est),stdEst = sd(est),minEst = mean(est) - sd(est),maxEst = mean(est) + sd(est)) %>%
      ggplot(aes(truth, muEst)) + geom_bar(stat = "identity", fill= paraColor) +
      theme_linedraw(base_size = 13) + xlab(capitalize(para)) + ylab("Estimation") + 
      geom_point(aes(truth, truth), shape = 17, size = 5, fill = "black")
    fileName = sprintf("figures/simParaRecovery/%s/%dTrial/%s_%s.png", modelName, nTrial, cond, para)
    ggsave(fileName, width = 4, height = 4)
  # back up version
    # ggplot(aes(truth, muEst)) + geom_bar(stat = "identity", fill= paraColor) +
    #   geom_errorbar(aes(ymin = minEst, ymax = maxEst), width = wid) +
    #   geom_errorbar(aes(ymin=truths, ymax=truths), colour="grey", linetype = 1, size = 1.5) +
    #   theme_linedraw(base_size = 13) + xlab(capitalize(para)) + ylab("Estimation") + 
    #   geom_point(aes(truth, truth), shape = 17, size = 5, fill = "black")
  }
}


