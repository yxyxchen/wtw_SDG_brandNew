library("ggplot2")
library("dplyr")
library("tidyr")
library("Hmisc")
library("coin")
source("subFxs/plotThemes.R")
source("subFxs/loadFxs.R") # load blockData and expPara
source("subFxs/helpFxs.R") # getParas
source("subFxs/analysisFxs.R") # plotCorrelation and getCorrelation
load("wtwSettings.RData")
modelName = "RL2"
# create output directories
dir.create("figures/expParaAnalysis")
saveDir = sprintf("figures/expParaAnalysis/%s", modelName)
dir.create(saveDir)


# load blockdata data
load("genData/expDataAnalysis/sessionData.RData")
load("genData/expDataAnalysis/kmOnGridSess.RData")
summaryData = sessionData

# load expPara
paras = getParas(modelName)
nPara = length(paras)
parentDir = "genData/expModelFitting"
dirName = sprintf("%s/%sdb",parentDir, modelName)
tempt = loadExpPara(paras, dirName)
useID = getUseID(tempt, paras)
expPara = merge(x=tempt[,c(paras, "id")],y=summaryData, by="id",all.x=TRUE)

# plot hist 
# paraNames = c("LR", "LP", expression(tau), expression(gamma), "P")
# paraNames = c("LR", "LP", expression(tau), "P")
paraNames = paras
expPara$condition = summaryData$condition[summaryData$id %in% expPara$id]
expPara %>% filter(id %in% useID) %>% select(c(paras, "condition")) %>%
  gather(-c("condition"), key = "para", value = "value") %>%
  mutate(para = factor(para, levels = paras, labels = paraNames ))%>%
  ggplot(aes(value)) + geom_histogram(bins = 8) +
  facet_grid(condition ~ para, scales = "free", labeller = label_parsed) + 
  myTheme + xlab(" ") + ylab(" ")
fileName = sprintf("%s/%s/hist.pdf", "figures/expParaAnalysis", modelName)
ggsave(fileName, width = 8, height = 4)


# load and merge trait data
personality = read.csv("data/SDGdataset.csv")
personality$id = personality$SubjectID
traits = c("Delay.of.Gratification", "Barratt.Impulsiveness",
           "Intolerance.of.Uncertainty", "Trait.Anxiety..STAIT.")
traitNames = c("DelayGra", "Impulsive",
               "Uncertain", "Anxiety")
nTrait = length(traits)
expPara = merge(x= expPara,y=personality[,c(traits, "id")], by="id",all.x=TRUE)
names(expPara)[which(names(expPara) %in% traits)] = traitNames
expPara$deltaPhi = expPara$phiP - expPara$phi

# calculate trait-para correlations
traitParaCorr = vector(mode = "list", length = (length(paras) + 1) * length(traits))
corrNo = matrix(1:((nPara + 1) * nTrait), nrow = nPara + 1, ncol = nTrait)
for(pIdx in 1 : (length(paras) + 1)){
  if(pIdx <= length(paras)){
    para = paras[pIdx]
  }else{
    para = "deltaPhi"
  }
  
  for(trIdx in 1 : length(traits)){
    trait = traits[trIdx]
    traitName = traitNames[trIdx]
    # combine persondatlaity, expPara and condtion. useID only
    input = data.frame(personality[summaryData$id %in% useID,trait], expPara[expPara$id %in% useID,para],
                       summaryData$condition[summaryData$id %in% useID])
    traitParaCorr[[corrNo[pIdx, trIdx]]]= getCorrelation(input)
    # plot
    # p = plotCorrelation(input, paraColor, T) 
    # p + ylab(capitalize(para)) + xlab(traitName) + saveTheme
    # parentDir = ifelse(dataType == "block", "figures/expParaAnalysis", "figures/expParaAnalysisSub/")
    # fileName = sprintf("%s/%s/%s_%s.png", parentDir, modelName, para, traitName)
    # ggsave(fileName, width = 6, height = 3)
  }
}
# paraNames = c("LR", "LP", "Tau", "Gamma", "P", "deltaL")
paraNames = c(paras, "deltaL")
dimNames = list(paraNames, c("DG", "IM", "IU", "AX"))
rhoTable = lapply(1:2, function(j) matrix(sapply(1: (nTrait * (nPara+1)), function(i) traitParaCorr[[i]]$rhos[j]),
                                          nrow = nPara+1, dimnames = dimNames))
pTable = lapply(1:2, function(j) matrix(sapply(1: (nTrait * (nPara+1)), function(i) traitParaCorr[[i]]$ps[j]),
                                        nrow = nPara+1, dimnames = dimNames))


# plot trait analysis
library("corrplot")
col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061")))
parentDir = ifelse(dataType == "sess", "figures/expParaAnalysisSub", "figures/expParaAnalysis")
for(i in 1 : 2){
  cond = conditions[i]
  fileName = sprintf("%s/%s/traitPara%s.png", "figures/expParaAnalysis/", modelName, cond)
  png(fileName)
  corrplot(rhoTable[[i]], 
           p.mat = pTable[[i]], 
           is.corr = T, 
           method = "color",
           #insig = "label_sig",
           tl.col = "black", tl.srt = 15, tl.cex = 1.5,
           col = col2(50)) 
  # mtext("Parameter", side = 2,  cex = 2)
  # mtext("Trait", side = 3, cex = 2)
  dev.off()
}

# look at trait linearly
for(cIdx in 1 : 2){
  cond = conditions[cIdx]
  expPara[expPara$condition == cond,c(paras, traitNames)] %>%
    gather(-c(paras), key = "trait", value = "traitValue") %>%
    gather(-c("trait", "traitValue"), key = "para", value = "paraValue") %>%
    mutate(trait = factor(trait, levels = traitNames), para = factor(para, levels = paras, labels = capitalize(paras))) %>%
    ggplot(aes(x = traitValue, y = paraValue)) + geom_point() +
    facet_grid(para ~ trait, scales = "free") + theme_linedraw(base_size = 13)+
    xlab("") + ylab("") + ggtitle(cond) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
  parentDir = "figures/expParaAnalysis"
  fileName = sprintf("%s/%s/lm_%s.png", parentDir, modelName, cond)
  ggsave(fileName, width = 6, height = 6)
}






