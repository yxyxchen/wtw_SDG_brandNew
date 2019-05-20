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

# input 
dataType = "sess"
modelName = "curiosityTrialSp"

# create output directories
dir.create("figures/expParaAnalysis")
saveDir = sprintf("figures/expParaAnalysis/%s", modelName)
dir.create(saveDir)

# load blockdata data
if(dataType == "block"){
  load("genData/expDataAnalysis/blockData.RData")
  load("genData/expDataAnalysis/kmOnGridBlock.RData")
  summaryData = blockData[blockData$blockNum == 1, ] # so summaryData only have something relevant
}else{
  load("genData/expDataAnalysis/sessionData.RData")
  load("genData/expDataAnalysis/kmOnGridSess.RData")
  summaryData = sessionData
}

# maybe I shoud truncate summaryData to make them the same as balabala 

# load expPara
paras = getParas(modelName)
nPara = length(paras)
parentDir = ifelse(dataType == "block", "genData/expModelFitting", "genData/expModelFittingSub")
dirName = sprintf("%s/%s",parentDir, modelName)
tempt = loadExpPara(paras, dirName)
useID = getUseID(tempt, paras)
expPara = merge(x=tempt[,c(paras, "id")],y=summaryData, by="id",all.x=TRUE)


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

# calculate trait-para correlations
traitParaCorr = vector(mode = "list", length = length(paras) * length(traits))
corrNo = matrix(1:(nPara * nTrait), nrow = nPara, ncol = nTrait)
for(pIdx in 1 : length(paras)){
  para = paras[pIdx]
  paraColor = paraColors[pIdx]
  for(trIdx in 1 : length(traits)){
    trait = traits[trIdx]
    traitName = traitNames[trIdx]
    # combine persondatlaity, expPara and condtion. useID only
    input = data.frame(personality[summaryData$id %in% useID,trait], expPara[expPara$id %in% useID,para],
                       summaryData$condition[summaryData$id %in% useID])
    traitParaCorr[[corrNo[pIdx, trIdx]]]= getCorrelation(input)
    # plot
    p = plotCorrelation(input, paraColor, T) 
    p + ylab(capitalize(para)) + xlab(traitName) + saveTheme
    fileName = sprintf("%s/%s/%s_%s.png", parentDir, modelName, para, traitName)
    ggsave(fileName, width = 6, height = 3)
  }
}
dimNames = list(paras, traitNames)
rhoTable = lapply(1:2, function(j) matrix(sapply(1: (nTrait * nPara), function(i) traitParaCorr[[i]]$rhos[j]),
                                          nrow = nPara, dimnames = dimNames))
pTable = lapply(1:2, function(j) matrix(sapply(1: (nTrait * nPara), function(i) traitParaCorr[[i]]$ps[j]),
                                        nrow = nPara, dimnames = dimNames))
# plot trait analysis
library("corrplot")
for(i in 1 : 2){
  cond = conditions[i]
  fileName = sprintf("traitPara%s_%s.png", cond, dataType)
  png(fileName)
  corrplot(rhoTable[[i]], 
           p.mat = pTable[[i]], 
           is.corr = T, 
           method = "color",
           # insig = "label_sig",
           tl.col = "black", tl.srt = 15, tl.cex = 1.5) 
  mtext("Parameter", side = 2,  cex = 2)
  mtext("Trait", side = 3, cex = 2)
  dev.off()
}


# look at trait linearly
expPara = expPara[!is.na(expPara$Anxiety),]
for(pIdx in 1 : nPara){
  para = paras[pIdx]
  paraColor = paraColors[pIdx]
  for(cIdx in 1 : 2){
    cond = conditions[cIdx]
    expPara[expPara$condition == cond,c(para, traitNames)] %>%
      gather(-c(para), key = "trait", value = "value") %>% 
      ggplot(aes_string(x = "value", y = para)) + geom_point() +
      facet_grid(~ trait, scales = "free") +
      theme_linedraw(base_size = 13) + ylab(capitalize(para)) + xlab("Trait value")
    parentDir = ifelse(dataType == "sess", "figures/expParaAnalysisSub", "figures/expParaAnalysis")
    fileName = sprintf("%s/%s/lm_%s_%s.png", parentDir, modelName, cond, para, cond)
    ggsave(fileName, width = 8, height = 3)
  }
}

trait = traits[2]
para = "gamma"
input = data.frame(x = personality[summaryData$id %in% useID,trait], y = expPara[expPara$id %in% useID,para],
           cond= summaryData$condition[summaryData$id %in% useID])
input = input[!is.na(input$x),]
ggplot(input, aes(x, y)) + geom_point() + facet_grid(~cond)

gapSize = 5
input$xGroup = cut2(input$x, cuts = seq(floor(min(input$x) / gapSize) * gapSize,
                                        ceiling(max(input$x) / gapSize) * gapSize, by = gapSize))
sumInput = group_by(input, cond, xGroup) %>% dplyr::summarise(meanY = mean(y))
ggplot(sumInput, aes(xGroup, meanY)) + geom_point() +facet_grid(~cond) + xlab("Impulsive")+ylab("Gamma")

# plot hist 
expPara$condition = summaryData$condition[summaryData$id %in% expPara$id]
for(i in 1 : length(paras)){
  para = paras[i]
  paraColor = paraColors[[i]]
  p = ggplot(expPara[expPara$id %in% useID,], aes_string(para)) + geom_histogram(bins = 6, fill = paraColor) + facet_grid(~condition)+
    saveTheme + ylab("Count") + xlab(capitalize(para))
  parentDir = ifelse(dataType == "sess", "figures/expParaAnalysisSub", "figures/expParaAnalysis")
  dir.create(parentDir)
  fileName = sprintf("%s/hist_%s.pdf", parentDir, para)
  ggsave(fileName, width = 6, height = 3)
}

