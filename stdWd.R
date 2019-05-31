traitParaCorr = vector(mode = "list", length = nPara)
for(pIdx in 1 : length(paras)){
  para = paras[pIdx]
  # plotParaAUC(expPara, para, blockData, useID)
  input = data.frame(x = expPara[[para]], y = expPara$stdWd, z = expPara$AUC,
                     cond = expPara$condition)
  input = input[expPara$id %in% useID,]
  traitParaCorr[[pIdx]] = getPartCorrelation(input)
}
rhoTable = sapply(1:2, function(j) sapply(1: nPara, function(i) traitParaCorr[[i]]$rhos[j]))
rownames(rhoTable) = c("LR", "LP", "Tau", "Gamma", "P")
colnames(rhoTable) = c("WTW_SD_HP", "WTW_SD_LP")
pTable = sapply(1:2, function(j) sapply(1: nPara, function(i) traitParaCorr[[i]]$ps[j]))

library("corrplot")
col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061")))

parentDir = ifelse(dataType == "sess", "figures/expParaAnalysisSub", "figures/expParaAnalysis")

fileName = sprintf("%s/%s/stdWdPara_%s.png", parentDir, modelName, dataType)
png(fileName)
corrplot(t(rhoTable), 
         p.mat = t(pTable), 
         is.corr = T, 
         method = "color",
         insig = "p-value",
         tl.col = "black", tl.srt = 0, tl.cex = 1.5,
         col = col2(50)) 
dev.off()



####
traitParaCorr = vector(mode = "list", length = nPara)
for(pIdx in 1 : length(paras)){
  para = paras[pIdx]
  # plotParaAUC(expPara, para, blockData, useID)
  input = data.frame(x = expPara[[para]], y = expPara$AUC, z = expPara$stdWd,
                     cond = expPara$condition)
  input = input[expPara$id %in% useID,]
  traitParaCorr[[pIdx]] = getPartCorrelation(input)
}
rhoTable = sapply(1:2, function(j) sapply(1: nPara, function(i) traitParaCorr[[i]]$rhos[j]))
rownames(rhoTable) = c("LR", "LP", "Tau", "Gamma", "P")
colnames(rhoTable) = c("WTW_AVE_HP", "WTW_AVE_LP")
pTable = sapply(1:2, function(j) sapply(1: nPara, function(i) traitParaCorr[[i]]$ps[j]))

library("corrplot")
col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061")))

parentDir = ifelse(dataType == "sess", "figures/expParaAnalysisSub", "figures/expParaAnalysis")

fileName = sprintf("%s/%s/AUCPara_%s.png", parentDir, modelName, dataType)
png(fileName)
corrplot(t(rhoTable), 
         p.mat = t(pTable), 
         is.corr = T, 
         method = "color",
         insig = "p-value",
         tl.col = "black", tl.srt = 0, tl.cex = 1.5,
         col = col2(50)) 
dev.off()
