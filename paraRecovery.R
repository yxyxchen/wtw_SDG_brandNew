load("expParas.RData")
library("tidyverse")
source("subFxs/plotThemes.R")
source("subFxs/loadFxs.R") # load blockData and expPara
source("subFxs/helpFxs.R") # getparaNames
library(latex2exp)
library(ggpubr)
theme_set(theme_pubr())

# model Name
modelName = "QL2"
paraNames = getParaNames(modelName)
nPara = length(paraNames)
paraLabels = c("$\\alpha_r$", "$\\alpha_u$", "$\\tau$", "$\\gamma$", "$\\eta$")

# output directories
saveDir = "figures/paraRecovery"
dir.create(saveDir)

# load expPara
paraNames = getParaNames(modelName)
nPara = length(paraNames)
dirName = sprintf("%s/%s","genData/expModelFit", modelName)
expParas = loadExpPara(paraNames, dirName)
passCheck = checkFit(paraNames, expParas)


# load repPara
dirName = sprintf("%s/%s/%s","genData/simModelFit", modelName, modelName)
repParas = loadExpPara(paraNames, dirName)


# filter
select = checkFit(paraNames, repParas) & checkFit(paraNames, expParas)
expParas = expParas[select,]
repParas = repParas[select,]

# plot
ps = list(length = nPara)
lows = c(0, 0, 0, 0.7, 0)
ups = c(0.3, 0.3, 22, 1, 7)
for(i in 1 : nPara){
  paraName = paraNames[i]
  lims = c(lows[i], ups[i])
  lims2 = lims
  lims2[1] = lims[1] - 0.1
  lims2[2] = lims[2] + 0.1
  paraLabel = paraLabels[i]
  ps[[i]] = data.frame(
    expPara = unlist(expParas[paraName]),
    repPara = unlist(repParas[paraName])) %>%
    ggplot(aes(expPara, repPara)) + geom_point() + myTheme + 
    scale_x_continuous(breaks = lims, labels = lims, limits = lims2) +
    scale_y_continuous(breaks = lims, labels = lims, limits = lims2) +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(TeX(paraLabels[i])) +
    ylab(TeX(sprintf("$\\hat{%s}$", substr(paraLabels[i],2, nchar(paraLabels[i]) - 1)))) +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
    geom_abline(slope = 1, color = "grey", linetype = 2)
}
ggarrange(ps[[1]], ps[[2]], ps[[3]], ps[[4]], ps[[5]],
                    ncol = nPara, nrow = 1)
ggsave(sprintf("figures/paraRecovery/%s_%s.png", modelName, modelName), width = 8, height = 2)
