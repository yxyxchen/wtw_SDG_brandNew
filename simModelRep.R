# libraries and scripts
library("ggplot2")
library("dplyr")
library("tidyr")
source("subFxs/helpFxs.R")
source("subFxs/loadFxs.R")
source("subFxs/plotThemes.R")

# make output dir
dir.create("figures")
dir.create("figures/simModelRep")

# load model names
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
ids = hdrData$id[hdrData$stress == "no_stress"]                 # column of subject IDs
nSub = length(ids) 

# check fit
modelName = "QL2"
paraNames = getParaNames(modelName)
nPara = length(paraNames)
expPara = loadExpPara(paraNames, sprintf("genData/expModelFit/%s", modelName))
simPara = loadExpPara(paraNames, sprintf("genData/simModelFit/%s/%s", modelName, modelName)) 
passCheckIds = simPara$id[checkFit(paraNames, simPara)]
nPassCheck = length(passCheckIds)

tempt = data.frame(rbind(expPara[,c(1:nPara, nPara*4 + 6)], simPara[,c(1:nPara,nPara*4 + 6)])) %>%
  filter(id %in% passCheckIds) %>% 
  mutate(source = rep(c("exp", "sim"), each = nPassCheck )) %>%
  gather(key = "paraName", value = "paraValue", -source, -id) %>%
  mutate(paraName = factor(paraName, levels = paraNames))


# reorganize to get paraData
paraData = data.frame(paraName = tempt$paraName[tempt$source == "exp"],
           expValue = tempt$paraValue[tempt$source == "exp"],
           simValue = tempt$paraValue[tempt$source == "sim"]) 

# calculate the wilcox correlations 
rankCorResults = lapply(1 : nPara, FUN = function(i) cor.test(
  paraData$expValue[paraData$paraName == paraNames[i]],
  paraData$simValue[paraData$paraName == paraNames[i]],
  method = "kendall"
))

# summary statistics for paraData
paraSumStats = paraData %>%
  gather(key = "source", value = "value", -paraName) %>%
  group_by(source, paraName) %>%
  summarise(upperq = quantile(value)[3], lowerq = quantile(value)[1]) %>%
  mutate(IQR = upperq - lowerq, 
         upperLimit = upperq + IQR * 1.5,
         lowerLimit = lowerq - IQR * 1.5)


# we use the 1.5 x the IQR for the larger of the two values 
# (generating or recovered parameter estimates) as the criteria 
# for exclude outliers
limits = data.frame(
  upperLimits = pmax(paraSumStats$upperLimit[paraSumStats$source == "simValue"],
                     paraSumStats$upperLimit[paraSumStats$source == "expValue"]),
  lowerLimits = pmin(paraSumStats$lowerLimit[paraSumStats$source == "simValue"],
                     paraSumStats$lowerLimit[paraSumStats$source == "expValue"])
)
rownames(limits) = unique(paraSumStats$paraName)
paraDataClean = 
  paraData %>% 
  filter(expValue <= limits[paraName, "upperLimits"],
         expValue >= limits[paraName, "lowerLimits"],
         simValue <= limits[paraName, "upperLimits"],
         expValue >= limits[paraName, "lowerLimits"])

# eliminate outliers and plot
# we only use 57 here 
for(i in 1 : nPara){
  paraName = paraNames[i]
  paraDataClean[paraDataClean$paraName == paraName,] %>%
    ggplot(aes(expValue, simValue)) +
    geom_point(color = themeColor, fill = "pink", shape = 21,
               size = 4) +
    myTheme
  ggsave(sprintf("figures/simModelRep/%s.eps", paraName),
         width = 4, height = 4)
}



  
  
  
  
  
