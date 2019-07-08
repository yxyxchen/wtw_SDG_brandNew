
likRatioTest = function(modelName1, modelName2, df, group = "all"){
  paras1 = getParas(modelName1)
  expPara1 = loadExpPara(paras1,
                         sprintf("genData/expModelFitting/%s", modelName1))
  useID1 = getUseID(expPara1, paras1)
  paras2 = getParas(modelName2)
  expPara2 = loadExpPara(paras2,
                         sprintf("genData/expModelFitting/%s", modelName2))
  useID2 = getUseID(expPara2, paras2)
  if(group == "HP" || group == "LP"){
    useID = allIDs[allIDs %in% useID1 & allIDs %in% useID2 & hdrData$condition == group]
  }
  
  sumLogEvi1 = filter(expPara1, id %in% useID) %>% summarise(sum(LL_all))
  sumLogEvi2 = filter(expPara2, id %in% useID) %>% summarise(sum(LL_all))
  
  delta = -2 * as.double(sumLogEvi1 - sumLogEvi2)
  p = pchisq(delta, df * length(useID), lower.tail=FALSE)
  return(p)
}