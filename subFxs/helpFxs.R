library("stringr")
getParas = function(modelName){
  if(str_detect(modelName, "curiosityTrial")) paras = c("phi", "tau", "gamma")
  else if(str_detect(modelName, "curiosityTrialR")) paras = c("phi", "tau", "phiR")
  else if(str_detect(modelName, "heuristicRL")) paras = c("phi", "tau", "threshd", "ini")
  else return("wrong model name")
  return(paras)
}

getUseID = function(blockData, expPara, paras){
  idList = unique(blockData$id)
  n = length(idList)
  RhatCols = which(str_detect(colnames(expPara), "hat"))[1 : length(paras)]
  EffeCols = which(str_detect(colnames(expPara), "Effe"))[1 : length(paras)]
  useID = idList[apply(expPara[,RhatCols] < 1.1, MARGIN = 1, sum) == length(paras) & 
                   apply(expPara[,EffeCols] >100, MARGIN = 1, sum) == length(paras)]
  return(useID)
}
  
getParaComb = function(paraTable){
  paraNames = names(paraTable)
  nPara = length(paraTable)
  nValue = unlist(lapply(1 : nPara, function(i) length(paraTable[[i]])))
  
  output = matrix(NA, nrow = prod(nValue), ncol = nPara)
  for(pIdx in 1 : nPara){
    eachRep = ifelse(pIdx >= nPara, 1, prod(nValue[(pIdx+1) : nPara])) # repetition number for each element in the seq
    seqRep = ifelse(pIdx <= 1, 1, prod(nValue[1 : (pIdx - 1)])) # repetition number for the seq
    output[,pIdx]= rep(rep(paraTable[[pIdx]], each = eachRep), seqRep)
  }
  colnames(output) = paraNames
  return(output)
}



