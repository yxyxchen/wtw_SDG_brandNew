library("stringr")
getParas = function(modelName){
  if(modelName == "curiosityTrial") paras = c("phi", "tau", "gamma")
  else if(modelName == "curiosityTrialR") paras = c("phi", "tau", "phiR")
  else if(modelName == "curiosityTrialSp") paras = c("phi", "tau", "gamma", "zeroPoint")
  else if(modelName == "curiosityTrialRSp") paras = c("phi", "tau", "phiR", "zeroPoint")
  else if(modelName == "fullModel") paras = c("phi", "tau", "gamma", "QwaitIni")
  else if(modelName == "full_model") paras = c("phi", "tau", "gamma", "QwaitIni")
  else return("wrong model name")
  return(paras)
}

getUseID = function(expPara, paras){
  idList = expPara$id
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



