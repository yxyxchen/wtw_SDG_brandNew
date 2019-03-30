library("stringr")
getParas = function(modelName){
  if(str_detect(modelName, "full_model")) paras = c("phi", "tau", "gamma", "QwaitIni")
  else if(str_detect(modelName, "QwaitIni")) paras = c("phi", "tau", "gamma")
  else if(str_detect(modelName, "reduce_one_phi")) paras = c("tau", "gamma", "QwaitIni")
  else if(str_detect(modelName, "reduce_one_gamma")) paras = c("phi", "tau", "QwaitIni")
  else if(modelName == "baseline") paras = c("waitRate")
  else if(str_detect(modelName, "tau")) paras = c("phi", "gamma", "QwaitIni")
  else if(str_detect(modelName, "no_learning_2")) paras = c("Tau")
  else if(str_detect(modelName, "R_learning2")) paras = c("phi1", "phi2", "tau")
  else if(str_detect(modelName, "monteRatio")) paras = c("phi", "tau", "gamma", "ratio")
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
  nPara = ncol(paraTable)
  nValue = nrow(paraTable)
  
  output = matrix(NA, nrow = nValue ^ nPara, ncol = nPara)
  for(i in 1 : nPara){
    output[,i]= rep(rep(paraTable[,i], each = nValue ^ (nPara - i)), nValue ^ (i- 1))
  }
  colnames(output) = paraNames
  return(output)
}



