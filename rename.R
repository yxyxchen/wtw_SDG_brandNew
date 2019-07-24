

modelNames = c("BL", "QL1", "QL2", "RL1", "RL2")
nModel = length(modelNames)
for(i in 1 : nModel){
  modelName = modelNames[i]
  dir = sprintf("genData/expModelFitting/%sdb", modelName)
  files = list.files(path = dir)
  nFile = length(files)
  for(f in 1 : nFile){
    file = files[[f]]
    address = sprintf("%s/%s", dir, file)
    idIndexs = str_locate(file, "s[0-9]+")
    oldID = as.double(substr(file, idIndexs[1]+1, idIndexs[2]))
    file.rename(from=address,to=sub(pattern= "s[0-9]+",
                                    replacement= sprintf("s%s", hdrData$ID[oldID]),
                                    address))
  } 
}

