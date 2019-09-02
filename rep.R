modelName = "RL2"
paraNames = getParaNames(modelName)
expPara = loadExpPara(paraNames, sprintf("genData/expModelFit/%s", modelName))


allData = loadAllData()
