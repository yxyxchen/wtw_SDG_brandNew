loadAllData = function() {
  # loads hdrData and trialData
  
  # outputs:
  # hdrData: exp info for each participant, like ID, condition, questionnaire measurements
  # trialData: a length = nSub list, each element containing trial-wise data for each partcipant. 
  
  # each elemant of trialData is formatted similarly to this example:
  # blockNum : [130x1 int]
  # trialNum : [130x1 int]
  # trialStartTime : [130x1 num] # when participant start a trial
  # sellTime : [130x1 num] # when participants sell a token
  # nKeyPresses : [130x1 int]
  # scheduledWait : [130x1 num] # delay durations for rewards
  # rewardTime : [130x1 num] # actual delay durations (constrainted by the screen update freq), NA for non-rewarded trials
  # timeWaited : [130x1 num] # persistence durations, namely sellTime - trialStartTime
  # trialEarnings : [130x1 int] trial-wise payments, either 10 or 0
  # totalEarnings : [130x1 int] cumulative payments
  
  # load experiment paras
  load('expParas.RData')
  
  # load hdrData
  hdrData = read.csv(file.path('data', 'hdrData.csv'), comment = "#")
  hdrData$stress = ifelse( hdrData$stress == 'stress', 'stress', 'no_stress')
  hdrData$condition = ifelse(hdrData$condition == 1, 'HP', 'LP')
  hdrData$id = as.character(hdrData$id)
  
  # load trialData
  trialData = list()
  nSub = nrow(hdrData)
  trialDataNames = c('blockNum', 'trialNum', 'trialStartTime', 'nKeyPresses', 'scheduledWait',
                      'rewardTime', 'timeWaited', 'sellTime', 'trialEarnings','totalEarnings')
  # loop over subjects
  for (sIdx in 1:nSub) {
    id = hdrData$id[sIdx]
    condition = hdrData$condition[hdrData$id == id]
    stress = hdrData$stress[hdrData$id == id]
    thisTrialData = vector(nBlock, mode = 'list')
    # loop over blocks
    for (bkIdx in 1 : nBlock){
      blockFile = list.files(path="data", pattern=(sprintf('wtw_stress_SDG%s_bk%d_1.txt',id, bkIdx)))
      if (length(blockFile) != 1) {
        cat('Could not identify a single data file for subject',thisID,' block', bkIdx, '\n')
        browser()
      }
      thisTrialData[[bkIdx]] = read.csv(file.path("data", blockFile), header = F, col.names = trialDataNames)
      thisTrialData[[bkIdx]]$blockNum = rep(bkIdx, nrow(thisTrialData[[bkIdx]]))
      thisTrialData[[bkIdx]]$condition = rep(condition, nrow(thisTrialData[[bkIdx]]))
    }
    trialData[[id]] = bind_rows(thisTrialData)
  } # end of loop over subjects
  # return outputs
  outputs = list(hdrData=hdrData, trialData=trialData)
  return(outputs)
} 




loadExpPara = function(paraNames, dirName){
  # number of paraNames 
  nE = length(paraNames) + 1
  # number of files
  fileNames = list.files(path= dirName, pattern=("*_summary.txt"))
  library("gtools")
  fileNames = mixedsort(sort(fileNames))
  n = length(fileNames) 
  sprintf("load %d files", n)
  
  # initialize the outout variable 
  expPara = matrix(NA, n, nE * 4)
  idList = vector(length = n)
  # loop over files
  for(i in 1 : n){
    fileName = fileNames[[i]]
    address = sprintf("%s/%s", dirName, fileName)
    junk = read.csv(address, header = F)
    idIndexs = str_locate(fileName, "s[0-9]+")
    idList[i] = substr(fileName, idIndexs[1]+1, idIndexs[2])
    # delete the lp__ in the old version
    if(nrow(junk) > nE){
      junk = junk[1:nE,]
    }
    expPara[i, 1:nE] = junk[,1]
    expPara[i, (nE + 1) : (2 * nE)] = junk[,3]
    expPara[i, (2*nE + 1) : (3 * nE)] = junk[,9]
    expPara[i, (3 * nE + 1) : (4 * nE)] = junk[,10]
  }
  # transfer expPara to data.frame
  expPara = data.frame(expPara, stringsAsFactors = F)
  junk = c(paraNames, "LL_all")
  colnames(expPara) = c(junk, paste0(junk, "SD"), paste0(junk, "Effe"), paste0(junk, "Rhat"))
  expPara$id = idList # ensure the levels are consistent, usually not that problematic though
  return(expPara)
}

# I also need to load 2.5% and 97.5%
loadCVPara = function(paraNames, dirName, pattern){
  # number of paraNames 
  nE = length(paraNames) + 1
  # number of files
  fileNames = list.files(path= dirName, pattern= pattern)
  library("gtools")
  fileNames = mixedsort(sort(fileNames))
  n = length(fileNames) 
  sprintf("load %d files", n)
  
  # initialize the outout variable 
  expPara = matrix(NA, n, nE * 6)
  idList = vector(length = n)
  # loop over files
  for(i in 1 : n){
    fileName = fileNames[[i]]
    address = sprintf("%s/%s", dirName, fileName)
    junk = read.csv(address, header = F)
    sIndexs = str_locate(fileName, "s[0-9]+")
    s = substr(fileName, sIndexs[1]+1, sIndexs[2])
    fIndexs = str_locate(fileName, "f[0-9]+")
    f = substr(fileName, fIndexs[1]+1, fIndexs[2])
    idList[i] = sprintf("s%s_f%s", s, f)
    # delete the lp__ in the old version
    if(nrow(junk) > nE){
      junk = junk[1:nE,]
    }
    expPara[i, 1:nE] = junk[,1]
    expPara[i, (nE + 1) : (2 * nE)] = junk[,3]
    expPara[i, (2*nE + 1) : (3 * nE)] = junk[,9]
    expPara[i, (3 * nE + 1) : (4 * nE)] = junk[,10]
    expPara[i, (4*nE + 1) : (5 * nE)] = junk[,4]
    expPara[i, (5 * nE + 1) : (6 * nE)] = junk[,8]
  }
  # transfer expPara to data.frame
  expPara = data.frame(expPara)
  junk = c(paraNames, "LL_all")
  colnames(expPara) = c(junk, paste0(junk, "SD"), paste0(junk, "Effe"), paste0(junk, "Rhat"),
                        paste0(junk, "2.5"),paste0(junk, "97.5"))
  expPara$id = idList
  return(expPara)
}

loadSimPara = function(paraNames, dirName){
  # number of paraNames 
  nE = length(paraNames) + 1
  # number of files
  fileNames = list.files(path= dirName, pattern=("*_summary.txt"))
  library("gtools")
  fileNames = mixedsort(sort(fileNames))
  n = length(fileNames) 
  sprintf("load %d files", n)
  
  # extract seqIdxs and cbIdxs from filenames
  seqIdxs = as.double(sapply(1:n, function(i){
    tempt = regexpr("_r[0-9]*", fileNames[i]) 
    match.length = attr(tempt, "match.length")
    start = tempt[[1]] + 2
    ending = tempt[[1]] + match.length - 1
    substring(fileNames[i], start, ending)
  })) + 1
  nSeq = length(unique(seqIdxs))
  
  cbIdxs =  as.double(sapply(1:n, function(i){
    tempt = regexpr("_s[0-9]*", fileNames[i]) 
    match.length = attr(tempt, "match.length")
    start = tempt[[1]] + 2
    ending = tempt[[1]] + match.length - 1
    substring(fileNames[i], start, ending)
  }))
  nComb = length(unique(cbIdxs))
  
  # initialize the outout variable 
  junk = c(paraNames, "LL_all")
  varNames = c(junk, paste0(junk, "SD"), paste0(junk, "Effe"), paste0(junk, "Rhat"))
  simParaHP = array(dim = c(nE * 4, nSeq, nComb),
                    dimnames = list(varNames, 1 : nSeq, 1 : nComb))
  simParaLP = array(dim = c(nE * 4, nSeq, nComb), dimnames = list(varNames, 1 : nSeq, 1 : nComb))
  # loop over files
  for(i in 1 : n){
    # read the file
    fileName = fileNames[[i]]
    address = sprintf("%s/%s", dirName, fileName)
    junk = read.csv(address, header = F)
    
    # determine the condition 
    if(grepl("HP", fileName)){
      simParaHP[, seqIdxs[i], cbIdxs[i]] = unlist(junk[,c(1,2,9,10)])
    }else{
      simParaLP[, seqIdxs[i], cbIdxs[i]] = unlist(junk[,c(1,2,9,10)])
    }
  }
  simPara = list(HP = simParaHP, LP =simParaLP)
  # transfer expPara to data.frame
  return(simPara)
}

