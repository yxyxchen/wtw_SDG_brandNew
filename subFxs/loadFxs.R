
loadAllData = function() {
  ##### step 1: organize hdrData
  # load summary data
  dataDir = 'data'
  fileName = sprintf('%s/SDGdataset.csv', dataDir)
  summaryData= read.csv(fileName)
  # adjust
  summaryData$stress = ifelse( summaryData$Condition == 'stress', 'stress', 'no_stress')
  summaryData$Task..1...unif..2...gp. = ifelse(summaryData$Task..1...unif..2...gp. == 1, 'HP', 'LP')
  
  # exclude AUC and totalEarnings from hdrData
  hdrData = summaryData[,-(4:17)]
  colnames(hdrData) = c('ID', 'stress', 'condition', 'cbal', 'perceivedStress',
                        'traitAnxiety', 'Gender', 'BDI', 'posAffect1', 'posAffect2', 
                        'negAffect1', 'negAffect2', 'uncertainty', 'delay',
                        'impulsive', 'postUnpleasant')
  hdrData$ID = as.character(hdrData$ID)
  # count number of subjects 
  nSubjects = nrow(summaryData)
  nBlocks = 3

  ##### step 2: load trialData
  # define data column names
  colnames = c('blockNum', 'trialNum', 'trialStartTime', 'nKeyPresses', 'scheduledWait',
               'rewardTime', 'timeWaited', 'sellTime', 'trialEarnings','totalEarnings')
  # 'timeWaited' is from trial onset to keypress (includes RT)
  # 'sellTime' is the keypress time relative to block onset.
  
  # initialize
  trialData = list()

  # loop over individual subjects
  for (sIdx in 1:nSubjects) {
    thisID = hdrData$ID[sIdx]
    thisCbal = hdrData$Cbal[hdrData$ID == thisID]
    thisCond = hdrData$condition[hdrData$ID == thisID]
    thisStress = hdrData$stress[hdrData$ID == thisID]
    # loop over blocks
    junk = vector(nBlocks, mode = 'list')
    for (bkIdx in 1 : nBlocks){
      thisFile = list.files(path=dataDir, pattern=(sprintf('wtw_stress_SDG%s_bk%d_1.txt',thisID, bkIdx)))
      if (length(thisFile) != 1) {
        cat('Could not identify a single data file for subject',thisID,' block', bkIdx, '\n')
        browser()
      }
      junk[[bkIdx]] = read.csv(file.path(dataDir,thisFile), header=FALSE, col.names=colnames)
      
      # adjust blockNum
      junk[[bkIdx]]$blockNum = junk[[bkIdx]]$blockNum * bkIdx
    }
    
    d = rbind(junk[[1]], junk[[2]], junk[[3]])
    # add info from hdrData
    d$condition = rep(thisCond, nrow(d))
    d$stress = rep(thisStress, nrow(d))
    
    
    # add to the list of all subjects' data
    trialData[[thisID]] = d
  } # end of loop over subjects
  
  # return the 2 data frames in a named list
  outputData = list(hdrData=hdrData, trialData=trialData)
  return(outputData)
  
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
  expPara = data.frame(expPara)
  junk = c(paraNames, "LL_all")
  colnames(expPara) = c(junk, paste0(junk, "SD"), paste0(junk, "Effe"), paste0(junk, "Rhat"))
  expPara$id = factor(idList, levels = levels(hdrData$ID)) # ensure the levels are consistent, usually not that problematic though
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
  expPara$id = factor(idList, levels = idList)
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

