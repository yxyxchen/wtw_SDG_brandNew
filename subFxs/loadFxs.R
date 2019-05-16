
# loadData.R
# varying-magnitude WTW

# each subject has:
#  header file (*hdr.txt)
#  main data file (ending both in .mat and .txt)
#  empty *keytimes.txt file, a vestige of the effort key-pressing version. 


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
    thisCbal = hdrData$Cbal[sIdx]
    thisCond = hdrData$condition[sIdx]
    thisStress = hdrData$stress[sIdx]
    # loop over blocks
    junk = vector(nBlocks, mode = 'list')
    for (bkIdx in 1 : nBlocks){
      thisFile = list.files(path=dataDir, pattern=(sprintf('wtw_stress_SDG%d_bk%d_1.txt',thisID, bkIdx)))
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


# l

loadExpPara = function(paras, dirName){
  # number of paras 
  nE = length(paras) + 1
  # number of files
  fileNames = list.files(path= dirName, pattern=("*_summary.txt"))
  library("gtools")
  fileNames = mixedsort(sort(fileNames))
  n = length(fileNames) 
  sprintf("load %d files", n)
  
  # load blockData   
  load("genData/expDataAnalysis/blockData.RData")
  blockData = blockData[blockData$blockNum == 1,]
  
  # initialize the outout variable 
  expPara = matrix(NA, n, nE * 4)
  idList = vector(length = n)
  # loop over files
  for(i in 1 : n){
    fileName = fileNames[[i]]
    address = sprintf("%s/%s", dirName, fileName)
    junk = read.csv(address, header = F)
    idList[i] = ifelse(is.na(as.numeric(substr(fileName, 2,4))),
                       ifelse(is.na(as.numeric(substr(fileName, 2,3))),  as.numeric(substr(fileName, 2,2)),
                              as.numeric(substr(fileName, 2,3))),
                       as.numeric(substr(fileName, 2,4)))
    # delete the lp__ in the old version
    if(nrow(junk) > nE){
      junk = junk[1:nE,]
    }
    expPara[i, 1:nE] = junk[,1]
    expPara[i, (nE + 1) : (2 * nE)] = junk[,2]
    expPara[i, (2*nE + 1) : (3 * nE)] = junk[,9]
    expPara[i, (3 * nE + 1) : (4 * nE)] = junk[,10]
  }
  # transfer expPara to data.frame
  expPara = data.frame(expPara)
  junk = c(paras, "LL_all")
  colnames(expPara) = c(junk, paste0(junk, "SD"), paste0(junk, "Effe"), paste0(junk, "Rhat"))
  # add id, condition, stress
  condition = blockData$condition[blockData$id %in% idList]
  stress = blockData$stress[blockData$id %in% idList]
  expPara[['id']] = idList  # needed for left_join
  expPara[["condition"]] = condition
  expPara[['stress']] = stress
  return(expPara)
}

getMode = function(x){
  ob = density(x)
  value = ob$x
  dense = ob$y
  return(value[which.max(dense)])
}

loadExpParaExtra = function(modelName, paras){
  nE = length(paras) + 2
  load("genData/expDataAnalysis/blockData.RData")
  idList = unique(blockData$id) 
  n = length(idList)
  nE = length(paras) + 1
  expParaMode = matrix(NA, n, nE)
  expParaMedian = matrix(NA, n, nE)
  for(i in 1 : n){
    ID = idList[i]
    fileName = sprintf("genData/expModelFitting/%s/s%d.txt", modelName, ID)
    junk = read.csv(fileName, header = F)
    
    expParaMode[i, ] = apply(junk[,1:nE], MARGIN = 2, getMode)
    expParaMedian[i, ] = apply(junk[,1:nE], MARGIN = 2, median)
  }
  expParaMode = data.frame(expParaMode)
  expParaMedian = data.frame(expParaMedian)
  
  junk = c(paras, "log_lik")
  colnames(expParaMode) =  paste0(junk, "Mode")
  expParaMode$id = idList  # needed for left_join
  colnames(expParaMedian) =  paste0(junk, "Median")
  expParaMedian$id = idList  # needed for left_join
  
  outputs = list("expParaMedian" = expParaMedian, "expParaMode" = expParaMode)
  return(outputs)
}

# for each sub
loadExpParaMedian = function(modelName, paras, id){
  nE = length(paras) + 2
  load("genData/expDataAnalysis/blockData.RData")
  nE = length(paras) + 1
  expParaMedian = vector(length = nE)
    fileName = sprintf("genData/expModelFitting/%s/s%d.txt", modelName, id)
    junk = read.csv(fileName, header = F)
    expParaMedian  = apply(junk[,1:nE], MARGIN = 2, median)
  return(expParaMedian)
}


loadSimPara_ = function(modelName, paras, cond){
  nE = length(paras) + 2
  load(sprintf("genData/simulation/%s/simParas.RData", modelName))
  n = nComb
  nE = length(paras)
  expPara = array(NA, dim = c(nComb, nRep, nE*5))
  for(i in 1 : n){
    expPara_ = matrix(NA, nRep, nE*5)
    for(r in 1 : nRep){
      fileName = sprintf("genData/simModelFitting/%s/%s_s%d_r%d_summary.txt",
                         modelName, cond, i, r)
      junk = read.csv(fileName, header = F)
      expPara_[r, 1:nE] = junk[1:nE,1]
      expPara_[r, (nE + 1) : (2 * nE)] = junk[1:nE,2]
      expPara_[r, (2*nE + 1) : (3 * nE)] = junk[1:nE,9]
      expPara_[r, (3 * nE + 1) : (4 * nE)] = junk[1:nE,10]
      expPara_[r, (4 * nE + 1) : (5 * nE)] = paraComb[i,1:4]
    }
    expPara[i, ,] = expPara_
  }
  return(expPara)
}
