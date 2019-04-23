# this script contains helper analysis functions 

# check the distribution of scheduled delays
# ...as measured in number of key presses (for the instrumental version of the task)
scheduledDelays <- function(blockData,label) {
  cat(sprintf('Scheduled delays for %s\n',blockLabel))
  bkDelays = blockData$scheduledWait
  print(summary(bkDelays))
  # empirical cumulative distribution of scheduled delays
  fn <- ecdf(blockData$scheduledWait)
  plot(fn, main = sprintf('Scheduled delays: %s',blockLabel), xlab='Scheduled delay (s)',
       ylab='Cumulative proportion', xlim=c(0,30))
  # autocorrelation function
  # acfOutput <- acf(bkDelays, lag.max=20, main = sprintf('Scheduled delays: %s',blockLabel))
}


# plot trialwise responses in detail
trialPlots <- function(blockData,label) {
  # vectors to be plotted
  rwdIdx = blockData$trialEarnings > loseValue
  quitIdx = blockData$trialEarnings <= loseValue
  rwdTrialNo = blockData$trialNum[rwdIdx]
  quitTrialNo = blockData$trialNum[quitIdx]
  rwdSchedDelay = blockData$scheduledWait[rwdIdx]
  quitSchedDelay = blockData$scheduledWait[quitIdx]
  waitDuration = blockData$timeWaited
  quitTime = waitDuration[quitIdx]
  # other parameters
  nTrials = length(blockData$trialEarnings)
  # make the plot and add series
  plotData = data.frame(trialNum = c(rwdTrialNo, quitTrialNo, quitTrialNo),
                        trialDuration = c(rwdSchedDelay, quitTime, quitSchedDelay),
                        condition = rep(c('reward', 'quit', 'quitSchedule'), time = 
                                        c(length(rwdTrialNo), length(quitTrialNo),
                                          length(quitTrialNo))))
  plotData$condition = factor( plotData$condition, levels = c('reward', 'quit', 'quitSchedule'))
  p = ggplot(plotData, aes(trialNum, trialDuration, color = condition)) + geom_point() +
  geom_line(data = plotData[plotData$condition != 'quitSchedule',],
            aes(trialNum, trialDuration, color = condition)) +
    scale_color_manual(values = c('blue', 'red', 'black')) + 
    xlab('Trial num') + ylab('Trial duration / s') + ggtitle(label) + displayTheme
  print(p)
  return(p)
}




# calculate kaplan-meier and area under the curve
kmsc <- function(blockData,tMax,blockLabel='',makePlot=FALSE,grid=0) {
  library(survival)
  waitDuration = blockData$timeWaited
  quitIdx = (blockData$trialEarnings == 0)
  # for rewarded trials, base the duration on the reward delivery time (not the subsequent response)
  waitDuration[!quitIdx] <- blockData$scheduledWait[!quitIdx]
  # fit the survival function
  kmfit <- survfit(Surv(waitDuration, quitIdx, type='right') ~ 1, 
                 type='kaplan-meier', conf.type='none', start.time=0, se.fit=FALSE)
  # extract elements of the survival curve object (?survfit.object)
  kmT = kmfit$time
  kmF = kmfit$surv
  # add a point at zero
  kmT = c(0, kmT)
  kmF = c(1, kmF)
  # keep only points up through tMax
  keepIdx = kmT<=tMax
  kmT <- kmT[keepIdx]
  kmF <- kmF[keepIdx]
  # extend the last value to exactly tMax
  kmT <- c(kmT, tMax)
  kmF <- c(kmF, tail(kmF,1))
  # calculate auc
  auc <- sum(diff(kmT) * head(kmF,-1))
  # plot if requested
  if (makePlot) {
   plotData = data.frame(kmT = kmT, kmF = kmF)
   p = ggplot(plotData, aes(kmT, kmF)) + geom_line() + xlab('Delay (s)') +
      ylab('Survival rate') + ylim(c(0,1)) + xlim(c(0,tMax)) +
        ggtitle(sprintf('KMSC: %s (AUC = %1.1f)',blockLabel,auc)) + 
        displayTheme
   print(p)
  }
  # put the survival curve on a standard grid
  kmOnGrid = vector()
  for (gIdx in 1:length(grid)) {
    g = grid[gIdx]
    # use the last point where t is less than or equal to the current grid value
    kmOnGrid[gIdx] = kmF[max(which(kmT<=g))]
  }
  return(list(kmT=kmT, kmF=kmF, auc=auc, kmOnGrid=kmOnGrid))
}

# this function can truncate trials in the simualtion object
# which enables us to zoom in and look and specific trials
truncateTrials = function(data, startTidx, endTidx){
  nVar = length(data)
  varNames = names(data)
  outputs = vector(mode = "list", length = nVar)
  for(i in 1 : nVar){
    junk = data[[varNames[i]]]
    outputs[[varNames[i]]] = junk[startTidx:endTidx]
  }
  return(outputs)
}

# correlation plot
# the first col of plotData is x, the second col is y, the third col is the group
plotCorrelation = function(data, dotColor, cor.method, isRank){
  conditions = c("HP", "LP")
  colnames(data) = c("x", "y", "cond")
  
  # calculate correlations
  corTests = lapply(1:2, function(i) cor.test(data[data$cond == conditions[i], "x"],
                                              data[data$cond == conditions[i], "y"],
                                              method = cor.method))
  
  rhos = sapply(1:2, function(i) round(corTests[[i]]$estimate, 3))
  ps = sapply(1:2, function(i) round(corTests[[i]]$p.value, 3))
  textColors = ifelse(ps < 0.05, "red", "blue")
  textData = data.frame(label = paste(rhos, "(p =", ps, ")"),
                        cond= c("HP", "LP"), color = textColors)
  # prepare rank 
  if(isRank){
    plotData = data %>% group_by(cond) %>% mutate(xRank = rank(x), yRank = rank(y))
  }

  
  # plot
  if(isRank){
    p0 = ggplot(plotData, aes(xRank, yRank)) + geom_point(size = 4, color = dotColor, fill = dotColor)
  }else{
    p0 = ggplot(plotData, aes(x, y)) + geom_point(size = 4, color = dotColor, fill = dotColor)
  }
  p = p0  + geom_text(data = textData,aes(x = -Inf,y = -Inf, label = label),
              hjust   = -0.2,vjust = -1,color = "blue",size = 5, fontface = 2, color = textColors) +
    facet_grid(~cond)
 return(p)
} 
