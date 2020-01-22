# experiment design paramters
conditions = c("HP", "LP")
     tMaxs = c(20, 40) # max trial durations in secs
    nBlock = 3
 blockMin = 7 # block duration in mins
 blockSec = blockMin * 60 # block duration in secs
       iti = 2 # iti duration in secs
tokenValue = 10 # value of the matured token
optimRewardRates = list( HP = 5 /6, LP = 0.930925) # unit : cent / sec
optimWaitThresholds = list(HP = 20, LP = 2.2) # unit : sec
# analyses parameters
tGrid = seq(0, blockSec * nBlock - 1, by = 1) # time grid for wtw time courses, open interval 
kmGrid = seq(0, min(tMaxs), by = 0.1) # time grid for Kaplan-Meier survival curves
save("conditions" = conditions,
     "tMaxs" = tMaxs,
     "blockMin" = blockMin,
     "blockSec" = blockSec,
     "nBlock" = nBlock,
     "iti" = iti,
     "tokenValue" = tokenValue,
     "optimRewardRates" = optimRewardRates,
     "optimWaitThresholds" = optimWaitThresholds,
     'tGrid' = tGrid,
     'kmGrid' = kmGrid,
     file = "expParas.RData")

