load("wtwSettings.RData")
blockSecs = 10 * 60
n = ceiling(blockSecs / iti)
scheduledWait.HP = vector(length = n)
seq = c();
for(i in 1 : n){
  tempt= drawSample("unif16", seq)
  scheduledWait.HP[i] = tempt[["delay"]]
  seq = tempt[["seq"]]
}

scheduledWait.LP = vector(length = n)
seq = c()
for(i in 1 : n){
  tempt= drawSample("logspace_1.75_32", seq)
  scheduledWait.LP[i] = tempt[["delay"]]
  seq = tempt[["seq"]]
}
save("scheduledWait.LP", "scheduledWait.HP", file = "scheduledWait.RData")
write.table(scheduledWait.HP,  file = "scheduledWaitHP.csv", sep = ",",
            col.names = F, row.names=FALSE)
write.table(scheduledWait.LP,  file = "scheduledWaitLP.csv", sep = ",",
            col.names = F, row.names=FALSE)