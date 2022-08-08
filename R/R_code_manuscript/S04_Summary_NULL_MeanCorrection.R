# S04 summary

#mE0016SUMMARY E0009V2 SUMMARY
# summary of E0009 longleaf results
setwd("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/E0016")
lfn = list.files()

tnull = NULL
for(i in 1:length(lfn)){
  tmp = read.table(lfn[i])
  tmp2 = t(tmp[c(1,2,3,nrow(tmp)-1, nrow(tmp)),])
  tnull = rbind(tnull, tmp2)
}

nsample = c(200,300,500)

t1etable = NULL
for(nn in nsample){
  tmp1 = as.data.frame(tnull[tnull[,4] == nn,])
  tmp1$detect = tmp1$V1 != paste("1;",nn, sep = "")
  t1e = table(tmp1$detect)[2] / table(tmp1$detect)[1]
  t1etable = rbind(t1etable, c(nn,t1e))
}
t1etable

# RESULT FROM 500 SIMULATIONS:
#             TRUE
# [1,] 200 0.026694045
# [2,] 300 0.014198783
# [3,] 500 0.004016064