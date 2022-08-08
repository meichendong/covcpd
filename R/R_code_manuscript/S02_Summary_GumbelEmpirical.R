# S02 -- SUMMARY

# RESULT FROM ONE SIMULATION. HOW TO INTEGRATE MULTIPLE?

# empirical distribution of mu and beta
library(evd)

nsample = 400
load(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data/E0012/E0012_GUMBEL_null_1_SAM",nsample,".rda", sep = ""))
stat.all$nratio = round(pmax(stat.all$V3, nsample-stat.all$V3)/pmin(stat.all$V3, nsample-stat.all$V3), 2)
# jvec = c(45,29,20,15,12,10,8,5)
jvec = c(96,62,45,35,29,24,20,17)
# jvec = c(146,95,71,55,45,38,33,28)
# jvec = c(196,129,95,75,62,52,45,40)
# jvec = c(246,162,120,95,79,67,58,51)

pdf(paste("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/covcpd/manuscript/submit/data/S02_N",nsample,"_GumbelQQplot.pdf", sep = ""), height = 10, width = 6)
par(mfrow = c(4,2))
for(j in jvec){
  tmp = stat.allperm[stat.allperm[,2] == stat.all[j,3],]
  med1 = median(tmp[,1])
  var1 = var(tmp[,1])
  beta = sqrt(6*var1/pi^2)
  mu   = med1 + beta*log(log(2))
  
  xgumbel = rgumbel(1000, loc = mu, scale = beta)
  
  qqplot(tmp[,1], xgumbel, xlab = "Mn calculated from permutation", ylab = paste("Gumbel(mu=",round(mu,1),",beta=",round(beta,1),")"), main = paste("Q-Q plot \n Nsample=",nsample,", nratio=",round(stat.all[j,4],1)))
  abline(a=0,b=1, col="red")
}
dev.off()
