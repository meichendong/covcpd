# S01 -- SIMULATION UNDER GAUSSIAN DISTRIBUTION

setwd("/pine/scr/m/e/meichen/hovestdta/data/E0012")
source("/pine/scr/m/e/meichen/hovestdta/scripts/FUNCTIONS.R") 
Rcpp::sourceCpp('/pine/scr/m/e/meichen/hovestdta/scripts/HDtest_CPP.cpp')
library(evd)


set.seed(simnum)
N = nsample
nperm = 1000

nivec.list.diff <- list(nivec= c(20,30,50),nivec2 = c(50,50))
diffblk = list( 1:2,1)
sigma.list.1 <- generateSigmaList(nivec.list = nivec.list.diff, structure = "Diff S, Identical W", diffblk = diffblk)
y1 = generateMVN(sigma = sigma.list.1[[1]], n=N)
X = t(y1)
n <- ncol(X)
cp.1 <- c(1,n) # NULL SET FOR START
m = length(cp.1) 
search_by = 2

start.p = 1
end.p = cp.1[2]
stat.all <- NULL
k = 10
searchx = seq((start.p + k),(end.p - k), by = search_by)
for (t1 in searchx){ 
  X1 <- t(X[,start.p:t1])
  X2 <- t(X[,(t1+1):end.p]) 
  stat0 = testCov_cpp(X=X1,Y=X2)
  # permutation step for each candidate point
  stat.all <- rbind(stat.all, c(stat0$CLX, stat0$pvalue, t1))
}
stat.all <- as.data.frame(stat.all)

Xsub = X[,start.p: end.p] 
stat.allperm = CLXstatPerm(t(Xsub), searchx = searchx, start = start.p, nperm = nperm)
reslist = list()

for (j in 1:nrow(stat.all)){
  tmp = stat.allperm[stat.allperm[,2] == stat.all[j,3],]
  med1 = median(tmp[,1])
  var1 = var(tmp[,1])
  beta = sqrt(6*var1/pi^2)
  mu   = med1 + beta*log(log(2))
  
  xgumbel = rgumbel(1000, loc = mu, scale = beta)
  reslist[[j]] = list(permstat = tmp[,1],
                      xgumbel = xgumbel,
                      beta = beta,
                      mu = mu) 
}
 

save.image("E0012_GUMBEL_null_simnum_SAMnsample.rda")
