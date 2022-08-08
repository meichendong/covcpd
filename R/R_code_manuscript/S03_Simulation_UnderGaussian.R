# S03

# simulation under Gaussian assumption
library(pheatmap)
library(umap)
library(huge)
library(JGL)
library(PRROC)
library(ecp)
library(Rcpp)
library(RcppArmadillo)

# setwd("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data")
# source("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/scripts/FUNCTIONS.R")
# Rcpp::sourceCpp('C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/scripts/FUNCTIONS_CPP_JGNSC.cpp')

setwd("/pine/scr/m/e/meichen/hovestdta/data/E0019")
source("/pine/scr/m/e/meichen/hovestdta/scripts/FUNCTIONS.R")
# sourceCpp("/pine/scr/m/e/meichen/hovestdta/scripts/ecp_cpp.cpp")
# Rcpp::sourceCpp('C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/scripts/HDtest_CPP.cpp')
Rcpp::sourceCpp('/pine/scr/m/e/meichen/hovestdta/scripts/HDtest_CPP.cpp')
Rcpp::sourceCpp('/pine/scr/m/e/meichen/hovestdta/scripts/FUNCTIONS_CPP_JGNSC.cpp')


# -----------------------------------
# generate MVN  data
# -----------------------------------
set.seed(simnum)
k0=50
if(nsample > 500){
  k0 = min(100, nsample/10)
}

# ===== NCP = 2, 3 CONDITIONS ======= #
if(NCHANGE ==1){
  if ("scenario" == "diffstructure"){
    # same mean, diff structure
    nivec.list.diff <- list(nivec= c(40,60,100),nivec2 = c(100,100))
    diffblk = list( 1:2,1)
    sigma.list.1 <- generateSigmaList(nivec.list = nivec.list.diff, structure = "Diff S, Identical W", diffblk = diffblk)
    # count1 <- CountMap5(sigma = sigma.list.1[[1]], ngene = 200, n = nsample)
    # count2 <- CountMap6(sigma = sigma.list.1[[2]], ngene = 200, n = nsample, scmeanvec = count1$thetaj, pijvec = count1$pij)
    
    count1 = generateMVN(sigma = sigma.list.1[[1]], n = nsample, meanvec = NULL)
    count2 = generateMVN(sigma = sigma.list.1[[2]], n = nsample, meanvec = NULL)
    
  } else if ("scenario" == "samestructure"){
    # same structure, different mean
    nivec.list <- list(nivec= c(100,100),nivec2 = c(100,100))
    sigma.list.1 <- generateSigmaList(nivec.list = nivec.list, structure = "Identical S, Diff W", diffblk = NULL)
    # count1 <- CountMap5(sigma = sigma.list.1[[1]], ngene = 200, n = nsample)
    # # count2 <- CountMap5(sigma = sigma.list.1[[1]], ngene = 200, n = nsample)
    # count2 <- CountMap6(sigma = sigma.list.1[[1]], ngene = 200, n = nsample, scmeanvec = count1$thetaj, pijvec = count1$pij, meanDiffPct = 0.2)
    mMeanvec = rep(0,200)
    mMeanvec2 = mMeanvec
    mMeanvec2[sample(1:200,0.2*200)] = runif(0.2*200,-2,2)
    count1 = generateMVN(sigma = sigma.list.1[[1]], n = nsample, meanvec = mMeanvec)
    count2 = generateMVN(sigma = sigma.list.1[[1]], n = nsample, meanvec = mMeanvec2)
    
  } else if ("scenario" == "diffstrdiffmean"){
    # diff mean, diff structure
    nivec.list.diff <- list(nivec= c(40,60,100),nivec2 = c(100,100))
    diffblk = list( 1:2,1)
    sigma.list.1 <- generateSigmaList(nivec.list = nivec.list.diff, structure = "Diff S, Identical W", diffblk = diffblk)
    # count1 <- CountMap5(sigma = sigma.list.1[[1]], ngene = 200, n = nsample)
    # # count2 <- CountMap5(sigma = sigma.list.1[[2]], ngene = 200, n = nsample)
    # # set % of means are diff: not all means are diff. they should be independent with structure? 
    # # set 20% ~ 30% means diff, the rest are the same
    # count2 <- CountMap6(sigma = sigma.list.1[[2]], ngene = 200, n = nsample, scmeanvec = count1$thetaj, pijvec = count1$pij, meanDiffPct = xdiffpct/100)
    mMeanvec = rep(0,200)
    mMeanvec2 = mMeanvec
    mMeanvec2[sample(1:200,0.2*200)] = runif(0.2*200,-2,2)
    count1 = generateMVN(sigma = sigma.list.1[[1]], n = nsample, meanvec = mMeanvec)
    count2 = generateMVN(sigma = sigma.list.1[[2]], n = nsample, meanvec = mMeanvec2)
  } 
  y = rbind(count1,
            count2) 
}



if("yform" =="raw"){
  yreplace = y
} else if("yform" =="log"){
  yreplace = logy
} else if ("yform" =="imputeNPN"){
  yreplace = yimpnpn
}


sumtable = NULL
# three metrics: 1. computing time; 2. annotation error; 3. F1 score
M_thresh = 0.1*nsample
T0 = seq(nsample, nsample*NCHANGE, nsample)
Kstar = NCHANGE
Tstar = rep(1:(1+NCHANGE), each = nsample)



# ------------- mean correction by ecp
library(ecp)
res.ecp = e.divisive(yreplace, alpha=1)

# res1 = res.ecp$estimates[c(-1,-length(res.ecp$estimates))]
nsubgroup = length(res.ecp$estimates) - 1
mcp = res.ecp$estimates
colnorm_y = NULL
for(i in 1:nsubgroup){
  sgroupi = yreplace[mcp[i]:(mcp[i+1]-1),]
  colnorm_sgroupi = apply(sgroupi,2, function(x){
    (x-mean(x)) #/sd(x)
  })
  colnorm_y = rbind(colnorm_y, colnorm_sgroupi)
}

# plot(sgroupi[,1], colnorm_sgroupi[,1])
check.ecp = e.divisive(colnorm_y, alpha = 1)
check.ecp$estimates


# -----------------------------------
# test performance of packages
# -----------------------------------

# ===============================
# my cai testing E0002 cpd_cov, test_cpd.R
# myres = cpd_cov(X = t(y),k=3, pthresh = 0.1)
res1 = NULL; tp=NULL; KHAT =NULL; prec=0; rec=0; F1score=0; AnnoErrPct = 1;
t1 = proc.time()
# my2 = cpd_cov(X = t(yreplace), k=k0, search_quantile = 0.5, maxseg = NCHANGE + 2, pp = 1/powercp)
my2 = covcpd(X = t(colnorm_y), k=k0, maxseg = 6, search_by = NULL, nperm = 1000, siglvl = 0.05)
t2 = proc.time()

res1 = my2$cps[-c(1,length(my2$cps))]
KHAT = length(res1)
tp=0
SE =0
for(i in 1:Kstar){
  tp0 = any(abs(res1 - T0[i]) <= M_thresh)
  tp = tp0+tp
  if(tp0){
    tmp = (min(abs(res1 - T0[i])))^2
  } else {
    tmp = nsample^2
  } 
  SE = SE + tmp
} 
if(tp){
  prec = tp/KHAT
  rec  = tp/Kstar
  F1score = 2*(prec*rec)/(prec + rec)
} else {
  F1score = 0
} 
MSE = sqrt(SE)

F1score = ifelse(is.na(2*(prec*rec)/(prec + rec)),0,2*(prec*rec)/(prec + rec))
AnnoErrPct = abs(KHAT - Kstar)/Kstar

n_est = diff(c(1, res1, (NCHANGE+1)*nsample))
n_est[1] = n_est[1] +1 
That = NULL
for(i in 1:(KHAT+1)){
  tmp = rep(i, n_est[i])
  That = c(That, tmp)
} 
ARI1 = ARI_mclust(Tstar, That)


test1 = data.frame(package = "covcpd",
                   processtime = (t2-t1)[3],
                   cp = paste(res1, collapse = "; "),
                   prec = prec,
                   recall = rec,
                   F1score = F1score,
                   AnnoErrPct = AnnoErrPct,
                   ARI = ARI1,
                   MSE = MSE)
sumtable = rbind.data.frame(sumtable, test1)


# -----------------------------------
# test performance of packages
# -----------------------------------
library(ecp)
res1 = NULL; tp=NULL; KHAT =NULL; prec=0; rec=0; F1score=0; AnnoErrPct = 1; ARI =0;
t1 = proc.time()
res.ecp = e.divisive(yreplace, alpha=1)
t2 = proc.time() 

res1 = res.ecp$estimates[c(-1,-length(res.ecp$estimates))]
tps=0
SE = 0
for(i in 1:Kstar){
  tp = any(abs(res1 - T0[i]) < M_thresh)
  if(tp){
    tmp = (min(abs(res1 - T0[i])))^2
  } else {
    tmp = nsample^2
  }
  tps = tp+tps
  SE = SE + tmp
}
MSE = sqrt(SE)
KHAT = length(res1)
prec = tps/KHAT
rec  = tps/Kstar
F1score = 2*(prec*rec)/(prec + rec)
AnnoErrPct = abs(KHAT - Kstar)/Kstar
That = res.ecp$cluster
ARI1 = ARI_mclust(Tstar, That)


test1 = data.frame(package = "ecp",
                   processtime = (t2-t1)["elapsed"],
                   cp = paste(res1, collapse = "; "),
                   prec = prec,
                   recall = rec,
                   F1score = F1score,
                   AnnoErrPct = AnnoErrPct,
                   ARI = ARI1,
                   MSE = MSE)
sumtable = rbind.data.frame(sumtable, test1)


# ===============================
library(ccid)
# t1 = proc.time()
# res.ccid = detect.th(y, approach = 'infinity', scales = -1)
# t2 = proc.time() 
res1 = NULL; tp=NULL; KHAT =NULL; prec=0; rec=0; F1score=0; AnnoErrPct = 1;
t1 = proc.time()
res.ccid = detect.ic(yreplace, approach = 'euclidean', scales = -1)
t2 = proc.time() 

res1 = res.ccid$changepoints
if(is.na(res1)){
  res1 = NULL
}
if(length(res1) >1){
  res1 = res1[seq(2,length(res1), by=2)]
}
tp=0
SE =0
for(i in 1:Kstar){
  tp0 = any(abs(res1 - T0[i]) < M_thresh)
  tp = tp0+tp
  if(tp0){
    tmp = (min(abs(res1 - T0[i])))^2
  } else {
    tmp = nsample^2
  } 
  SE = SE + tmp
} 
MSE = sqrt(SE)
KHAT = length(res1)
prec = ifelse(is.na(tp/KHAT),0,tp/KHAT)
rec  = tp/Kstar
F1score = ifelse(is.na(2*(prec*rec)/(prec + rec)),0,2*(prec*rec)/(prec + rec))
AnnoErrPct = abs(KHAT - Kstar)/Kstar

n_est = diff(c(1, res1, (NCHANGE+1)*nsample))
n_est[1] = n_est[1] +1 
That = NULL
for(i in 1:(KHAT+1)){
  tmp = rep(i, n_est[i])
  That = c(That, tmp)
} 
ARI1 = ARI_mclust(Tstar, That)

test1 = data.frame(package = "ccid -- euclidean",
                   processtime = (t2-t1)[3],
                   cp = paste(res1, collapse = "; "),
                   prec = prec,
                   recall = rec,
                   F1score = F1score,
                   AnnoErrPct = AnnoErrPct,
                   ARI = ARI1,
                   MSE = MSE)
sumtable = rbind.data.frame(sumtable, test1)


# ===============================
library(changepoint.geo)
res1 = NULL; tp=NULL; KHAT =NULL; prec=0; rec=0; F1score=0; AnnoErrPct = 1;
##Variance change in all series
t1 = proc.time()
res.geo <- geomcp(yreplace)
t2 = proc.time() 
z1 = res.geo@dist.cpts #mean change
z2 = res.geo@ang.cpts #variance change

res1 = z1
tp=0
SE =0
for(i in 1:Kstar){
  tp0 = any(abs(res1 - T0[i]) < M_thresh)
  tp = tp0+tp
  if(tp0){
    tmp = (min(abs(res1 - T0[i])))^2
  } else {
    tmp = nsample^2
  } 
  SE = SE + tmp
} 
MSE = sqrt(SE)
KHAT = length(res1)
prec = ifelse(is.na(tp/KHAT),0,tp/KHAT)
rec  = tp/Kstar
F1score = ifelse(is.na(2*(prec*rec)/(prec + rec)),0,2*(prec*rec)/(prec + rec))
AnnoErrPct = abs(KHAT - Kstar)/Kstar

n_est = diff(c(1, res1, (NCHANGE+1)*nsample))
n_est[1] = n_est[1] +1 
That = NULL
for(i in 1:(KHAT+1)){
  tmp = rep(i, n_est[i])
  That = c(That, tmp)
} 
ARI1 = ARI_mclust(Tstar, That)
test1 = data.frame(package = "changepoint.geo: mean change",
                   processtime = (t2-t1)[3],
                   cp = paste(z1, collapse = "; "),
                   prec = prec,
                   recall = rec,
                   F1score = F1score,
                   AnnoErrPct = AnnoErrPct,
                   ARI = ARI1,
                   MSE = MSE)
sumtable = rbind.data.frame(sumtable, test1)

res1 = NULL; tp=NULL; KHAT =NULL; prec=0; rec=0; F1score=0; AnnoErrPct = 1;
res1 = z2
tp=0
SE =0
for(i in 1:Kstar){
  tp0 = any(abs(res1 - T0[i]) < M_thresh)
  tp = tp0+tp
  if(tp0){
    tmp = (min(abs(res1 - T0[i])))^2
  } else {
    tmp = nsample^2
  } 
  SE = SE + tmp
} 
MSE = sqrt(SE)
KHAT = length(res1)
prec = ifelse(is.na(tp/KHAT),0,tp/KHAT)
rec  = tp/Kstar
F1score = ifelse(is.na(2*(prec*rec)/(prec + rec)),0,2*(prec*rec)/(prec + rec))
AnnoErrPct = abs(KHAT - Kstar)/Kstar

n_est = diff(c(1, res1, (NCHANGE+1)*nsample))
n_est[1] = n_est[1] +1 
That = NULL
for(i in 1:(KHAT+1)){
  tmp = rep(i, n_est[i])
  That = c(That, tmp)
} 
ARI1 = ARI_mclust(Tstar, That)
test1 = data.frame(package = "changepoint.geo: cov change",
                   processtime = (t2-t1)[3],
                   cp = paste(z2, collapse = "; "),
                   prec = prec,
                   recall = rec,
                   F1score = F1score,
                   AnnoErrPct = AnnoErrPct,
                   ARI = ARI1,
                   MSE = MSE)
sumtable = rbind.data.frame(sumtable, test1)

# ---------------

sumtable$truecp = paste(T0, collapse = ";")
sumtable$simu  = simnum
sumtable$ncp   = NCHANGE
sumtable$type  = "scenario"
sumtable$normalize = "yform"
# sumtable$diffpct = xdiffpct


if(simnum == 1 & NCHANGE == 1 & "scenario" == "samestructure" & nsample == 200){
  save.image("E0019_test.rda")
}

write.table(sumtable, "E0019_CP_NCHANGE_SIM_simnum_NSAMnsample_scenario_yform_covcpd.txt")

