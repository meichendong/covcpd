# S04
# ORIGINAL CODE: E0009, E0016

library(pheatmap)
library(umap)
library(huge)
library(JGL)
library(PRROC)
library(ecp)
library(Rcpp)
library(RcppArmadillo)
# 
# setwd("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/data")
# source("C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/scripts/FUNCTIONS.R")
# Rcpp::sourceCpp('C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/scripts/FUNCTIONS_CPP_JGNSC.cpp')
# Rcpp::sourceCpp('C:/Users/meichen/OneDrive - University of North Carolina at Chapel Hill/994/medulloblastoma/scripts/HDtest_CPP.cpp')

setwd("/pine/scr/m/e/meichen/hovestdta/data/E0016")
source("/pine/scr/m/e/meichen/hovestdta/scripts/FUNCTIONS.R")
# sourceCpp("/pine/scr/m/e/meichen/hovestdta/scripts/ecp_cpp.cpp")
Rcpp::sourceCpp('/pine/scr/m/e/meichen/hovestdta/scripts/HDtest_CPP.cpp')
Rcpp::sourceCpp('/pine/scr/m/e/meichen/hovestdta/scripts/FUNCTIONS_CPP_JGNSC.cpp')



# -----------------------------------
# generate MVN  data
# -----------------------------------
# NCHANGE=1
# simnum = 1
# nsample = 100

set.seed(simnum)
k0=50
# if(nsample > 199){
#   k0 = 99
# } else if(nsample > 99){
#   k0=40
# }

# ===== NCP = 2, 3 CONDITIONS ======= #
if(NCHANGE ==1){
  if ("scenario" == "diffstructure"){
    # same mean, diff structure
    nivec.list.diff <- list(nivec= c(40,60,100),nivec2 = c(100,100))
    diffblk = list(1:2, 1)
    # diffblk = list( 2:3,2:5)
    sigma.list.1 <- generateSigmaList(nivec.list = nivec.list.diff, structure = "Diff S, Identical W", diffblk = diffblk)
    count1 <- CountMap5(sigma = sigma.list.1[[1]], ngene = 200, n = nsample)
    count2 <- CountMap6(sigma = sigma.list.1[[2]], ngene = 200, n = nsample, scmeanvec = count1$thetaj, pijvec = count1$pij)
  } else if ("scenario" == "samestructure"){
    # same structure, different mean
    nivec.list <- list(nivec= c(100,100),nivec2 = c(100,100))
    sigma.list.1 <- generateSigmaList(nivec.list = nivec.list, structure = "Identical S, Diff W", diffblk = NULL)
    count1 <- CountMap5(sigma = sigma.list.1[[1]], ngene = 200, n = nsample)
    count2 <- CountMap5(sigma = sigma.list.1[[2]], ngene = 200, n = nsample)
  } else if ("scenario" == "diffstrdiffmean"){
    # diff mean, diff structure
    nivec.list.diff <- list(nivec= c(40,60,100),nivec2 = c(100,100))
    diffblk = list( 1:2,1)
    sigma.list.1 <- generateSigmaList(nivec.list = nivec.list.diff, structure = "Diff S, Identical W", diffblk = diffblk)
    count1 <- CountMap5(sigma = sigma.list.1[[1]], ngene = 200, n = nsample)
    count2 <- CountMap5(sigma = sigma.list.1[[2]], ngene = 200, n = nsample)
  } 
  
  # if ("cpxxx" == "null"){
  y = count1$count 
  # } else if("cpxxx" == "cp1"){
  #   y = rbind(count1$count,
  #             count2$count)
  # }
}


if("yform" =="raw"){
  yreplace = y
} else if("yform" =="log"){
  yreplace = logy
} else if ("yform" =="imputeNPN"){
  yreplace = yimpnpn
}


# -----------------------------

# ------------- mean correction by ecp
library(ecp)

mcp = c(1, sample(30:(nsample-30),1), nsample+1)
colnorm_y = NULL
for(i in 1:2){
  sgroupi = yreplace[mcp[i]:(mcp[i+1]-1),]
  colnorm_sgroupi = apply(sgroupi,2, function(x){
    (x-mean(x)) #/sd(x)
  })
  colnorm_y = rbind(colnorm_y, colnorm_sgroupi)
}


# ===============================
# my cai testing E0002 cpd_cov, test_cpd.R
# myres = cpd_cov(X = t(y),k=3, pthresh = 0.1)
res1 = NULL; tp=NULL; KHAT =NULL; prec=0; rec=0; F1score=0; AnnoErrPct = 1;
t1 = proc.time()
# my2 = cpd_cov(X = t(yreplace), k=k0, search_quantile = 0.5, maxseg = NCHANGE + 2, pp = 1/powercp)
res = covcpd(X = t(colnorm_y), k=k0, maxseg = 3, search_by = NULL, nperm = 1000, siglvl = 0.05)
t2 = proc.time()





res1 = c(paste(res$cps, collapse = ";"), 
         paste(round(res$pvals, 8), collapse = ";"),
         paste(unlist(res$candtau),";"), 
         quantile(res$minpdist, seq(5, 95, by=5)/100),
         simnum, nsample, 
         mcp[2])

write.table(res1, "E0016_simnum_nsample.txt")
