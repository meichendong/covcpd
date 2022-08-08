#' covcpd: detect covariance-based change points for scRNA-seq data
#'
#' @param X raw count matrix, genes by cells
#' @param nperm number of permutations
#' @param siglvl significance level, default is 0.05
#' @param k the number of minimum cells in each segment, default is 30.
#' @param maxseg the maximum number of segments
#' @param verbose logical, whether to print out processing information. Default is TRUE.
#' @param search_by the number of cells between two consecutive candidate change points.
#' @param learnMode to learn the curve for Gumbel distribution parameters alpha and beta, we provide either "linear" or "GAM" method 
#' @return cps: the estimated change points location.
#'         cps_order: the estimated change points location ordered by the order they were detected.
#'         pvals: the p-values for the detected cps
#'         candtau: candidate change point locations that were tested
#'         pthresh: p-value threshold learned by permutation
#'         minpdist: minimum p-value distribution learned by permuttaion.
#' @export

covcpd <- function(X, nperm = 500, siglvl = 0.05, k=30, maxseg = NULL, 
                   verbose = T, search_by = NULL, learnMode = c("linear","GAM")[1]){
  n <- ncol(X)
  cp.1o = cp.1 <- c(1,n) # NULL SET FOR START
  pvec = NULL
  candtau = NULL
  
  # find a cp from the whole interval first
  seg <- 1
  m = length(cp.1)
  if(n/k < 3){
    k = round(n/3)
  }
  if(is.null(maxseg)){
    maxseg = max(3, ceiling(n/k))
  }
  if (is.null(search_by)){
    if(n>1000){
      search_by = ceiling(n/100)
    } else if(n>500){
      search_by = 5
    } else if(n <=500){
      search_by = 2
    }
  }
  
  #~~~~~~~~~~~~ 1. learn mu and beta
  if (verbose){
    message("Start learning Gumbel distr. parameters mu and beta... \n")
  }
  by1 = 0.3 
  
  nr =seq(1,(n-k)/k, by1)
  searchx = unique(c(round(n/(1+nr)), n - round(n/(1+nr))))
  searchx = searchx[order(searchx)]
  
  learngum.perm = CLXPermHier(t(X), searchx = searchx, start = 1,nt2 = 4, minband = k, nperm = nperm)
  learngum.perm = as.data.frame(learngum.perm) 
  learngum.perm$nratio = pmax(learngum.perm$V2, learngum.perm$V4-learngum.perm$V2)/pmin(learngum.perm$V2, learngum.perm$V4 - learngum.perm$V2)
  
  dtnullmed1 = aggregate(learngum.perm$V1, list(learngum.perm$nratio, learngum.perm$V4), median)
  dtnullvar1 = aggregate(learngum.perm$V1, list(learngum.perm$nratio, learngum.perm$V4), var)
  dtnullbeta1 = sqrt(6*dtnullvar1$x/pi^2)
  dtnullmu1  = dtnullmed1$x + dtnullbeta1*log(log(2))
  
  dtx1 = data.frame(mu = dtnullmu1, beta=dtnullbeta1, nratio = dtnullmed1$Group.1, n1n2 = dtnullmed1$Group.2)
  dtx1$nratio2 = dtx1$nratio^2
  dtx1$n1n2sq = (dtx1$n1n2)^2
  
  if (learnMode == "linear"){
    lm1 = lm(mu ~ nratio + nratio2 + n1n2 + n1n2sq, data = dtx1)
    lm2 = lm(beta ~ nratio + nratio2 + n1n2 + n1n2sq, data = dtx1)
  } else if (learnMode == "GAM"){ 
    gam1 = mgcv::gam(mu ~ s(nratio) + s(n1n2) + nratio*n1n2, data = dtx1) 
    gam2 = mgcv::gam(beta ~ s(nratio) + s(n1n2) + nratio*n1n2, data = dtx1) 
  }
  
  
  # ~~~~~~~~~~~~~~~ 2. learn p-threshold
  if (verbose){
    message("Start learning type-I error threshold by permutation... \n")
  } 
  searchx = seq(k,(n-1 - k), by = search_by)
  learnp.perm = CLXPermHier(t(X), searchx = searchx, start = 1,nt2 = 4, minband = k, nperm = nperm)
  learnp.perm = as.data.frame(learnp.perm)
  learnp.perm$nratio = pmax(learnp.perm$V2, learnp.perm$V4-learnp.perm$V2)/pmin(learnp.perm$V2, learnp.perm$V4 - learnp.perm$V2)
  newdt = data.frame(nratio = learnp.perm$nratio, n1n2 = learnp.perm$V4)
  
  if (learnMode == "linear"){
    newdt$nratio2 = newdt$nratio^2
    newdt$n1n2sq = (newdt$n1n2)^2
    est_mu = predict(lm1, newdata = newdt)
    est_beta = predict(lm2, newdata = newdt)
  } else if (learnMode == "GAM"){
    est_mu = mgcv::predict.gam(gam1, newdata = newdt)
    est_beta = mgcv::predict.gam(gam2, newdata = newdt)
  }
  
  learnp.perm$pgum = 1- evd::pgumbel(learnp.perm$V1, loc = est_mu, scale = est_beta)
  learnp_acrosst2 = aggregate(learnp.perm$pgum, list(t1 = learnp.perm$V2, permindex = learnp.perm$V3), min)
  learnp_acrossperm = aggregate(learnp_acrosst2$x, list(learnp_acrosst2$permindex), min)
  pthresh = quantile(learnp_acrossperm$x, siglvl)   
  
  # start search
  # ~~~~
  if (verbose){
    message("Start searching for change points... \n")
  }
  while (seg < maxseg){ 
    newcps = NULL
    for (jj in 1:(m-1)){ 
      if (cp.1[jj] > 1){
        start.p = cp.1[jj] +1
      } else {
        start.p = 1
      }
      end.p = cp.1[jj+1] 
      if ((end.p - start.p) <= 2*k){
        message("skip interval (", start.p, ", ", end.p, ")\n")
      } else { # keep finding cp
        
        stat.all <- NULL 
        searchx = seq((start.p + k),(end.p-1 - k), by = search_by) 
        
        for (t1 in searchx){ 
          t2seq = c(seq(t1+1+k, end.p, by = max(round((end.p - t1-1-k)/4), 1)), end.p)
          pt2   = 1
          pgums = NULL
          mns = NULL
          for(tt in 1:length(t2seq)){
            X1 <- as.matrix(t(X[,start.p:t1]))
            X2 <- as.matrix(t(X[,(t1+1):t2seq[tt]]))  
            X1 = getMatNPN(X1, nrow(X1))
            X2 = getMatNPN(X2, nrow(X2))
            stat0 = testCov_cpp(X=X1,Y=X2)
            nr = max(nrow(X1),nrow(X2))/min(nrow(X1),nrow(X2))
            
            tmpdt = data.frame(nratio = nr, n1n2 = nrow(X1)+ nrow(X2))
            
            if (learnMode == "linear"){
              tmpdt$nratio2 = tmpdt$nratio^2
              tmpdt$n1n2sq = (tmpdt$n1n2)^2
              est_mu_tmp = predict(lm1, newdata = tmpdt)
              est_beta_tmp = predict(lm2, newdata = tmpdt)
            } else if (learnMode == "GAM"){
              est_mu_tmp = mgcv::predict.gam(gam1, tmpdt)
              est_beta_tmp = mgcv::predict.gam(gam1, tmpdt)
            }
            pgum = 1- evd::pgumbel(stat0$CLX, loc = est_mu_tmp, scale = est_beta_tmp)
            pgums = rbind(pgums, c(pgum, stat0$pvalue, t2seq[tt]))
            mns = rbind(mns, c( stat0$CLX, est_mu_tmp))
            pt2 = min(pgum, pt2)
          }     
          # permutation step for each candidate point
          stat.all <- rbind(stat.all, c(max(abs(mns[,1] - mns[,2])),pt2, t1))
        }
        stat.all <- as.data.frame(stat.all) 
        
        if(length(searchx)>5){ 
          minp = min(stat.all$V2)
          taus = stat.all$V3[stat.all$V2 == minp] 
          if (length(taus)>1){ 
            stmp = stat.all[stat.all$V3 %in% taus,]
            tau0 <- stmp$V3[which.max(stmp$V1)]
          } else {
            tau0 = taus
          }
          
          if (search_by > 2){
            candidates = seq(tau0 - search_by, tau0 + search_by, by=1)
            stat.all.cand = NULL
            for (t1 in candidates){ 
              t2seq = c(seq(min(t1+1+k, end.p), end.p, by = max(round((end.p - t1-1-k)/5), 2)), end.p)
              pt2   = 1
              pgums = NULL
              mns = NULL
              for(tt in 1:length(t2seq)){
                X1 <- as.matrix(t(X[,start.p:t1]))
                X2 <- as.matrix(t(X[,(t1+1):t2seq[tt]]))   
                X1 = getMatNPN(X1, nrow(X1))
                X2 = getMatNPN(X2, nrow(X2))
                stat0 = testCov_cpp(X=X1,Y=X2)
                nr = max(nrow(X1),nrow(X2))/min(nrow(X1),nrow(X2))
                
                tmpdt = data.frame(nratio = nr, n1n2 = nrow(X1)+ nrow(X2))
                
                if (learnMode == "linear"){
                  tmpdt$nratio2 = tmpdt$nratio^2
                  tmpdt$n1n2sq = (tmpdt$n1n2)^2
                  est_mu_tmp = predict(lm1, newdata = tmpdt)
                  est_beta_tmp = predict(lm2, newdata = tmpdt)
                } else if (learnMode == "GAM"){
                  est_mu_tmp = mgcv::predict.gam(gam1, tmpdt)
                  est_beta_tmp = mgcv::predict.gam(gam1, tmpdt)
                }
                pgum = 1- evd::pgumbel(stat0$CLX, loc = est_mu_tmp, scale = est_beta_tmp)
                pgums = rbind(pgums, c(pgum, stat0$pvalue, t2seq[tt]))
                mns = rbind(mns, c( stat0$CLX, est_mu_tmp))
                pt2 = min(pgum, pt2)
              }     
              stat.all.cand <- rbind(stat.all.cand, c(max(abs(mns[,1] - mns[,2])),pt2, t1))
            }
            
            stat.all.cand = as.data.frame(stat.all.cand)
            minp = min(stat.all.cand$V2)
            taus = stat.all.cand$V3[stat.all.cand$V2 == minp] 
            
            if (length(taus)>1){ 
              stmp = stat.all.cand[stat.all.cand$V3 %in% taus,]
              tau0 <- stmp$V3[which.max(stmp$V1)]
            } else {
              tau0 = taus
            }
          }  
          
        } else { 
          minp = min(stat.all$V2)
          tau0 = stat.all$V3[which.min(stat.all$V2)]
        } 
        
        pvec = c(pvec, minp)
        candtau = c(candtau, tau0)
        if (minp < pthresh){
          newcps = c(newcps, tau0)  
        }  
        
      }
      
    } 
    
    cp.1o <- append(cp.1o, newcps)
    cp.1 <- unique(cp.1o)
    cp.1 <- cp.1[order(cp.1)]
    
    if (length(cp.1) > m){
      m = length(cp.1)
      seg = m - 1
    } else {
      break
    }
  }
  return(list(cps = cp.1,
              cps_order = cp.1o,
              pvals = pvec,
              candtau = candtau,
              pthresh = pthresh,
              minpdist = learnp_acrossperm$x))
  
}