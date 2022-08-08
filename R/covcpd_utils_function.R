# covcpd supporting functions

#' Generate a list of covariance matrices.
#' @return a list of covariance matrices.
#' @export

generateSigmaList <- function(nivec.list, u=c(-60:60)/100 , ud = c(-100:-60, 60:100)/100,
                              structure = "Identical S, Identical W", diffblk = NULL){
  # structure S, Weight W: "Identical S, Identical W", "Identical S, Diff W", "Diff S, Identical W", "Diff S, Diff W"
  # if using the structure "Diff S, Identical W", setting some structure different while keep the rest of structure and weights the same, specify diffblk = list(diffblk1=1, diffblk2=c(1,2),...)

  sigma.list <- list()
  # check if identical structure
  checkI <- check_ident_list(nivec.list)
  if (structure =="Identical S, Identical W"){
    if (checkI){ #------------------------------------------- if structures are identical
      blklist <- list()
      sigma <- matrix(0, nrow = 0, ncol = sum(nivec.list[[1]]))
      nblk <- length(nivec.list[[1]])
      for (b in 1:nblk){
        ni <- nivec.list[[1]][b]
        blklist[[b]] <- generateBlki(ni=ni, ud=ud)
        zeroleft <- matrix(0, nrow = ni, ncol = sum(nivec.list[[1]][0:(b-1)]))
        zeroright <- matrix(0, nrow = ni, ncol = ifelse(b<nblk,sum(nivec.list[[1]][(b+1):nblk]),0))
        temp <- blklist[[b]]$sigmam
        sigma <- rbind(sigma, cbind(zeroleft, temp, zeroright))
      }
      gnames = paste("gene",1:sum(nivec.list[[1]]), sep = "")
      rownames(sigma) <- gnames
      colnames(sigma) <- gnames
      for(ss in 1:length(nivec.list)){
        sigma.list[[ss]] <- sigma
      }
    } else {
      message("nivec.list and the selected network structure do NOT match ...\n")
      break()
    }

  } else if (structure =="Identical S, Diff W"){
    if (checkI){ #------------------------------------- if structures are identical
      for(ss in 1:length(nivec.list)){
        blklist <- list()
        sigma <- matrix(0, nrow = 0, ncol = sum(nivec.list[[ss]]))
        # sigma <- matrix(0, nrow = ni*nblk, ncol = ni*nblk)
        nblk <- length(nivec.list[[ss]])
        for (b in 1:nblk){
          ni <- nivec.list[[ss]][b]
          blklist[[b]] <- generateBlki(ni=ni, ud=ud)
          zeroleft <- matrix(0, nrow = ni, ncol = sum(nivec.list[[ss]][0:(b-1)]))
          zeroright <- matrix(0, nrow = ni, ncol = ifelse(b<nblk,sum(nivec.list[[ss]][(b+1):nblk]),0))
          temp <- blklist[[b]]$sigmam
          sigma <- rbind(sigma, cbind(zeroleft, temp, zeroright))
        }
        gnames = paste("gene",1:sum(nivec.list[[ss]]), sep = "")
        rownames(sigma) <- gnames
        colnames(sigma) <- gnames
        sigma.list[[ss]] <- sigma
      }
    } else {
      message("nivec.list and the selected network structure do NOT match ...\n")
      break()
    }

  } else if (structure =="Diff S, Identical W"){
    # for the part that structures are the same, weights are the same
    # assume the first ndiff rows have different structure.
    blklist <- list()
    sigma <- matrix(0, nrow = 0, ncol = sum(nivec.list[[1]]))
    nblk <- length(nivec.list[[1]])
    for (b in 1:nblk){
      ni <- nivec.list[[1]][b]
      blklist[[b]] <- generateBlki(ni=ni, ud=ud)
      zeroleft <- matrix(0, nrow = ni, ncol = sum(nivec.list[[1]][0:(b-1)]))
      zeroright <- matrix(0, nrow = ni, ncol = ifelse(b<nblk,sum(nivec.list[[1]][(b+1):nblk]),0))
      temp <- blklist[[b]]$sigmam
      sigma <- rbind(sigma, cbind(zeroleft, temp, zeroright))
    }
    gnames = paste("gene",1:sum(nivec.list[[1]]), sep = "")
    rownames(sigma) <- gnames
    colnames(sigma) <- gnames
    sigma.list[[1]] <- sigma

    for(ss in 2:length(nivec.list)){ # from the second matrix, only change the diffblk[[ss]] block part
      blklist <- list()
      sigma2 <- sigma
      nblk <- length(nivec.list[[ss]])
      temps <- matrix(0, nrow = 0, ncol = sum(nivec.list[[ss]]))
      for (b in diffblk[[ss]]){
        ni <- nivec.list[[ss]][b]
        blklist[[b]] <- generateBlki(ni=ni, ud=ud)
        zeroleft <- matrix(0, nrow = ni, ncol = sum(nivec.list[[ss]][0:(b-1)]))
        zeroright <- matrix(0, nrow = ni, ncol = ifelse(b<nblk,sum(nivec.list[[ss]][(b+1):nblk]),0))
        temp <- blklist[[b]]$sigmam
        temps <- cbind(zeroleft, temp, zeroright)
        diffid <- (sum(nivec.list[[ss]][0:(b-1)])+1) : sum(nivec.list[[ss]][0:b])
        sigma2[diffid,] <- temps
      }
      sigma.list[[ss]] <- sigma2
    }


  } else if (structure == "Diff S, Diff W"){
    for(ss in 1:length(nivec.list)){
      blklist <- list()
      sigma <- matrix(0, nrow = 0, ncol = sum(nivec.list[[ss]]))
      # sigma <- matrix(0, nrow = ni*nblk, ncol = ni*nblk)
      nblk <- length(nivec.list[[ss]])
      for (b in 1:nblk){
        ni <- nivec.list[[ss]][b]
        blklist[[b]] <- generateBlki(ni=ni, ud=ud)
        zeroleft <- matrix(0, nrow = ni, ncol = sum(nivec.list[[ss]][0:(b-1)]))
        zeroright <- matrix(0, nrow = ni, ncol = ifelse(b<nblk,sum(nivec.list[[ss]][(b+1):nblk]),0))
        temp <- blklist[[b]]$sigmam
        sigma <- rbind(sigma, cbind(zeroleft, temp, zeroright))
      }
      gnames = paste("gene",1:sum(nivec.list[[ss]]), sep = "")
      rownames(sigma) <- gnames
      colnames(sigma) <- gnames
      sigma.list[[ss]] <- sigma
    }
  } else {
    message("Please choose the joint network structure from: 'Identical S, Identical W', 'Identical S, Diff W', 'Diff S, Identical W', 'Diff S, Diff W'")
  }
  return(sigma.list)
}




check_ident_list <- function(xlist){
  if (length(unique(unlist(lapply(xlist, length)))) == 1){ #length of sublists are the same
    nsame = 0 # identical pairs of sublist
    for (cc in 1:length(xlist)){
      for (dd in 1:length(xlist)){
        if (cc < dd){
          nsame = nsame + all(xlist[[cc]] == xlist[[dd]])
        }
      }
    }
    if (nsame == choose(length(xlist), 2)){
      y = T
    } else {
      message("The elements in sublists are not all identical... \n")
      y = F
    }
  } else {
    message("The length of list elements are different... \n")
    y = F
  }
  return(y)
}

generateBlki <- function(ni, ud){
  mati <- matrix(0, ni, ni)
  for (i in 1:ni){
    for (j in i:ni){
      mati[i,j] <- ifelse(i!=j & runif(1) >0.6, sample(ud, 1),  #half chance: coexpression of i and j. prob in [-1,-0.6] U [0.6,1]
                          ifelse(i==j,1,0))
    }
  }
  mati0 = mati + t(mati) - diag(1, nrow = ni, ncol = ni)
  # eg = eigen(mati0)
  # str(eg)
  # delta = min(eg$values[eg$values>0])
  mati1 <- mati0 + diag(1, nrow = ni, ncol = ni)
  for (i in 1:20){
    # cat(i,".. \n")
    if (matrixcalc::is.positive.definite(mati1)){
      mati1 <- mati1
      break
    } else {
      mati1 <- mati1 + diag(1, nrow = ni, ncol = ni)
    }
  }

  binv = solve(mati1)
  bii <- diag(binv)
  sigmam <- binv
  for (i in 1:ni){
    for (j in 1:ni){
      sigmam[i,j] <- binv[i,j]/sqrt(bii[i]*bii[j])
    }
  }
  res <- list(sigmam = sigmam,
              Bm = mati1)
  return(res)
}


getCountList <- function(sigma.list, nvec = c(500, 500), ngene = NULL, a3 = 1.07, b3 = 2.28){
  if (is.null(ngene)){
    ngene = ncol(sigma.list[[1]])
  }
  count.list <- list()
  for (c in 1:length(sigma.list)){
    counti <- CountMap2(sigma = sigma.list[[c]], ngene = ngene, n=nvec[c], a3 = a3, b3 = b3)
    count.list[[c]] <- counti
  }
  return(count.list)
}

#' Generate a scRNA-seq count matrix.
#' See C0301, C0302 for function usages.
#' @return a list of covariance matrices.
#' @export
CountMap5 <- function(sigma, ngene, n, a1 = 3,
                      b1 = 1, a20 = 2,  b20 = 3, a30 = 1, b30 = 10){
  precision1 <- solve(sigma)
  mu <- rep(0, ngene)
  adj <- abs(sign(precision1))
  diag(adj) <- 0
  x <- mvtnorm::rmvnorm(n, mu, sigma)
  y <- x

  # -------------------------------------------------
  # generate count data from posterior distribution
  # one subject at a time
  # subject i
  z <- matrix(1, nrow = n, ncol = ngene)
  pij.vec <- rep(1, ngene)
  ycomplete <- x

  alphavec <- rep(0, ngene)
  betavec <- rep(0,ngene)
  scmeanvec <- rep(0,ngene)
  # -----------------------------
  for (j in 1:ngene) {
    dat <- x[, j]
    mu_v <- mean(dat)
    sd_v <- sd(dat)
    p_v <- pnorm(dat, mu_v, sd_v)

    # ---------------------
    # simulate sc from zip -- 1 gene
    # 3. alpha, beta, thetaij
    alphaj <- rgamma(1, shape = a20, rate = b20)
    betaj <- rgamma(1, shape = a30, rate = b30)
    scmeanj <- rgamma(1, shape = alphaj, rate = betaj)
    if (any(scmeanj >=1000) | any(scmeanj <1) ){
      scmeanj <- runif(1, 1,10)
    }
    sctemp <- rpois(n, scmeanj)
    alphavec[j] <- alphaj
    betavec[j] <- betaj
    scmeanvec[j] <- scmeanj

    ytemp <- quantile(sctemp, p_v)
    ycomplete[,j] <- ytemp
    pij <- rbeta(1, shape1 = a1, shape2 = b1)
    zijc <- sapply(1:n, function(x){
      px = ifelse(scmeanvec[j] > 10, 1, 1- ((1-pij) + pij * exp(-scmeanvec[j])))
      rbinom(1, size = 1, prob = px)
    })

    z[,j] <- zijc
    pij.vec[j] <- round(pij,3)
  }

  yout = ycomplete*z

  colnames(yout) <- paste("gene",1:ngene, sep = "")
  rownames(yout) <-  paste("cell",1:n,sep = "")
  colnames(ycomplete) <- paste("gene",1:ngene, sep = "")
  rownames(ycomplete) <- paste("cell",1:n,sep = "")
  colnames(z) <- paste("gene",1:ngene, sep = "")
  rownames(z) <- paste("cell",1:n,sep = "")

  result <- list()
  result$count <- yout
  result$count.nodrop <- ycomplete
  result$zijc <- z
  result$thetaj <- round(scmeanvec,2)
  result$pij <- pij.vec
  result$sigma <- sigma
  result$precision <- precision1
  return(result)
}


#' Generate a scRNA-seq count matrix, allowing the input of mean expression levels and dropout probabilities.
#' See C0301, C0302 for function usages.
#' @return a list of covariance matrices.
#' @export

CountMap6 <- function(sigma, ngene, n, a1 = 3,
                      b1 = 1, a20 = 2,  b20 = 3, a30 = 1, b30 = 10,
                      scmeanvec = NULL, pijvec = NULL, meanDiffPct = 0){
  if(meanDiffPct >0){
    ngdiff = ngene * meanDiffPct
    gdiff = sample(1:ngene, ngdiff)
  }

  precision1 <- solve(sigma)
  mu <- rep(0, ngene)
  adj <- abs(sign(precision1))
  diag(adj) <- 0
  x <- mvtnorm::rmvnorm(n, mu, sigma)
  y <- x

  # -------------------------------------------------
  # generate count data from posterior distribution
  # one subject at a time
  # subject i
  z <- matrix(1, nrow = n, ncol = ngene)
  pij.vec <- rep(1, ngene)
  ycomplete <- x

  alphavec <- rep(0, ngene)
  betavec <- rep(0,ngene)
  # -----------------------------
  if (is.null(scmeanvec)){
    scmeanvec <- rep(0,ngene)
    for (j in 1:ngene) {
      # i=1
      # j=1
      dat <- x[, j]
      mu_v <- mean(dat)
      sd_v <- sd(dat)
      p_v <- pnorm(dat, mu_v, sd_v)

      # ---------------------
      # simulate sc from zip -- 1 gene
      # 3. alpha, beta, thetaij
      alphaj <- rgamma(1, shape = a20, rate = b20)
      betaj <- rgamma(1, shape = a30, rate = b30)
      scmeanj <- rgamma(1, shape = alphaj, rate = betaj)
      if (any(scmeanj >=1000) | any(scmeanj <1) ){
        scmeanj <- runif(1, 1,10)
      }
      sctemp <- rpois(n, scmeanj)
      alphavec[j] <- alphaj
      betavec[j] <- betaj
      scmeanvec[j] <- scmeanj

      ytemp <- quantile(sctemp, p_v)
      ycomplete[,j] <- ytemp
      if(is.null(pijvec)){
        pij <- rbeta(1, shape1 = a1, shape2 = b1)
      } else {
        pij <- pijvec[j]
      }
      zijc <- sapply(1:n, function(x){
        px = ifelse(scmeanvec[j] > 10, 1, 1- ((1-pij) + pij * exp(-scmeanvec[j])))
        rbinom(1, size = 1, prob = px)
      })

      z[,j] <- zijc
      pij.vec[j] <- round(pij,3)
    }
  } else if (meanDiffPct > 0){
    for (j in 1:ngene) {
      dat <- x[, j]
      mu_v <- mean(dat)
      sd_v <- sd(dat)
      p_v <- pnorm(dat, mu_v, sd_v)

      # ---------------------
      # simulate sc from zip -- 1 gene
      # 3. alpha, beta, thetaij
      if (j %in% gdiff){
        scmeanj = scmeanvec[j] * runif(1, 0.5,2) # allow the mean to vary from 0.5 times to 2 times of the original mean
        scmeanvec[j] <- scmeanj

      } else {
        scmeanj <- scmeanvec[j]
      }

      sctemp <- rpois(n, scmeanj)

      ytemp <- quantile(sctemp, p_v)
      ycomplete[,j] <- ytemp
      if(is.null(pijvec)){
        pij <- rbeta(1, shape1 = a1, shape2 = b1)
      } else if (j %in% gdiff){
        pij <- rbeta(1, shape1 = a1, shape2 = b1)
      } else {
        pij <- pijvec[j]
      }
      zijc <- sapply(1:n, function(x){
        px = ifelse(scmeanvec[j] > 10, 1, 1- ((1-pij) + pij * exp(-scmeanvec[j])))
        rbinom(1, size = 1, prob = px)
      })

      z[,j] <- zijc
      pij.vec[j] <- round(pij,3)
    }
  } else {
    for (j in 1:ngene) {
      dat <- x[, j]
      mu_v <- mean(dat)
      sd_v <- sd(dat)
      p_v <- pnorm(dat, mu_v, sd_v)

      # ---------------------
      # simulate sc from zip -- 1 gene
      # 3. alpha, beta, thetaij
      scmeanj <- scmeanvec[j]
      sctemp <- rpois(n, scmeanj)

      ytemp <- quantile(sctemp, p_v)
      ycomplete[,j] <- ytemp
      if(is.null(pijvec)){
        pij <- rbeta(1, shape1 = a1, shape2 = b1)
      } else {
        pij <- pijvec[j]
      }
      zijc <- sapply(1:n, function(x){
        px = ifelse(scmeanvec[j] > 10, 1, 1- ((1-pij) + pij * exp(-scmeanvec[j])))
        rbinom(1, size = 1, prob = px)
      })

      z[,j] <- zijc
      pij.vec[j] <- round(pij,3)
    }
  }


  # -------------------------------

  yout = ycomplete*z

  colnames(yout) <- paste("gene",1:ngene, sep = "")
  rownames(yout) <-  paste("cell",1:n,sep = "")
  colnames(ycomplete) <- paste("gene",1:ngene, sep = "")
  rownames(ycomplete) <- paste("cell",1:n,sep = "")
  colnames(z) <- paste("gene",1:ngene, sep = "")
  rownames(z) <- paste("cell",1:n,sep = "")

  result <- list()
  result$count <- yout
  result$count.nodrop <- ycomplete
  result$zijc <- z
  result$thetaj <- round(scmeanvec,2)
  result$pij <- pij.vec
  result$sigma <- sigma
  result$precision <- precision1
  return(result)
}

