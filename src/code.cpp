// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' rowmeans cpp function
//'
//' @param X a matrix
//' @export
// [[Rcpp::export]]
arma::vec arma_rowMeans(const arma::mat X){
  int       nRows = X.n_rows;
  arma::vec out(nRows);
  for(int i = 0; i < nRows; i++){
    out(i) = mean(X.row(i));
  }

  return(out);
}

//' colmeans cpp function
//'
//' @param X a matrix.
//' @export
// [[Rcpp::export]]
arma::vec arma_colMeans(const arma::mat X){
  int       nCols = X.n_cols;
  arma::vec out(nCols);
  for(int i = 0; i < nCols; i++){
    out(i) = mean(X.col(i));
  }
  return(out);
}

//' calculate the Mn statistic from Cai et al. (2013)
//'
//' @param X samples by covariates (genes)
//' @param Y samples by covariates (genes)
//' @param alpha significance level
//' @export
// [[Rcpp::export]]
List testCov_cpp(arma::mat X, arma::mat Y, double alpha = 0.05){
  int p = X.n_cols, n1 = X.n_rows, n2 = Y.n_rows;
  // tmp and scalev
  double pi     = 3.1415926,
    //qalpha = -log(8*pi) - 2*log(log(1/(1-alpha))),
    //cri    = 4*log(p) - log(log(p)) + qalpha,
    //atmp, btmp,
    CLX, CLX_rev, CLX_p;
  arma::mat Sx, Sy, xa(n1, p), ya(n2, p), Vx, Vy ; //, vx1, vx2, vx3, numo, deno, Tnm, ts1, ts2, ts1_, ts2_;
  arma::vec Xcolm, Ycolm; //, g1(n1), g2(n2), ts(J, fill::zeros);

  Sx = cov(X)*(n1-1)/n1;
  Sy = cov(Y)*(n2-1)/n2;
  Xcolm = arma_colMeans(X);
  Ycolm = arma_colMeans(Y);
  for (int i=0; i < n1; i++){
    xa.row(i) = X.row(i) - Xcolm.t();
  }
  for (int i=0; i < n2; i++){
    ya.row(i) = Y.row(i) - Ycolm.t();
  }

  Vx = (pow(xa.t(),2) * pow(xa,2) / n1) - ((xa.t() * xa) % Sx * 2/n1) + pow(Sx,2);
  Vy = (pow(ya.t(),2) * pow(ya,2) / n2) - ((ya.t() * ya) % Sy * 2/n2) + pow(Sy,2);

  CLX = (pow((Sx - Sy),2) / (Vx/n1 + Vy/n2)).max();
  CLX_rev = CLX - (4*log(p) - log(log(p)));
  CLX_p = 1 - exp(-exp(-CLX_rev/2) / (sqrt(8*pi)));

  return List::create(
    _["pvalue"] = CLX_p,
    _["CLX"] = CLX,
    _["CLX_type1"] = CLX_rev
  );
}

//' Permutation to learn data-specific Gumbel distribution parameters
//'
//' @param X: samples by covariates (genes)
//' @param searchx: candidate change point locations that will be tested.
//' @param start: the start location of the current segment
//' @param nperm: number of permutations
//' @export
// [[Rcpp::export]]
arma::mat CLXstatPerm(arma::mat X, arma::vec searchx, int start, unsigned int nperm){ //
  arma::mat Xperm = X;
  arma::mat X1, X2, statmat(nperm*searchx.size(), 3, fill::zeros);
  List stat_perm;
  arma::uvec permv;
  int t1;
  // X should be cells by genes
  // first permute the matrix
  for(unsigned int j=0; j < nperm; j++){
    permv = randperm(X.n_rows);
    Xperm = X.rows(permv);
    for (unsigned int i =0; i < searchx.size(); i++){
      t1 = searchx(i) - 1;
      X1 = Xperm.rows(0 , t1 - start + 1);
      X2 = Xperm.rows(t1 - start + 2 , X.n_rows - 1);
      stat_perm = testCov_cpp(X1, X2);
      statmat(i + j*searchx.size(),0) = stat_perm[1]; // return the CLX
      statmat(i + j*searchx.size(),1) = t1+1;
      statmat(i + j*searchx.size(),2) = j;
    }
  }

  // return List::create(
  //   _["permv"] = permv,
  //   _["xperm"] = Xperm,
  //   _["t1"] = t1,
  //   _["X1"] = X1,
  //   _["X2"] = X2,
  //   _["statperm"] = stat_perm,
  //   _["statmat"] = statmat
  // );
  return statmat;
}


// [[Rcpp::export]]
arma::uvec calc_ranks(const arma::vec& da_ta) {
  return (arma::sort_index(arma::sort_index(da_ta)) + 1);
}  // end calc_ranks


//' Hierarchical Permutation to learn data-specific Gumbel distribution parameters
//'
//' @param X: samples by covariates (genes)
//' @param searchx: candidate change point locations that will be tested.
//' @param start: the start location of the current segment
//' @param nperm: number of permutations
//' @param nt2: number of varied segment 2. Default is 4, which means we fix the length of segment 1 and vary the length of segment 2 for 4 times. The lengths for the varied segment 2s are equally spaced.
//' @param minband: minimum bandwidth of a segment. Default is 20.
//' @export
// [[Rcpp::export]]
arma::mat CLXPermHier(arma::mat X, arma::vec searchx, int start,unsigned int nt2 = 4, int minband = 20, unsigned int nperm = 100){
  arma::mat Xperm = X;
  arma::mat X1, X2, statmat(nperm*searchx.size()*(nt2+1), 4, fill::zeros);
  List stat_perm;
  arma::uvec permv;
  unsigned int t1, t2;
  // X should be cells by genes
  // first permute the matrix
  for(unsigned int j=0; j < nperm; j++){
    permv = randperm(X.n_rows);
    Xperm = X.rows(permv);
    for (unsigned int i =0; i < searchx.size(); i++){
      t1 = searchx(i) - 1;
      for(unsigned int s=0; s < nt2 + 1; s++){
        t2 = t1 + minband + std::max((int) round((X.n_rows + start - t1 - minband )/nt2), 2) * s;
        if(t2 > X.n_rows-1){
          t2 = X.n_rows-1;
        };
        X1 = Xperm.rows(0 , t1 - start + 1);
        X2 = Xperm.rows(t1 - start + 2 , t2 - start + 1);
        stat_perm = testCov_cpp(X1, X2);
        statmat(s + i*(nt2+1) + j*searchx.size()*(nt2+1),0) = stat_perm[1]; // return the CLX statistic
        statmat(s + i*(nt2+1) + j*searchx.size()*(nt2+1),1) = t1+1;
        statmat(s + i*(nt2+1) + j*searchx.size()*(nt2+1),2) = j;
        statmat(s + i*(nt2+1) + j*searchx.size()*(nt2+1),3) = t2 + 1;
        // std::cout << j << i << s <<"\t"<< t2<< std::endl;
      } // end s
    } // end i
  } // end j

  // return List::create(
  //   _["permv"] = permv,
  //   _["xperm"] = Xperm,
  //   _["t1"] = t1,
  //   _["X1"] = X1,
  //   _["X2"] = X2,
  //   _["statperm"] = stat_perm,
  //   _["statmat"] = statmat
  // );
  return statmat;
}

// [[Rcpp::export]]
double myQnorm(double mp){
  double out = R::qnorm5(mp, 0.0, 1.0, 1, 0);
  return out;
}

// [[Rcpp::export]]
arma::mat getMatNPN(arma::mat input, int N)
{
  arma::vec newX;
  arma::vec tempCol(N);
  for(int icol =0; icol< (int)input.n_cols; icol++){
    newX = sort(input.col(icol));
    std::map<int, int> ranks;

    int rank = 1;
    for(int index = 0; index < N; index++)
    {
      int element = newX[index];
      // Update rank of element
      if (ranks[element] == 0)
      {
        ranks[element] = rank;
        rank++;
      }
    }
    // Assign ranks to elements
    tempCol = input.col(icol);
    for(int index = 0; index < N; index++)
    {
      //int element = input.col(icol).row(index);
      //input.col(icol).row(index) = ranks[input.col(icol).row(index)];
      tempCol[index] = ranks[tempCol[index]];
    }

    arma::vec res(N);
    for(int i = 0; i < N; i++) {
      double mp = tempCol[i]/(N+1);
      res[i] = myQnorm(mp);
    }
    double msd = stddev(res);
    for (int i =0; i<N; i++){
      res[i] = res[i] / msd;
    }

    input.col(icol) = res;

    tempCol.clear();
    ranks.clear();
    res.clear();
  }


  return input;

}

