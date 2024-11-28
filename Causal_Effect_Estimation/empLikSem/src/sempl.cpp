#include "sempl.h"

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double sempl_lrt(arma::vec weights, arma::mat y, arma::mat b_r, arma::mat o_r, double innerTol, int maxInnerIter, double aelWeight) {

  el_sem graph = el_sem(weights, y, b_r, o_r, innerTol, maxInnerIter, aelWeight);
  return graph.get_lrstat();

}


//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::List semplC(arma::vec weights, arma::mat y, arma::mat b_r, arma::mat o_r, double innerTol, int maxInnerIter, double aelWeight) {

  el_sem graph = el_sem(weights, y, b_r, o_r, innerTol, maxInnerIter, aelWeight);
  arma::mat B = mat(size(b_r), fill::zeros);
  B((find(b_r == 1))) = weights;

  arma::mat resid = graph.get_resid();
  resid.each_col() %= sqrt(graph.get_weights());
  arma::mat Omega = resid.t() * resid;

  return Rcpp::List::create(Rcpp::Named("lrt") = graph.get_lrstat(),
                            Rcpp::Named("weights") = graph.get_weights(),
                            Rcpp::Named("df") = graph.get_df(),
                            Rcpp::Named("dual") = graph.get_dual(),
                            Rcpp::Named("B") = B,
                            Rcpp::Named("Omega") = Omega);

}

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
arma::vec sempl_grad(arma::vec weights, arma::mat y, arma::mat b_r, arma::mat o_r, double innerTol, int maxInnerIter, double aelWeight) {

  el_sem graph = el_sem(weights, y, b_r, o_r, innerTol, maxInnerIter, aelWeight);
  return graph.get_grad();

}


//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double naive_sempl_lrt(arma::vec weights, arma::mat y, arma::mat b_r, arma::mat o_r, double innerTol, int maxInnerIter, double aelWeight) {

  el_sem_naive graph = el_sem_naive(weights, y, b_r, o_r, innerTol, maxInnerIter, aelWeight);
  return graph.get_lrstat();

}


//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
arma::vec naive_sempl_grad(arma::vec weights, arma::mat y, arma::mat b_r, arma::mat o_r, double innerTol, int maxInnerIter, double aelWeight) {

  el_sem_naive graph = el_sem_naive(weights, y, b_r, o_r, innerTol, maxInnerIter, aelWeight);
  return graph.get_grad();

}


//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::List naive_semplC(arma::vec weights, arma::mat y, arma::mat b_r, arma::mat o_r, double innerTol, int maxInnerIter, double aelWeight) {

  el_sem_naive graph = el_sem_naive(weights, y, b_r, o_r, innerTol, maxInnerIter, aelWeight);
  arma::mat B = mat(size(b_r), fill::zeros);
  arma::mat Omega = mat(size(b_r), fill::zeros);

  B((find(b_r == 1))) = weights.subvec(0, accu(b_r) - 1);
  Omega(find(trimatl(o_r) == 1)) = weights.subvec(accu(b_r), weights.n_elem - 1);
  Omega = symmatl(Omega);

  return Rcpp::List::create(Rcpp::Named("lrt") = graph.get_lrstat(),
                            Rcpp::Named("weights") = graph.get_weights(),
                            Rcpp::Named("df") = graph.get_df(),
                            Rcpp::Named("dual") = graph.get_dual(),
                            Rcpp::Named("B") = B,
                            Rcpp::Named("Omega") = Omega);


}

