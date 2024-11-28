#ifndef EL_SEM_NAIVE_H
#define EL_SEM_NAIVE_H


#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>
#include "el_sem_settings.h"

using namespace Rcpp ;
using namespace arma;

class el_sem_naive
{
public:
  // Initialize sem object:
  el_sem_naive(arma::vec weights, arma::mat y, arma::mat b_r, arma::mat o_r, double tol, int maxIter, double aelWeight);

  int update_dual();
  vec get_weights();
  vec get_dual();
  vec get_grad();
  int update_BO(arma::vec weights);
  double get_lrstat();
  double get_convCrit();
  int get_lower_ind(int i);
  int get_df();


protected:
  vec d_; //vector of denominators of p_n (ie p_n = 1 / d_(n))
  vec dual_; // dual variables
  double tol_;
  int maxIter_;
  double lr_;
  mat b_;
  mat o_;
  mat constraints_;
  mat o_ind_;
  double aelWeight_;

  uvec b_ones_;
  uvec o_ones_;
  uvec b_rows_, b_cols_, o_rows_, o_cols_;


  int v_;
  int n_;
  int degree_restrict_;
  mat y_;
  int counter_;
  double conv_crit_;

  void set_gradient_hessian(vec &grad, mat &hessian);
  int backtracking(vec update);
  double step_size_;



private:
};
#endif // ELSEM_H

