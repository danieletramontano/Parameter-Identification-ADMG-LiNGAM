#include "elSEM_naive.h"


el_sem_naive::el_sem_naive(arma::vec weights, arma::mat y, arma::mat b_r, arma::mat o_r, double tol, int maxIter, double aelWeight)
{
  y_ = y; // n_ by v_ matrix where each row is an observation, each column is a variable
  v_ = y.n_cols; //number of variables

  aelWeight_ = aelWeight;
  if(aelWeight_ > 0){
    n_ = y.n_rows + 1;
  } else {
    n_ = y.n_rows;
  }

  counter_ = 0; //number of iterations
  tol_ = tol; // tol for the inner minimization
  maxIter_ = maxIter; // max iterations for the inner minimization
  conv_crit_ = 1.0; //convergence critiera for inner maximization (norm of gradient)

  //Find the free elements of B and the constrained elements of Omega
  b_ones_ = find(b_r == 1);
  o_ones_ = find(trimatl(o_r) == 1);
  b_ = mat(v_, v_, fill::zeros);
  o_ = mat(v_, v_, fill::zeros);


  //Check
  b_rows_ = b_ones_ - floor(b_ones_ / v_) * v_;
  b_cols_ = floor(b_ones_ / v_);

  o_rows_ = o_ones_ - floor(o_ones_ / v_) * v_;
  o_cols_ = floor(o_ones_ / v_);

  constraints_ = arma::mat(n_, v_ + v_ * (v_ + 1) / 2, fill::zeros);
  dual_ = vec(constraints_.n_cols, fill::zeros);
  d_ = vec(n_, fill::ones); // the denominator (1 + lambda * g(x_i, B))



  update_BO(weights);
}

int el_sem_naive::update_BO(arma::vec weights)
{
  counter_ = 0; //number of iterations
  conv_crit_ = 1.0; //convergence critiera for inner maximization (norm of gradient)
  constraints_.zeros();
  dual_.zeros();
  d_.ones();

  b_(b_ones_) = weights.subvec(0, b_ones_.n_elem - 1);
  o_(o_ones_) = weights.subvec(b_ones_.n_elem, weights.n_elem - 1);

  arma::mat resid = y_ -  y_ * b_.t();

  // Filling in the constraint matrix
  constraints_.zeros();


  int z;
  int k;
  int m;


  // mean constraints
  if(aelWeight_ > 0){
    constraints_.cols(0, v_ - 1).rows(0, n_ - 2) = resid;

    // covariance constraints
    int count = v_;
    for(k = 0; k < v_ ; k++) {
      for(m = k; m < v_ ; m++) {
        constraints_.col(count++).rows(0, n_ - 2) = resid.col(m) % resid.col(k) - o_(m, k);
      }
    }
    arma::rowvec u = mean(constraints_.rows(0, n_ - 2), 0);

    constraints_.row(n_ - 1) = -aelWeight_ * u;

  } else { // no adjustment

  // mean constraints
  constraints_.cols(0, v_ - 1) = resid;

  // covariance constraints
  // k indexes columns
  // m indexes rows

  int count = v_;
  for(k = 0; k < v_ ; k++) {
    for(m = k; m < v_ ; m++) {
      constraints_.col(count++) = resid.col(m) % resid.col(k) - o_(m, k);
    }
  }
  }

  update_dual();
  return 0;
}


int el_sem_naive::update_dual()
{
  dual_.zeros();
  d_.ones();
  vec grad(constraints_.n_cols, fill::zeros);
  vec update(constraints_.n_cols, fill::zeros);
  mat hessian(constraints_.n_cols, constraints_.n_cols, fill::zeros);

  // backtracking parameters
  int back_tracking_counter;
  int max_back_track = 20; // steps required to scale update by 1e-8

  conv_crit_ = 1.0;
  while(conv_crit_ > tol_) {

    // build in backtracking if necessary

    set_gradient_hessian(grad, hessian);

    update = solve(hessian, grad);


    back_tracking_counter = 0;
    while(!backtracking(update) && back_tracking_counter < max_back_track) {
      update *= BACKTRACKING_SCALING_CONST;
      if(++back_tracking_counter > max_back_track) {
        // if back tracking does not return because of max back tracking iterations
        return INFEASIBLE_RETURN;
      }
    }
    // else
    // update the dual
    dual_ -= update;

    // update the convergence criteria
    conv_crit_ = norm(grad ,2);

    // if max iterations hit, give infeasible return
    // In practice, this shouldn't be used often because
    // we have a strictly convex problem, but condition is built in for error checking
    if(counter_++ > maxIter_) {
      return INFEASIBLE_RETURN;
    }

  }

  d_ = (constraints_ * dual_) + 1.0;
  return 0;
}


// check whether update remains in feasible space
int el_sem_naive::backtracking(vec update)
{
  arma::vec new_d = (constraints_ * (dual_ - update)) + 1.0;
  if(all( new_d > (1.0 / n_)) & (2 * accu(log((new_d))) > get_lrstat() )){
    return 1;
  }
  return 0;

}


// set gradient and hessian (wrt to dual variables) to appropriate values
void el_sem_naive::set_gradient_hessian(vec &grad, mat &hessian)
{

  d_ = (constraints_ * dual_) + 1.0;
  int k;

  for(k = 0; k < constraints_.n_cols; k++){
    grad(k) = -accu(constraints_.col(k) / d_);
  }

  // Rcout <<"Set Grad" <<std::endl;
  hessian.zeros();
  int i;
  for(i = 0; i < n_ ; i ++) {
    hessian += constraints_.row(i).t() * constraints_.row(i) / pow(d_(i), 2);
  }

  // If the hessian is poorly conditioned take a gradient step instead
  if(rcond(hessian) < 1e-10){
    hessian = eye(size(hessian));
  }

}

//return scaled d
arma::vec el_sem_naive::get_weights()
{
  return 1.0 / (d_ * n_);
}

//return scaled d
arma::vec el_sem_naive::get_dual()
{
  return dual_;
}

double el_sem_naive::get_convCrit(){
  return conv_crit_;
}


//return scaled d
double el_sem_naive::get_lrstat()
{
  return 2 * accu(log((d_)));
}




arma::vec el_sem_naive::get_grad()
{
  // jacobian of the constraints wrt to the parameters
  arma::mat dg_dTheta(constraints_.n_cols, b_ones_.n_elem + o_ones_.n_elem, fill::zeros);
  arma::vec dL_dTheta(b_ones_.n_elem + o_ones_.n_elem, fill::zeros);
  arma::vec mod_weights =  (1.0 / d_.rows(0, n_ - 2 )) - 1.0 / d_(n_ - 1) * aelWeight_ / (n_ - 1.0);



  int u, v, k, m, n;
  uvec v_parents, u_parents;


  // Gradient for Mean Constraints
  for(v = 0; v < v_; v++){
    v_parents = find(b_rows_ == v);
    for(k = 0; k < v_parents.n_elem; k++){
      if(aelWeight_ > 0 ){
        dg_dTheta(v, v_parents(k)) = -accu(y_.col(b_cols_(v_parents(k))).rows(0, n_ - 2) % mod_weights);
      } else {
        dg_dTheta(v, v_parents(k)) = -accu(y_.col(b_cols_(v_parents(k))) / d_);
      }


    }
  }

  vec weighted_residu, weighted_residv;

  int count = v_;
  arma::uvec o_ind;
  // u indexes columns
  // v indexes rows
  for(u = 0; u < v_ ; u++ ) {
    for(v = u; v < v_ ; v++ ) {
      u_parents = find(b_rows_ == u);
      v_parents = find(b_rows_ == v);

      if(aelWeight_ > 0 ){
        weighted_residu = constraints_.col(u).rows(0, n_ - 2) %  mod_weights;
        weighted_residv = constraints_.col(v).rows(0, n_ - 2) %  mod_weights;
      } else {
        weighted_residu = constraints_.col(u) / d_;
        weighted_residv = constraints_.col(v) / d_;
      }
      for(m = 0; m < v_parents.n_elem; m++){
        if(aelWeight_ > 0){
          dg_dTheta(count, v_parents(m)) += -accu(y_.col(b_cols_(v_parents(m))).rows(0, n_ - 2) % weighted_residu);
        } else {
          dg_dTheta(count, v_parents(m)) += -accu(y_.col(b_cols_(v_parents(m))) % weighted_residu);
        }

      }

      for(m = 0; m < u_parents.n_elem; m++){
        if(aelWeight_ > 0){
          dg_dTheta(count, u_parents(m)) += -accu(y_.col(b_cols_(u_parents(m))).rows(0, n_ - 2) % weighted_residv);
        } else {
          dg_dTheta(count, u_parents(m)) += -accu(y_.col(b_cols_(u_parents(m))) % weighted_residv);
        }

      }

      if(o_(v, u) != 0){
        o_ind = find((o_rows_ == v) % (o_cols_ == u)) ;
        // Possible change required here
        dg_dTheta(count, b_rows_.n_elem +  o_ind(0) ) = -n_;
      }

      count++;
    }
  }



  dL_dTheta = -(dual_.t() * dg_dTheta).t();

  return -2 * dL_dTheta;
}

int el_sem_naive::get_df(){
  return constraints_.n_cols - b_rows_.n_elem - o_rows_.n_elem;
}


