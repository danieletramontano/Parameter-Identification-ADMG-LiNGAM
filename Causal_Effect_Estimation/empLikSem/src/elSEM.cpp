#include "elSEM.h"



el_sem::el_sem(arma::vec weights, arma::mat y, arma::mat b_r, arma::mat o_r, double tol, int maxIter, double aelWeight)
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


  b_ = mat(v_, v_);

  //Find the free elements of B and the constrained elements of Omega
  b_ones_ = find(b_r != 0);


  arma::mat temp(v_, v_, fill::ones);

  o_zeros_ = find((o_r + trimatu(temp))  == 0 );


  //Check
  b_rows_ = b_ones_ - floor(b_ones_ / v_) * v_;
  b_cols_ = floor(b_ones_ / v_);

  o_rows_ = o_zeros_ - floor(o_zeros_ / v_) * v_;
  o_cols_ = floor(o_zeros_ / v_);


  constraints_ = arma::mat(n_, v_ + o_zeros_.n_elem, fill::zeros);
  dual_ = vec(constraints_.n_cols, fill::zeros);
  d_ = vec(n_, fill::ones); // the denominator (1 + lambda * g(x_i, B))

  update_B(weights);
}

int el_sem::update_B(arma::vec b_weights)
{
  counter_ = 0; //number of iterations
  conv_crit_ = 1.0; //convergence critiera for inner maximization (norm of gradient)
  constraints_.zeros();
  dual_.zeros();
  d_.ones();

  b_.zeros();
  b_(b_ones_) = b_weights;
  arma::mat resid = y_ -  y_ * b_.t();


  // Filling in the constraint matrix
  constraints_.zeros();

  int count = v_;
  int z;
  int k;
  int m;


  if(aelWeight_ > 0){
    // mean constraints
    constraints_.cols(0, v_ - 1).rows(0, n_ - 2) = resid;

  // covariance constraints
  for(k = 0; k < o_zeros_.n_elem; k++) {
      constraints_.col(k + v_).rows(0, n_ - 2) = resid.col(o_rows_(k)) % resid.col(o_cols_(k));
  }

  constraints_.row(n_ - 1) = -aelWeight_ * mean(constraints_.rows(0, n_ - 2), 0);

  } else {
    // No adjusted points
    constraints_.cols(0, v_ - 1) = resid;

    // covariance constraints
    for(k = 0; k < o_zeros_.n_elem; k++) {
      constraints_.col(k + v_) = resid.col(o_rows_(k)) % resid.col(o_cols_(k));
    }
  }

  update_dual();
  return 0;
}


int el_sem::update_dual()
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
    
    double epsilon = 1e-5;
    
    // Create a diagonal matrix for regularization
    MatrixXd regularization = epsilon * MatrixXd::Identity(hessian.rows(), hessian.cols());
    hessian = hessian + regularization
    cout << "Hello World!"
    
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
int el_sem::backtracking(vec update)
{
  arma::vec new_d = (constraints_ * (dual_ - update)) + 1.0;
  if(all( new_d > (1.0 / n_)) & (2 * accu(log((new_d))) > get_lrstat() )){
    return 1;
  }
  return 0;

}


// set gradient and hessian (wrt to dual variables) to appropriate values
void el_sem::set_gradient_hessian(vec &grad, mat &hessian)
{

  d_ = (constraints_ * dual_) + 1.0;
  int k;

  for(k = 0; k < constraints_.n_cols; k++){
    grad(k) = -accu(constraints_.col(k) / d_);
  }

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
arma::vec el_sem::get_weights()
{
  return 1.0 / (d_ * n_);
}

//return scaled d
arma::vec el_sem::get_dual()
{
  return dual_;
}

double el_sem::get_convCrit(){
  return conv_crit_;
}


//return scaled d
double el_sem::get_lrstat()
{
  return 2 * accu(log((d_)));
}




arma::vec el_sem::get_grad()
{
  // jacobian of the constraints wrt to the parameters
  arma::mat dg_dTheta(constraints_.n_cols, b_ones_.n_elem, fill::zeros);
  arma::vec dL_dTheta(b_ones_.n_elem, fill::zeros);
  arma::vec mod_weights =  (1.0 / d_.rows(0, n_ - 2 )) - 1.0 / d_(n_ - 1) * aelWeight_ / (n_ - 1.0);

  int v, k, m, n;
  double weighted_mean;

  uvec v_parents, u_parents;

  // Gradient for Mean Constraints
  for(v = 0; v < v_; v++){

    v_parents = find(b_rows_ == v);
    for(k = 0; k < v_parents.n_elem; k++){
      if(aelWeight_ > 0){
        dg_dTheta(v, v_parents(k)) = -accu(y_.col(b_cols_(v_parents(k))).rows(0, n_ - 2) % mod_weights);
      } else {
        dg_dTheta(v, v_parents(k)) = -accu(y_.col(b_cols_(v_parents(k))) / d_);
      }
    }
  }


  //cov restrictions
  arma::vec weighted_residu, weighted_residv;
  for(k = 0; k < o_zeros_.n_elem; k++){
    u_parents = find(b_rows_ == o_rows_(k));
    v_parents = find(b_rows_ == o_cols_(k));

    if(aelWeight_ > 0 ){
      weighted_residu = constraints_.col(o_rows_(k)).rows(0, n_ - 2) % mod_weights;
      weighted_residv = constraints_.col(o_cols_(k)).rows(0, n_ - 2) % mod_weights;
    } else {
      weighted_residu = constraints_.col(o_rows_(k)) / d_;
      weighted_residv = constraints_.col(o_cols_(k)) / d_;
    }


    for(m = 0; m < v_parents.n_elem; m++){
        if(aelWeight_ > 0 ){
          dg_dTheta(v_ + k, v_parents(m)) = -accu(y_.col(b_cols_(v_parents(m))).rows(0, n_ - 2) % weighted_residu);
        } else {
          dg_dTheta(v_ + k, v_parents(m)) = -accu(y_.col(b_cols_(v_parents(m))) % weighted_residu);
        }
    }

    for(m = 0; m < u_parents.n_elem; m++){
      if(aelWeight_ > 0 ){
        dg_dTheta(v_ + k, u_parents(m)) = -accu(y_.col(b_cols_(u_parents(m))).rows(0, n_ - 2) % weighted_residv);
      } else {
        dg_dTheta(v_ + k, u_parents(m)) = -accu(y_.col(b_cols_(u_parents(m))) % weighted_residv);
      }
    }
  }


  dL_dTheta = -(dual_.t() * dg_dTheta).t();

  return -2 * dL_dTheta;
}

int el_sem::get_df() {
  return constraints_.n_cols - b_rows_.n_elem;
}


arma::mat el_sem::get_resid() {
  return constraints_.cols(0, v_ - 1);
}

arma::mat el_sem::get_constraints() {
  return constraints_;
}
