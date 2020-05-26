#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Gaussian kernel 
// [[Rcpp::export]]
vec GKernel(vec x){
  vec res(x.size());
  res=exp(-pow(x,2)/2);
  return(res);
}

// Ep kernel 
// [[Rcpp::export]]
vec EKernel(vec x){
  vec res(x.size());
  res=(1-pow(x,2))*3/4;
  res.elem(find(pow(x,2)>1)).zeros();
  return(res);
}





// [[Rcpp::export]]
mat VCM_cpp(mat X, vec U, vec Y, double h, vec nx, vec W_sample ,int kernel){
  int p = X.n_cols;
  mat res(2*p, nx.size());
  vec Uu(U.size());
  vec W_u(U.size());
  mat Uu_X(size(X)); 
  mat Gamma_u(size(X)[0],2*size(X)[1]);
  Gamma_u.cols(0,X.n_cols-1)=X;
  mat W_u_Gamma_u(size(Gamma_u));
  for(int k=0; k<nx.size(); k++){
    double u=nx[k];
    Uu=U-u;
    if(kernel==1) //Gaussian kernel
    {
      W_u=1/h*GKernel(Uu/h);
    }
    if(kernel==2) //Ep kernel
    {
      W_u=1/h*EKernel(Uu/h);
    }
    W_u = W_u%W_sample;
    Uu_X=X;
    Uu_X.each_col() %= Uu;  
    Gamma_u.cols(X.n_cols,2*X.n_cols-1)=Uu_X;
    // Gamma_u=join_cols(X,Uu_X);
    W_u_Gamma_u=Gamma_u;
    W_u_Gamma_u.each_col() %= W_u;
    res.col(k) = solve(Gamma_u.t()*W_u_Gamma_u,W_u_Gamma_u.t()*Y);
  }
  return (res);
}



// [[Rcpp::export]]
mat deleRows(mat x, int k1, int k2){
  x.shed_rows(k1-1,k2-1);
  return(x);
}



// [[Rcpp::export]]
double VCM_cv_cpp(mat X, vec U, vec Y, double h, vec W_sample, vec W_cv,IntegerVector start_ind, IntegerVector end_ind, int kernel){
  int p = X.n_cols;
  vec res(start_ind.size());
  for (int index=0; index<start_ind.size(); index++){
    vec u_grid(end_ind(index)-start_ind(index)+1);
    u_grid=U.rows(start_ind[index]-1,end_ind[index]-1);
    mat X_rest=deleRows(X,start_ind[index],end_ind[index]);
    vec U_rest=deleRows(U,start_ind[index],end_ind[index]);
    vec Y_rest=deleRows(Y,start_ind[index],end_ind[index]);
    vec W_sample_rest=deleRows(W_sample,start_ind[index],end_ind[index]);
    mat beta = VCM_cpp(X_rest, U_rest, Y_rest, h, u_grid, W_sample_rest, kernel);
    beta.shed_rows(p,2*p-1);
    mat X_keep=X.rows(start_ind[index]-1,end_ind[index]-1);
    vec Y_keep=Y.rows(start_ind[index]-1,end_ind[index]-1);
    vec W_cv_keep=W_cv.rows(start_ind[index]-1,end_ind[index]-1);
    res[index]=sum(W_cv_keep%pow(Y_keep-sum(X_keep%beta.t(),1),2));
  }
  return (mean(res));
}


