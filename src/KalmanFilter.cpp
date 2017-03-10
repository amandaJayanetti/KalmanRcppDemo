#include <RcppEigen.h>

using namespace std;

//' Kalman Filter
//'
//' @param Y    Best estimate of measurements/state at times k,k+dt,k+2dt,..
//' @param X    Initial state
//' @param F    Prediction matrix used to obtain next state
//' @param H    Matrix to model the signal
//' @param Q    Covariance matrix to model uncertainity associated with untracked influences 
//' @param R    Covariance matrix to model uncertainity eg. sensor noise
//' @param P    Error covariance matrix of the state estimate 
//' 
//' @examples
//' \dontrun{
//' Y<-matrix(c(0.39,0.50,0.48,0.29,0.25,0.32,0.34,0.48,0.41,0.45),10,1,byrow=TRUE)
//' X<-matrix(c(0.39,0,7),3,1,byrow=TRUE)
//' F<-matrix(c(1,0.05,0,0,1,0.05,0,0,1),3,3,byrow = TRUE)
//' H<-matrix(c(1,1,0),1,3,byrow = TRUE)
//' Q<-matrix(c(0.0025, 0.0025, 0, 0, .0025, 0, 0, 0.0025, 0),3,3,byrow=TRUE)
//' R<-matrix(c(3),1,1,byrow=TRUE)
//' P<-matrix(c(0.05,0.05,5,0.5,10000,10,0.5,10,100),3,3,byrow = TRUE)
//' kalmanUpdate(Y,X,F,H,Q,R,P)
//' }
//'
//' @return    updated measurements/state at times k,k+dt,k+2dt,..
 // [[Rcpp::export]]
Eigen::MatrixXd kalmanUpdate(Eigen::MatrixXd Y, Eigen::MatrixXd X,Eigen::MatrixXd F,Eigen::MatrixXd H,Eigen::MatrixXd Q, Eigen::MatrixXd R, Eigen::MatrixXd P){
  int rows=Y.rows();
  int cols=X.rows();
  Eigen::MatrixXd X_,K,Result(rows,cols);
  
  try{
    for(int i=0;i<Y.size();i++){
      Eigen::VectorXd y=Y.row(i);
      X_=F*X;
      P=F*P*F.transpose() + Q;
      K=P*H.transpose()*(H*P*H.transpose()+R).inverse(); 
      X_=X_ + K*(y-(H*X_));
      X=X_;
      P=P-(K*H*P);
      Eigen::VectorXd B(Eigen::Map<Eigen::VectorXd>(X_.data(), X_.cols()*X_.rows()));
      Result.row(i)=B.transpose(); 
    }
  }
  catch(...){
    Rcpp::stop("Unknown Exception");
  }
  return Result;
}
 

