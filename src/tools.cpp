#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  VectorXd res(4);
  
  // check the validity of the estimations and ground_truth
  if((estimations.size() == 0) or (estimations.size() != ground_truth.size()))
  {
    std::cout << "Invalid estimations or ground_truth data" << std::endl;
    return rmse;
  }
  // accumulate squared residuals
  for(unsigned int i=0; i < estimations.size(); ++i){
    res = estimations[i] - ground_truth[i];
    //coefficient-wise multiplication
    res = res.array()*res.array();
    rmse += res;
  }
  
  // calculate the mean
  rmse = rmse / estimations.size();
  
  // calculate the squared root
  rmse = rmse.array().sqrt();
  
  return rmse;

}

// Normalise Angle
void Tools::NormAngle(double &x)
{
  while (x> M_PI) x-=2.*M_PI;
  while (x<-M_PI) x+=2.*M_PI;
}

// calculate mean predicted measurement
VectorXd Tools::CalZpred(const int &n_z, const int &n_sigma, const VectorXd &weights, const MatrixXd & Zsig)
{
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < n_sigma; i++) {
      z_pred = z_pred + weights(i) * Zsig.col(i);
  }
  return z_pred;
}

// Calculate S and Tc
void Tools::CalSandTc(MatrixXd &S, MatrixXd &Tc, const int &n_sigma, const MatrixXd &Zsig,
                      const VectorXd &z_pred, const MatrixXd &Xsig_pred, const VectorXd &x,
                      const VectorXd &weights, const MatrixXd &R)
{
  //calculate cross correlation matrix and measurement covariance matrix
  S.fill(0.0);
  Tc.fill(0.0);
  for (int i = 0; i < n_sigma; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    
    //angle normalization
    Tools::NormAngle(z_diff(1));
    
    MatrixXd z_diff_t = z_diff.transpose();
    
    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    
    //angle normalization
    Tools::NormAngle(x_diff(3));
    
    S = S + weights(i) * z_diff * z_diff_t;
    
    Tc = Tc + weights(i) * x_diff * z_diff_t;
    
  }

  //add measurement noise covariance matrix
  S = S + R;
  
}