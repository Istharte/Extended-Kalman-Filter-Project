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
