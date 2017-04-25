#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);
  
  // Normalize Angle
  void NormAngle(double &x);
  
  // calculate 
  Eigen::VectorXd CalZpred(const int &n_z, const int &n_sigma, const Eigen::VectorXd &weights, const Eigen::MatrixXd & Zsig);
  
  void CalSandTc(Eigen::MatrixXd &S, Eigen::MatrixXd &Tc, const int &n_sigma, const Eigen::MatrixXd &Zsig,
                 const Eigen::VectorXd &z_pred, const Eigen::MatrixXd &Xsig_pred,
                 const Eigen::VectorXd &x, const Eigen::VectorXd &weights,
                 const Eigen::MatrixXd &R);

};

#endif /* TOOLS_H_ */
