#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() : n_laser_(2), n_radar_(3) {
  is_initialized_ = false;
  
  time_us_ = 0;
  
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  //set state dimention
  n_x_ = 5;
  
  //set augumented dimention
  n_aug_ = n_x_ + 2;
  
  // initial state vector
  x_ = VectorXd(n_x_);
  
  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;
  
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 2;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  
  //define spreading parameter
  lambda_ = 3 - n_aug_;
  
  // set sigma points dimention
  n_sigma_ = 2 * n_aug_ + 1;
  
  // Initialize predicted sigma points matrix.
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_);
  
  // Initialize weights of sigma points
  weights_ = VectorXd(n_sigma_);
  
  // set weights
  double weight_0 = lambda_/(lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i < n_sigma_; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_ + lambda_);
    weights_(i) = weight;
  }
  
  // Initialize measurement noise covariance matrix for laser
  R_laser_ = MatrixXd(n_laser_, n_laser_);
  R_laser_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;
  
  // Initialize measurement noise covariance matrix for radar
  R_radar_ = MatrixXd(n_radar_, n_radar_);
  R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0,std_radrd_*std_radrd_;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float px0 = meas_package.raw_measurements_(0) * cos(meas_package.raw_measurements_(1));
      float py0 = meas_package.raw_measurements_(0) * sin(meas_package.raw_measurements_(1));
      
      // Since vertical velocity is not measured by radar, 
      // aproximate initial velocity and yaw by the velocity in the rho direction.
      float v0 = meas_package.raw_measurements_(2);
      float yaw0 = meas_package.raw_measurements_(1);
      
      x_ << px0, py0, v0, yaw0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0, 0, 0;
    }
    
    time_us_ = meas_package.timestamp_;
    
    // done initializing, no need to predict or update
    is_initialized_ = true;
    
    return;
  }
  
  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  if((use_laser_ == true & meas_package.sensor_type_ == MeasurementPackage::LASER) or 
     (use_radar_ == true & meas_package.sensor_type_ == MeasurementPackage::RADAR))
  {
    float dt = (meas_package.timestamp_ -time_us_) / 1.0e6;  //dt - expressed in seconds
    time_us_ = meas_package.timestamp_;
    
    if(dt>1e-5)
    {
      // Predict
      UKF::Prediction(dt);
    }
    
    /*****************************************************************************
     *  Update
     ****************************************************************************/
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Radar updates
      UKF::Update(meas_package, n_radar_);
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // Lidar updates
      UKF::Update(meas_package, n_laser_);
    }
    
    // print the output
    cout << "x_ = " << x_ << endl;
    cout << "P_ = " << P_ << endl;     
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
    //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();
  
  //set first column of sigma point matrix
  Xsig.col(0)  = x_;
  
  //set remaining sigma points
  for (int i = 0; i < n_x_; i++)
  {
    Xsig.col(i+1)     = x_ + sqrt(lambda_+n_x_) * A.col(i);
    Xsig.col(i+1+n_x_) = x_ - sqrt(lambda_+n_x_) * A.col(i);
  }
  
    //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma_);
  
  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  MatrixXd L_n = sqrt(lambda_ + n_aug_) * L;
  
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + L_n.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - L_n.col(i);
  }
  
  //predict sigma points
  for (int i = 0; i< n_sigma_; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;
    
    // yawrate after delta_t 
    double yaw_dt = yaw + yawd*delta_t;
    
    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw_dt) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw_dt) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw_dt;
    double yawd_p = yawd;
    double p_a = 0.5 * nu_a * delta_t * delta_t;

    //add noise
    px_p = px_p + p_a * cos(yaw);
    py_p = py_p + p_a * sin(yaw);
    v_p = v_p + nu_a*delta_t;
    
    double nu_dt = nu_yawdd * delta_t;
    yaw_p = yaw_p + 0.5 * nu_dt * delta_t;
    yawd_p = yawd_p + nu_dt;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  
  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < n_sigma_; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
  //angle normalization
  tools.NormAngle(x_(3));
  
  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_sigma_; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    
    //angle normalization
    tools.NormAngle(x_diff(3));
    
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
  
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */

void UKF::Update(MeasurementPackage meas_package, const int &n_z) {
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sigma_);
  
  VectorXd z = VectorXd(n_z);
  cout << n_z <<"\n";
  //transform sigma points into measurement space
  if (n_z == 3) {   //radar sigma points
      for (int i = 0; i < n_sigma_; i++) {  //2n+1 simga points
        // extract values for better readibility
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v  = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);
        
        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;
    
        // measurement model
        Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);            //r
        if (p_x != 0) {
          Zsig(1,i) = atan2(p_y,p_x);                   //phi
        } else {
          Zsig(1,i) = M_PI / 2;
        }
        
        if (fabs(Zsig(0,i)) > 1e-5)
        {
          Zsig(2,i) = (p_x*v1 + p_y*v2 ) / Zsig(0,i);   //r_dot
        } else {
          Zsig(2,i) = v;
        }
      }
      
      z << meas_package.raw_measurements_(0),
           meas_package.raw_measurements_(1),
           meas_package.raw_measurements_(2);
  
  } else if (n_z == 2) {    //laser sigma points
    //transform sigma points into measurement space
    for (int i = 0; i < n_sigma_; i++) {  //2n+1 simga points
      // extract values for better readibility
      double p_x = Xsig_pred_(0,i);
      double p_y = Xsig_pred_(1,i);
      
      // measurement model
      Zsig(0,i) = p_x;                        //px
      Zsig(1,i) = p_y;                        //py
    }
    
    z << meas_package.raw_measurements_(0),
         meas_package.raw_measurements_(1);
  }
  //mean predicted measurement
  VectorXd z_pred = tools.CalZpred(n_z, n_sigma_, weights_, Zsig);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  
  if (n_z == 3) {
   tools.CalSandTc(S, Tc, n_sigma_, Zsig, z_pred, Xsig_pred_, x_,
                        weights_, R_radar_);
  } else if (n_z == 2) {
    tools.CalSandTc(S, Tc, n_sigma_, Zsig, z_pred, Xsig_pred_, x_,
                        weights_, R_laser_);
  }
  
  MatrixXd Si = S.inverse();
  
  //Kalman gain K;
  MatrixXd K = Tc * Si;
  
  //residual
  VectorXd z_diff_p = z - z_pred;

  //angle normalization
  tools.NormAngle(z_diff_p(1));
  
  //update state mean and covariance matrix
  x_ = x_ + K * z_diff_p;
  P_ = P_ - K*S*K.transpose();
  
  //angle normalization
  tools.NormAngle(x_(3));
  if (n_z == 3) {
    NIS_radar_ = z_diff_p.transpose() * Si * z_diff_p;
  } else if (n_z == 2) {
    NIS_laser_ = z_diff_p.transpose() * Si * z_diff_p;
  }
  
}