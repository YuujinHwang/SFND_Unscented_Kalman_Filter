#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 2;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  is_initialized_ = false;

  time_us_ = 0.0;

  n_x_ = 5;

  n_aug_ = n_x_ + 2;

  lambda_ = 3 - n_aug_;

  weights_ = VectorXd(2 * n_aug_ + 1);

  weights_.fill(0.5/(lambda_ + n_aug_));
  weights_(0) = lambda_/(lambda_ + n_aug_);

  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);


}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if(!is_initialized_)
  {
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      double rho = meas_package.raw_measurements_[0]; //range
      double phi = meas_package.raw_measurements_[1]; //angle
      double v = meas_package.raw_measurements_[2];   //velocity
      double x = rho * cos(phi);
      double y = rho * sin(phi);
      double vx = v * cos(phi);
      double vy = v * sin(phi);
      x_ << x, y, v, phi, 0;
      P_ << 0.1, 0, 0, 0, 0,
            0, 0.1, 0, 0, 0,
            0, 0, 0.1, 0, 0,
            0, 0, 0, 0.1, 0,
            0, 0, 0, 0, 1; 
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      double x = meas_package.raw_measurements_[0];
      double y = meas_package.raw_measurements_[1];
      double phi = atan2(y,x);

      x_ << x, y, 0, phi, 0;
      P_ << 0.1, 0, 0, 0, 0,
            0, 0.1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 0.1, 0,
            0, 0, 0, 0, 1; 
    }
  std::cout << "initialized : "<< std::endl << x_ << std::endl;
  is_initialized_ = true;
  return;
  }
  double dt  = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;
  // std::cout << dt << std::endl;

  Prediction(dt);

  if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }
  else if(meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    UpdateLidar(meas_package);
  }

}



void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  // create argumented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  // create argumented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  MatrixXd L = P_aug.llt().matrixL();
  Xsig_aug.col(0) = x_aug;

  for(int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  // predict sigma points

  for(int i = 0; i < 2*n_aug_+1; i++)
  {
    x_aug = Xsig_aug.col(i);
    double px, py, vel, psi, psid;

    // model integration
    if (fabs(x_aug(4)) < 0.001)
    {
      px = x_aug(0) + x_aug(2)*delta_t*cos(x_aug(3));
      py = x_aug(1) + x_aug(2)*delta_t*sin(x_aug(3));
    }
    else
    {
      px = x_aug(0) + x_aug(2) * (sin(x_aug(3) + x_aug(4)*delta_t) - sin(x_aug(3))) / x_aug(4);
      py = x_aug(1) + x_aug(2) * (cos(x_aug(3)) - cos(x_aug(3) + x_aug(4)*delta_t)) / x_aug(4);
    }
    
    // add noise
    px = px + x_aug(5)*delta_t*delta_t*cos(x_aug(3))/2.0;
    py = py + x_aug(5)*delta_t*delta_t*sin(x_aug(3))/2.0;

    vel = x_aug(2) + x_aug(5)*delta_t;

    psi = x_aug(3) + x_aug(4)*delta_t + x_aug(6)*delta_t*delta_t/2.0;
    psid = x_aug(4) + x_aug(6)*delta_t;

    Xsig_pred_(0,i) = px;
    Xsig_pred_(1,i) = py;
    Xsig_pred_(2,i) = vel;
    Xsig_pred_(3,i) = psi;
    Xsig_pred_(4,i) = psid;
  }

  // set state
  x_.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++)
  {
    x_ = x_ + weights_(i)*Xsig_pred_.col(i);
  }
  // std::cout<< Xsig_pred_ << std::endl;
  // std::cout << x_ << std::endl;
  // std::cout << weights_ << std::endl;
  // set covariance matrix
  P_.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while(x_diff(3) <-M_PI) x_diff(4) += 2.*M_PI;

    P_ = P_ + weights_(i)*x_diff*x_diff.transpose();
  }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

  VectorXd z_ = meas_package.raw_measurements_;

  int n_z_ = 2;

  MatrixXd Zsig = MatrixXd(n_z_, 2*n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z_);
  MatrixXd S = MatrixXd(n_z_, n_z_);
  for(int i = 0; i < 2*n_aug_+1; i++)
  {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double pz = Xsig_pred_(2,i);

    Zsig(0,i) = px;
    Zsig(1,i) = py;
  }

  z_pred.fill(0.0);

  for(int i = 0; i < 2*n_aug_+1; i++)
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  S.fill(0.0);

  for(int i = 0; i < 2*n_aug_+1; i++)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S = S + weights_(i)*z_diff*z_diff.transpose(); 
  }

  MatrixXd R = MatrixXd(n_z_, n_z_);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;

  S = S + R;

  MatrixXd Tc = MatrixXd(n_x_, n_z_);
   Tc.fill(0.0);
   for(int i = 0; i < 2*n_aug_+1; i++){
     VectorXd x_diff = Xsig_pred_.col(i) - x_;

     VectorXd z_diff = Zsig.col(i) - z_pred;

     Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
   }

   MatrixXd K = Tc * S.inverse();

   VectorXd z_diff = z_ - z_pred;

   x_ = x_ + K*z_diff;

   P_ = P_ - K*S*K.transpose();
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  VectorXd z_ = meas_package.raw_measurements_;

  int n_z_ = 3;

  MatrixXd Zsig = MatrixXd(n_z_, 2*n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z_);
  MatrixXd S = MatrixXd(n_z_, n_z_);

  for(int i = 0; i < 2*n_aug_+1; i++)
  {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double psi = Xsig_pred_(3,i);

    double vx = v*cos(psi);
    double vy = v*sin(psi);

    Zsig(0,i) = sqrt(px*px + py*py);
    Zsig(1,i) = atan2(py, px);
    Zsig(2,i) = (px*vx + py*vy)/Zsig(0,i);
  }

  z_pred.fill(0.0);

  for(int i = 0; i < 2*n_aug_+1; i++)
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  S.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    while(z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while(z_diff(1) <-M_PI) z_diff(1) += 2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R  = MatrixXd(n_z_, n_z_);
  R << std_radr_*std_radr_, 0, 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;

  S = S + R;

  MatrixXd Tc = MatrixXd(n_x_, n_z_);
   Tc.fill(0.0);
   for(int i = 0; i < 2*n_aug_+1; i++){
     VectorXd z_diff = Zsig.col(i) - z_pred;
     while(z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
     while(z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

     VectorXd x_diff = Xsig_pred_.col(i) - x_;
     while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
     while(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

     Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
   }

   MatrixXd K = Tc * S.inverse();

   VectorXd z_diff = z_ - z_pred;
   while(z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
   while(z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

   x_ = x_ + K*z_diff;

   P_ = P_ - K*S*K.transpose();
}