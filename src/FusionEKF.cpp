#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  /* init ekf_.P_ */
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;


  /* init H_laser_ */
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}


void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  static double previous_x, previous_y = 0;
  /**
   * Initialization
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR){
    std::cout << "RADAR" << std::endl;
  } else {
    std::cout << "LASER" << std::endl;
  }
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */
    std::cout <<"initializing " << std::endl;
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      double rhodot = measurement_pack.raw_measurements_[2];
    if (phi < -M_PI)phi += 2*M_PI;
    if (M_PI < phi)phi -= 2*M_PI;
      double px = rho * cos(phi); // rho * cos(phi)
      double py = rho * sin(phi);
      double vx = rhodot * cos(phi);
      double vy = rhodot * sin(phi);
      ekf_.x_ << px, py, vx, vy;
      previous_timestamp_ = measurement_pack.timestamp_;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      ekf_.x_ << measurement_pack.raw_measurements_[0], 
              measurement_pack.raw_measurements_[1], 
              0, 
              0;
      previous_timestamp_ = measurement_pack.timestamp_;
    }

    previous_x = ekf_.x_(0);
    previous_y = ekf_.x_(1);
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  double noise_ax = 9, noise_ay = 9;
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //laser
  /* update F */
  MatrixXd tmp_f = MatrixXd::Identity(4,4);
  tmp_f(0, 2) = dt;
  tmp_f(1, 3) = dt;
  ekf_.F_ = tmp_f;
  /* update Q */
  MatrixXd tmp_q(4,4);
  tmp_q <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
       0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
       dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
       0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
  ekf_.Q_ = tmp_q;

  ekf_.Predict();
  

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    std::cout <<"radar" << std::endl;
    double rho = measurement_pack.raw_measurements_[0];
    double phi = measurement_pack.raw_measurements_[1];
    double rhodot = measurement_pack.raw_measurements_[2];
    if (phi < -M_PI)phi += 2*M_PI;
    if (M_PI < phi)phi -= 2*M_PI;
    assert(-M_PI < phi && phi < M_PI);
    double px = rho * cos(phi); // rho * cos(phi)
    double py = rho * sin(phi);
    double vx = ekf_.x_(2);//0;//(px - previous_x) / dt;
    double vy = ekf_.x_(3);//0;//(py - previous_y) / dt;
    std::cout << "vx = " << px <<" " << previous_x <<" " << dt<<std::endl;
    std::cout << "vy = " << py <<" " << previous_y <<" " << dt << std::endl;
    std::cout << "px,py,vx,vy = " << px <<" " << py <<" " << vx <<" " << vy << endl;
    VectorXd tmp_z(4);
    tmp_z << px, py, vx, vy;
    Hj_ = tools.CalculateJacobian(tmp_z);
    ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, Hj_, R_radar_, ekf_.Q_);
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // TODO: Laser updates
    ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, H_laser_, R_laser_, ekf_.Q_);
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  previous_x = ekf_.x_(0);
  previous_y = ekf_.x_(1);
  previous_timestamp_ = measurement_pack.timestamp_;

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
