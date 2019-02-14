#include "kalman_filter.h"
#include<iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
    // KF Prediction step
    x_ = F_*x_;// + u_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
    // KF Measurement update step
    MatrixXd y_ = z - H_ * x_;
    MatrixXd S_ = H_ * P_ * H_.transpose() + R_;
    MatrixXd K_ = P_ * H_.transpose() * S_.inverse();
    MatrixXd I = MatrixXd::Identity(4, 4);

    // new state
    x_ = x_ + K_ * y_;
    P_ = (I - K_ * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
    // KF Measurement update step
    double px = x_(0);
    double py = x_(1);
    double vx = x_(2);
    double vy = x_(3);

    MatrixXd tmp_z(3, 1);
    tmp_z << sqrt(px*px + py*py),  atan2(py, px), (px*vx + py*vy) / sqrt(px*px + py*py);

    MatrixXd y_ = z - tmp_z;
    MatrixXd S_ = H_ * P_ * H_.transpose() + R_;
    MatrixXd K_ = P_ * H_.transpose() * S_.inverse();
    MatrixXd I = MatrixXd::Identity(4, 4);

    if (y_(1) < -M_PI)y_(1) += 2*M_PI;
    if (M_PI < y_(1))y_(1) -= 2*M_PI;
    // new state
    x_ = x_ + K_ * y_;
    P_ = (I - K_ * H_) * P_;
}
