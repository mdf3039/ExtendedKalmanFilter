#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

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

void KalmanFilter::Predict(const VectorXd &acc, const float dt) {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  //the x_ is obtained; add in the effect of acc to the four components
  x_[0] = x_[0]+0.5*acc[0]*dt*dt;
  x_[1] = x_[1]+0.5*acc[1]*dt*dt;
  x_[2] = x_[2]+acc[0]*dt;
  x_[3] = x_[3]+acc[1]*dt;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &y, const MatrixXd &H, const MatrixXd &R) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P_ * Ht + R;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &y_e, const MatrixXd &Hj, const MatrixXd &R_e) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  MatrixXd Hj_t = Hj.transpose();
  MatrixXd S = Hj * P_ * Hj_t + R_e;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Hj_t;
  MatrixXd K = PHt * Si;
  //new estimate
  x_ = x_ + (K * y_e);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj) * P_;

}
