#include <math.h>
#include "kalman_filter.h"
#include <iostream> //.for std

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;


//debug
#define USERDEBUG

#ifdef USERDEBUG
#define Debug(x) cout << x
#else
#define Debug(x)
#endif


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

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */

  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
  //P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  //debug info
  Debug( "[kalman_filter]: Update for RADAR " << endl);
  Debug( "[kalman_filter]: x_ = " << x_ << endl);

  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  float rho = sqrt(px * px + py * py);
  float phi = atan2(py, px);
  float rho_dot = 0.0f;

  if (fabs(rho) >= 0.000001f){
    rho_dot = (px * vx + py * vy) / rho;
  }

  Debug( "[kalman_filter]: rho_dot = " << rho_dot << endl);

  VectorXd z_pred = VectorXd(3);
  z_pred << rho, phi, rho_dot;

  Debug( "[kalman_filter]: z_pred = " << z_pred << endl);

  VectorXd y = z - z_pred;
  Debug( "[kalman_filter]: y = " << y << endl);
  Debug( "[kalman_filter]: M_PI = " << M_PI << endl);

  while((y(1) > M_PI) || (y(1) < -M_PI)){
    if(y(1) > M_PI){
      y(1) -= 2.0f * M_PI;
    }
    else{
      y(1) += 2.0f * M_PI;
    }
  }

  Debug( "[kalman_filter]: y = " << y << endl);

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  Debug( "[kalman_filter]: K = " << K << endl);

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

  Debug( "[kalman_filter]: x_ = " << x_ << endl);
  Debug( "[kalman_filter]: P_ = " << P_ << endl);
}
