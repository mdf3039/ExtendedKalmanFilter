#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

MatrixXd CalculateJacobian(const VectorXd& x_state) {
    //obtain negations when needed
    float a = 1.0;
    if (x_state(0)<0){
        a = 1.0;
    }
    float b = 1.0;
    if (x_state(1)<0){
        b = 1.0;
    }
	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2)*a;
	float vy = x_state(3)*b;

	//pre-compute a set of terms to avoid repeated calculation
	float c1 = px*px+py*py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//check division by zero
	if((fabs(c1) < 0.00001) && (c1>=0)){
		c1 = 0.00001;
	}
	if((fabs(c1) < 0.00001) && (c1<0)){
		c1 = -0.00001;
	}

	//compute the Jacobian matrix
	Hj << (px/c2), (py/c2), 0, 0,
		  -(py/c1), (px/c1), 0, 0,
		  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

    /*if ((x_state(0)<0)&&(x_state(1)<0)){
        Hj << (px/c2), (py/c2), 0, 0,
            (py/c1), -(px/c1), 0, 0,
            py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
    }*/

	return Hj;
}

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  //measurement covariance
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  //measurement matrix
  ekf_.H_ = MatrixXd(2,4);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  //create a 4D state vector, we don't know yet the values of the x state
  ekf_.x_ = VectorXd(4);
  //state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  //the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  float noise_ax = 9.0;
  float noise_ay = 9.0;

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_ << 1,0,0,0,
        0,1,0,0;

  ekf_.F_ << 1,0,1,0,
        0,1,0,1,
        0,0,1,0,
        0,0,0,1;

  ekf_.P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;

  Hj_ << 1,1,0,0,
        -1,1,0,0,
        1,1,1,1;

  ekf_.H_ << 1,0,0,0,
            0,1,0,0;


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0]*sin(measurement_pack.raw_measurements_[1]),
            measurement_pack.raw_measurements_[0]*cos(measurement_pack.raw_measurements_[1]),
            0,0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0],
                measurement_pack.raw_measurements_[1], 0, 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  /*if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    return;
  }*/
  float noise_ax = 9;
  float noise_ay = 9;
  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  //set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
             0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
             dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
             0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    //ekf_.UpdateEKF(measurement_pack.raw_measurements_);

    //compute y, using h(x)
    VectorXd h_of_x;
    h_of_x = VectorXd(3);

    //calculate angle, given the quadrant X and Y values are in
    float pi = 3.14159265;
    float angle_r = atan(ekf_.x_[1]/ekf_.x_[0]);
    float range_rate = (ekf_.x_[0]*ekf_.x_[2]+ekf_.x_[1]*ekf_.x_[3])/sqrt(pow(ekf_.x_[0],2)+pow(ekf_.x_[1],2));
    if ((ekf_.x_[0]<0)&&(ekf_.x_[1]>0)){
        float angle_r = pi+atan(ekf_.x_[1]/ekf_.x_[0]);
        range_rate = (ekf_.x_[0]*-1.0*ekf_.x_[2]+ekf_.x_[1]*ekf_.x_[3])/sqrt(pow(ekf_.x_[0],2)+pow(ekf_.x_[1],2));
    }
    else if ((ekf_.x_[0]<0)&&(ekf_.x_[1]<0)){
        float angle_r = -1.0*pi+atan(ekf_.x_[1]/ekf_.x_[0]);
        range_rate = (ekf_.x_[0]*-1.0*ekf_.x_[2]+ekf_.x_[1]*-1.0*ekf_.x_[3])/sqrt(pow(ekf_.x_[0],2)+pow(ekf_.x_[1],2));
    }

    h_of_x << sqrt(pow(ekf_.x_[0],2)+pow(ekf_.x_[1],2)),
                angle_r,
                range_rate;
    //obtain the vector of raw measurements
    VectorXd z;
    z = VectorXd(3);
    z << measurement_pack.raw_measurements_[0],measurement_pack.raw_measurements_[1],
        measurement_pack.raw_measurements_[2];
    //make sure the angle is between -pi and pi
    while(z[1]<-3.14159265){
        z[1]+=2*3.14159265;
    }
    while(z[1]>3.14159265){
        z[1]-=2*3.14159265;
    }
    //if the difference between the two angles is more than 2.5 (arbitrary number close to pi),
    //add pi
    if abs(z[1]-h_of_x[1])>2.5{
        h_of_x[1] = h_of_x[1]+pi;
    }

    VectorXd y = z - h_of_x;
    //find K, using the H_jacobian
    Hj_ = CalculateJacobian(ekf_.x_);

    ekf_.UpdateEKF(y,Hj_,R_radar_);


  }
  else {
    // Laser updates
    //ekf_.Update(measurement_pack.raw_measurements_);
    VectorXd y = measurement_pack.raw_measurements_ - ekf_.H_*ekf_.x_;
    ekf_.Update(y,ekf_.H_,R_laser_);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
