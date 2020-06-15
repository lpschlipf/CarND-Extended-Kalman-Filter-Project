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
  // Observation Transformations
  H_laser_<< 1, 0, 0, 0,
             0, 1, 0, 0;
  Hj_ << 0, 0, 0, 0,
         0, 0, 0, 0,
         0, 0, 0, 0;

  noise_ax = 9.0;
  noise_ay = 9.0;

  // Set up all the other matrices
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.Q_ = MatrixXd(4, 4);
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * Initialize state vector based on first measurement.
     * To get some rough initial speed eastimate, the assumptions are that the car initially started in the 
     * coordinate origion at time 0 and thus the radial speed can be used.
     */

    // first measurement
    cout << "EKF: ...Initialization" << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = measurement_pack.raw_measurements_(0);
      double phi = measurement_pack.raw_measurements_(1);
      double rhodot = measurement_pack.raw_measurements_(2);
      ekf_.x_ << rho * cos(phi),
                 rho * sin(phi),
                 5.2,
                 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_ << measurement_pack.raw_measurements_(0),
                 measurement_pack.raw_measurements_(1),
                 5.2,
                 0;
    }
    // Initial covariance Matrix
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1000, 0,
               0, 0, 0, 1000;

    // Set timestamp correctly, i.e. we start now at zero.
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */
  // Calculate elapsed time. Convert to seconds.
  float dt = (float)(measurement_pack.timestamp_ - previous_timestamp_) * 0.000001;
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;

  ekf_.Q_ << dt * dt * dt * dt * noise_ax / 4, 0, dt * dt * dt * noise_ax / 2, 0,
             0, dt * dt * dt * dt * noise_ay / 4, 0, dt * dt * dt * noise_ay / 2,
             dt * dt * dt * noise_ax / 2, 0, dt * dt * noise_ax, 0,
             0, dt * dt * dt * noise_ay / 2, 0, dt * dt * noise_ay;

  ekf_.Predict();

  /**
   * Update
   * Perform a Kalman or an Extended Kalman update depending on the Sensor that is used.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }
}
