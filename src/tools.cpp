#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * Calculate the RSME of each value eastimated in the state vector.
   */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
  if (estimations.size() != ground_truth.size()) {
    std::cout << "Estimate and ground truth list size did not match in RSME calculation." << std::endl;
  }
  else {
    // TODO: accumulate squared residuals
    for (int i = 0; i < estimations.size(); ++i) {
      if (estimations[i].size() == 0 || ground_truth[i].size() == 0) {
        std::cout << "States were empty in RSME calculation." << std::endl;
      }
      else if (estimations[i].size() != ground_truth[i].size()) {
        std::cout << "Unequal state sizes encountered in RMSE calculation." << std::endl;
      }
      else {
        VectorXd diff = estimations[i] - ground_truth[i];
        VectorXd sqdiff = diff.array() * diff.array();
        rmse = rmse + sqdiff;
      }
    }
  }
  // TODO: calculate the mean
  rmse = rmse.array() / estimations.size();
  // TODO: calculate the squared root
  rmse = rmse.array().sqrt();
  // return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3, 4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // check division by zero
  if ((px == 0) || (py == 0)) {
    std::cout << "Divide by zero encountered in Jacobian calculation." << std::endl;
    return Hj;
  }
  else {
    float P = px * px + py * py;
    Hj << px / sqrt(P), py / sqrt(P), 0, 0,
          -py / P, px / P, 0, 0,
          py * (vx*py - vy*px) / pow(P, 1.5), px * (vy*px - vx*py) / pow(P, 1.5), px / sqrt(P), py / sqrt(P);
    return Hj;
  }
}
