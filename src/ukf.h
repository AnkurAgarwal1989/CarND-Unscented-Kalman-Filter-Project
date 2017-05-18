#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* augmented sigma points matrix (holds large augmented data)
  MatrixXd Xsig_aug_;

  ///* predicted sigma points matrix (sigma points after prediction)
  MatrixXd Xsig_pred_;

  ///* predicted sigma points matrix for radar (sigma points after prediction)
  MatrixXd Zsig_radar_pred_;

  ///* Lidar Measurement Function
  MatrixXd H_laser_;

  ///* Radar measuerement vector [rho, phi, rho_dot] (predicted)
  VectorXd z_radar_pred_;

  ///* Laser measuerement vector [px, py] (predicted)
  VectorXd z_laser_pred_;

  ///* Radar covariance matrix
  MatrixXd S_radar_;

  ///* Radar Kalman Gain
  MatrixXd K_radar_;

  ///* Radar Noise
  MatrixXd R_radar_;

  ///* Laser Noise
  MatrixXd R_laser_;


  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Radar meas dimension
  int n_z_radar_;

  ///* Laser meas dimension
  int n_z_laser_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  //////////////////////////////////////////////////UKF Filter Steps

  /**
  * Initialize Filter
  */
  void Init(MeasurementPackage &meas_package);

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
  * Prediction Predicts sigma points, the state, and the state covariance
  * matrix
  * @param delta_t Time between k and k+1 in s
  */
  void Prediction(double delta_t);

  /**
  * Updates the state and the state covariance matrix using a laser measurement
  * @param meas_package The measurement at k+1
  */
  void UpdateLidar(VectorXd z);

  /**
  * Updates the state and the state covariance matrix using a radar measurement
  * @param meas_package The measurement at k+1
  */
  void UpdateRadar(VectorXd z);
  /////////////////////////////////////////////////UKF Filter Steps


  //////////////////////////////////////////////UKF Filter Implementation

  /**
   * NormalizeAngle
   */
  void NormalizeAngle(double& angle);

  /**
   * GenerateSigmaPoints
  */
  void GenerateSigmaPoints();

  /**
  * PredictSigmaPoints
  */
  void PredictSigmaPoints(double dt);

  /**
  * CalculateWeights
  */
  void CalculateWeights();

  /**
  * Calculate Mean and Covariance for the prediction step
  */
  void CalculateMeanCovariance_State();

  /**
  * Calculate Mean and Covariance for the Radar measurements
  */
  void CalculateMeanCovariance_Radar();

  /**
  * Calculate Kalman Gain for Radar
  */
  void CalculateKalmanGain_Radar();
  
  /**
  * A helper method to calculate NIS. (Normalized Innovation Squared)
  */
  double CalculateNIS(const VectorXd &diff_meas, const MatrixXd &S);
  /////////////////////////////////////////////////UKF Filter Implementation
};

#endif /* UKF_H */
