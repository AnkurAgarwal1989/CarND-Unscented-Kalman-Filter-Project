#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

#define EPS 1e-6

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  //init check
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  //state sizes
  n_x_ = 5;
  n_aug_ = 7;

  n_z_radar_ = 3;
  n_z_laser_ = 2;

  // Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  //Lidar Matrices
  z_laser_pred_ = VectorXd(n_z_laser_);
  H_laser_ = MatrixXd(n_z_laser_, n_x_);
  R_laser_ = MatrixXd(n_z_laser_, n_z_laser_);

  //Radar matrices
  z_radar_pred_ = VectorXd(n_z_radar_);
  S_radar_ = MatrixXd(n_z_radar_, n_z_radar_);
  R_radar_ = MatrixXd(n_z_radar_, n_z_radar_);

  //Predicted Sigma Points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Zsig_radar_pred_ = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);

  //Augmented Sigma Points matrix
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //Vector of weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  //Calculate and stpre the weights
  CalculateWeights();
  //cout << weights_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.2;  //1.0      //data1:0.3 //data2: 0.15

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0; //10.0  //data1:0.8 //data2: 0.01

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
}

UKF::~UKF() {}

/************************************************************************
UKF Filter Implementation
Begin
/************************************************************************/

/**
 * Generates a matrix of augmeneted sigma points from the current state.
 * The process matrix will be augmented to add noise.
 * 
*/
void UKF::GenerateSigmaPoints() {
  //Create a new augmented Process matrix. Q is acceleration noise variance
  /*|P  0|
    |0  Q|*/
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_*std_a_;
  P_aug(6, 6) = std_yawdd_*std_yawdd_;

  //Augment the state
  VectorXd x_aug = VectorXd(n_aug_);

  //State becomes first sigma point
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //Sqrt of matrix
  MatrixXd A = P_aug.llt().matrixL(); //A is an n_aug_ x n_aug_ matrix
  double temp = sqrt(lambda_ + n_aug_);
  Xsig_aug_.col(0) = x_aug;
  //For every column in A matrix, we create 2 sigma points
  for (int i = 0; i < n_aug_; ++i) {
    VectorXd sig = temp*A.col(i);
    Xsig_aug_.col(i + 1) = x_aug + sig;
    Xsig_aug_.col(i + 1 + n_aug_) = x_aug - sig;
    //cout << "Generating sigma POints" << endl;
    //cout << Xsig_aug_ << endl;
  }
}

/*
 * Given a set of augmented sigma points, predict sigma points
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 * These sigma points will be reused by the Update steps (radar measurement update)
 */
void UKF::PredictSigmaPoints(double dt) {
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    double px = Xsig_aug_(0, i);
    double py = Xsig_aug_(1, i);
    double v = Xsig_aug_(2, i);
    double phi = Xsig_aug_(3, i);
    double phi_d = Xsig_aug_(4, i);
    double nu_a = Xsig_aug_(5, i);
    double nu_phi_dd = Xsig_aug_(6, i);

    //Our state update formula varies if rho_d ==0
    if (abs(phi_d) < 0.001) {
      px += v*cos(phi)*dt;
      py += v*sin(phi)*dt;
    }
    else {
      px += (v / phi_d) * (sin(phi + phi_d*dt) - sin(phi));
      py += (v / phi_d) * (-cos(phi + phi_d*dt) + cos(phi));
    }
    phi += phi_d*dt;

    //Add Noise
    double dt2 = dt*dt;
    px += 0.5*(dt2)*cos(phi)*nu_a;
    py += 0.5*(dt2)*sin(phi)*nu_a;
    v += dt*nu_a;
    phi += 0.5*(dt2)*cos(phi)*nu_a;
    phi_d += dt*nu_phi_dd;

    //Wrte to predicted sigma matrix
    Xsig_pred_(0, i) = px;
    Xsig_pred_(1, i) = py;
    Xsig_pred_(2, i) = v;
    Xsig_pred_(3, i) = phi;
    Xsig_pred_(4, i) = phi_d;
  }

}

/**
* Calculate Weight Vector
* if n_aug and lambda values don't change, we need to call this only once.
* Each sigma point has some weihht associated with it.
* w0 = lambda/(n_sigma + lambda)
* W_others = 1/2*(n_sigma + lambda)

*/
void UKF::CalculateWeights() {
  //Number of sigma points = 2*n_aug + 1
  //Caclulate weight vector
  double denom = n_aug_ + lambda_; //n_aug_ + lambda_;
  weights_(0) = lambda_ / denom;
  double temp = 0.5 / denom;
  for (int i = 1; i < 2 * n_aug_ + 1; ++i) {
    weights_(i) = temp;
  }
}

/**
* Calculates Mean and Covariance from the predicted sigma points.
* Note: for the Prediction Step
* Calculates x and P
*/
void UKF::CalculateMeanCovariance_State() {
  //Compute new state vector (x) and covariance matrix (P)
  //We can reuse previous memory
  //Calculate Mean State Vector (n_x_)
  x_.fill(0.0);
  P_.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; ++i) {
    x_ += weights_(i)*Xsig_pred_.col(i);
  }
  //Calculate Covariance Matrix (n_x_ , n_x)
  //variance is squared diff from mean.
  VectorXd x_diff_temp = VectorXd(n_x_);
  for (int i = 0; i < 2*n_aug_ + 1; ++i) {
    x_diff_temp = Xsig_pred_.col(i) - x_;
    //Normalize angle [-pi, +pi]
    //As simple trick is to use atan2(as it is always returns [-pi, +pi]
    //phi is ind 3
    //angle normalization
    while (x_diff_temp(3)> M_PI) x_diff_temp(3) -= 2.*M_PI;
    while (x_diff_temp(3)<-M_PI) x_diff_temp(3) += 2.*M_PI;
    P_ += weights_(i) * x_diff_temp * x_diff_temp.transpose();
  }
  //cout << "Mean State " << endl;
  //cout << x_ << endl;

  //cout << "Covariance " << endl;
  //cout << P_ << endl;

}

/**
* Calculates measurement Mean and Covariance from the predicted sigma points.
* Note: for the Radar Measurement Step
* Calculates z_radar and S_radar
*/
void UKF::CalculateMeanCovariance_Radar()
{
  //Use Xsig_pred_
  //Radar sigma points
  z_radar_pred_.fill(0.0);
  S_radar_.fill(0.0);

  for (int i = 0; i < 2*n_aug_ + 1; ++i) {
    //Radar measurement function
    //rho =sqrt(px^2+py^2)
    //phi = atan2(py, px)
    //phi_dot = px*vx + py*vy/(rho)
    // measurement model
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double phi = Xsig_pred_(3, i);
    double v_x = v*cos(phi);
    double v_y= v*sin(phi);

    Zsig_radar_pred_(0, i) = sqrt(p_x*p_x + p_y*p_y); //rho
    Zsig_radar_pred_(1, i) = atan2(p_y, p_x); //phi
    Zsig_radar_pred_(2, i) = (p_x*v_x + p_y*v_y) / max(sqrt(p_x*p_x + p_y*p_y), EPS);   //rho_dot                     
    
    z_radar_pred_ += weights_(i)*Zsig_radar_pred_.col(i);
  }
  //cout << "Z Pred " << endl;
  //cout << z_radar_pred_ << endl;
  //Calculate Covariance Matrix (n_z_ , n_z_)
  //variance is squared diff from mean.
  VectorXd z_diff_temp = VectorXd(n_z_radar_);
  for (int i = 0; i < 2*n_aug_  + 1; ++i) {
    z_diff_temp = Zsig_radar_pred_.col(i) - z_radar_pred_;
    //angle normalization
    while (z_diff_temp(1)> M_PI) z_diff_temp(1) -= 2.*M_PI;
    while (z_diff_temp(1)<-M_PI) z_diff_temp(1) += 2.*M_PI;
    S_radar_ += weights_(i) * z_diff_temp * z_diff_temp.transpose();
  }

  S_radar_ += R_radar_;
}

/**
* Calculate Kalman Gain for the Radar. This is diffrent from the Laser calculation as we do not have the H matrix
* We first calculate the Correlation matrix Tc and then calculate K
*/
void UKF::CalculateKalmanGain_Radar()
{ //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_radar_);
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

                                             //residual
    VectorXd z_diff = Zsig_radar_pred_.col(i) - z_radar_pred_;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }


  K_radar_ = Tc * S_radar_.inverse();
}


double UKF::CalculateNIS(const VectorXd &diff_meas, const MatrixXd &S) {
  //NIS = diff_meas.T * S_inv * diff_meas
  return diff_meas.transpose() * S.inverse() * diff_meas;
}

/************************************************************************
UKF Filter Implementation
End
/************************************************************************/





/************************************************************************
UKF Filter Steps
Begin
/************************************************************************/

void UKF::Init(MeasurementPackage &meas_package)
{
  time_us_ = meas_package.timestamp_;
  //TODO
  //Init x_ and P_
  H_laser_.fill(0.0);
  H_laser_(0, 0) = 1;
  H_laser_(1, 1) = 1;

  R_laser_ << std_laspx_*std_laspx_, 0,
    0, std_laspy_*std_laspy_;

  R_radar_ << std_radr_*std_radr_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0, std_radrd_*std_radrd_;

  x_.fill(0.0);
  P_.setIdentity();
  if (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER) {
    x_(0) = meas_package.raw_measurements_(0);
    x_(1) = meas_package.raw_measurements_(1);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR) {
    double phi = meas_package.raw_measurements_(1);
    x_(0) = meas_package.raw_measurements_(0)*cos(phi);
    x_(1) = meas_package.raw_measurements_(0)*sin(phi);
  }

}

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

  if (!is_initialized_) {
    is_initialized_ = true;
    cout << "Initializing filter" << endl;
    Init(meas_package);
    return;
  }

  double dt = (meas_package.timestamp_ - time_us_);
  dt /= 1000000.0;

  if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER){
    cout << "Processing LASER measurement" << endl;
    Prediction(dt);
    time_us_ = meas_package.timestamp_;
    UpdateLidar(meas_package.raw_measurements_);
  }

  else if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR) {
    cout << "Processing RADAR measurement" << endl;
    Prediction(dt);
    time_us_ = meas_package.timestamp_;
    UpdateRadar(meas_package.raw_measurements_);
  }
  
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  //Generate a set of augmented Sigma Points
  GenerateSigmaPoints();
  //cout << "Generated sigma Points" << endl;
  //cout << Xsig_aug_ << endl;
  //Predict Sigma Points
  PredictSigmaPoints(delta_t);

  //Calculate Mean and Covariance
  //Values overwrite previous state vector and covariance matrix
  CalculateMeanCovariance_State();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(VectorXd z) {
  
  //We know Lidar's H(measurement function)
  //Perform calc exactly like kalman from past project
  VectorXd y = z - H_laser_*x_;
  MatrixXd Ht = H_laser_.transpose();
  MatrixXd PHt = P_*Ht;
  MatrixXd S = H_laser_*PHt + R_laser_;
  MatrixXd K = PHt*S.inverse();

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  x_ = x_ + K*y;
  P_ = (I - K*H_laser_)*P_;

  NIS_laser_ = CalculateNIS(y, S);
  cout << "Laser NIS: " << NIS_laser_ << endl << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(VectorXd z) {
  //Radar measurements are non-linear
  //We need to perform the sigma point calculation very similar to the Prediction steps
  //We can reuse the sigma points calculated previously for this.

  //Calculate measurement and covariance
  CalculateMeanCovariance_Radar();
 
  //Residual
  VectorXd y = z - z_radar_pred_;

  CalculateKalmanGain_Radar();

  //Update the state
  MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
  x_ = x_ + K_radar_*y;
  P_ = P_ - K_radar_*S_radar_*K_radar_.transpose();

  //Normalize angle phi
  while (y(1)> M_PI) y(1) -= 2.*M_PI;
  while (y(1)<-M_PI) y(1) += 2.*M_PI;
  //Calculate NIS
  NIS_radar_ = CalculateNIS(y, S_radar_);
  cout << "Radar NIS: " << NIS_radar_ << endl << endl;
}

/************************************************************************
UKF Filter Steps
End
/************************************************************************/
