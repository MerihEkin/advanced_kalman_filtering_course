// ------------------------------------------------------------------------------- //
// Advanced Kalman Filtering and Sensor Fusion Course - Extended Kalman Filter
//
// ####### STUDENT FILE #######
//
// Usage:
// -Rename this file to "kalmanfilter.cpp" if you want to use this code.

#include "kalmanfilter.h"
#include "utils.h"

// -------------------------------------------------- //
constexpr double ACCEL_STD = 1.0;
constexpr double GYRO_STD = 0.01/180.0 * M_PI;
constexpr double INIT_VEL_STD = 10.0;
constexpr double INIT_PSI_STD = 45.0/180.0 * M_PI;
constexpr double GPS_POS_STD = 3.0;
constexpr double LIDAR_RANGE_STD = 3.0;
constexpr double LIDAR_THETA_STD = 0.02;

constexpr double ACCEL_VAR = ACCEL_STD * ACCEL_STD;
constexpr double GYRO_VAR = GYRO_STD * GYRO_STD;
constexpr double INIT_VEL_VAR = INIT_VEL_STD * INIT_VEL_STD;
constexpr double INIT_PSI_VAR = INIT_PSI_STD * INIT_PSI_STD;
constexpr double GPS_POS_VAR = GPS_POS_STD * GPS_POS_STD;
constexpr double LIDAR_RANGE_VAR = LIDAR_RANGE_STD * LIDAR_RANGE_STD;
constexpr double LIDAR_THETA_VAR = LIDAR_THETA_STD * LIDAR_THETA_STD;
// -------------------------------------------------- //

const Matrix4d Qa = (Matrix4d() << 
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, GYRO_VAR, 0.0,
    0.0, 0.0, 0.0, ACCEL_VAR).finished();

const Matrix2d R = (Matrix2d() <<
    LIDAR_RANGE_VAR, 0.0, 0.0, LIDAR_THETA_VAR).finished();

void KalmanFilter::handleLidarMeasurements(const std::vector<LidarMeasurement>& dataset, const BeaconMap& map)
{
    // Assume No Correlation between the Measurements and Update Sequentially
    for(const auto& meas : dataset) {handleLidarMeasurement(meas, map);}
}

void KalmanFilter::handleLidarMeasurement(LidarMeasurement meas, const BeaconMap& map)
{
    if (isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        BeaconData map_beacon = map.getBeaconWithId(meas.id);
        if (meas.id != -1 && map_beacon.id != -1)
        {   
            double x = state(0);
            double y = state(1);
            double psi = state(2);

            double Lx = map_beacon.x;
            double Ly = map_beacon.y;

            double dx = Lx - x;
            double dy = Ly - y;

            double meas_range = meas.range;
            double meas_theta = meas.theta;

            double range_predict = sqrt(dx * dx + dy * dy);
            double theta_predict = atan2(dy, dx) - psi;
            theta_predict = wrapAngle(theta_predict);

            MatrixXd H(2, 4);
            double d = range_predict;

            H << -dx / d,         -dy / d,         0.0, 0.0,
                dy / (d * d),   -dx / (d * d),  -1.0, 0.0;

            Vector2d nu{
                meas_range - range_predict,
                wrapAngle(meas_theta - theta_predict)
            };

            MatrixXd S = H * cov * H.transpose() + R;
            MatrixXd K = cov * H.transpose() * S.inverse();

            state = state + K * nu;
            cov = (Matrix4d::Identity() - K * H) * cov;
        }

        setState(state);
        setCovariance(cov);
    }
}

void KalmanFilter::predictionStep(GyroMeasurement gyro, double dt)
{
    if (isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        // Calculate jacobian around the linearization point x
        Matrix4d F = Matrix4d::Identity();

        // state is in the form [PX, PY, PSI, V]
        double v   = state(3);
        double psi = state(2);

        F(0, 2) -= dt * v * std::sin(psi);
        F(1, 2) += dt * v * std::cos(psi);
        F(0, 3) += dt * std::cos(psi);
        F(1, 3) += dt * std::sin(psi);
        
        
        // State prediction:
        Eigen::Vector4d u;
        u << v * std::cos(psi),
            v * std::sin(psi),
            gyro.psi_dot,
            0.0;
        state += dt * u;
        state(2) = wrapAngle(state(2));

        // Covariance update
        cov = F * (cov * F.transpose()) + (dt*dt) * Qa;

        setState(state);
        setCovariance(cov);
    } 
}

void KalmanFilter::handleGPSMeasurement(GPSMeasurement meas)
{
    // All this code is the same as the LKF as the measurement model is linear
    // so the EKF update state would just produce the same result.
    if(isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        VectorXd z = Vector2d::Zero();
        MatrixXd H = MatrixXd(2,4);
        MatrixXd R = Matrix2d::Zero();

        z << meas.x,meas.y;
        H << 1,0,0,0,0,1,0,0;
        R(0,0) = GPS_POS_STD*GPS_POS_STD;
        R(1,1) = GPS_POS_STD*GPS_POS_STD;

        VectorXd z_hat = H * state;
        VectorXd y = z - z_hat;
        MatrixXd S = H * cov * H.transpose() + R;
        MatrixXd K = cov*H.transpose()*S.inverse();

        state = state + K*y;
        cov = (Matrix4d::Identity() - K*H) * cov;

        setState(state);
        setCovariance(cov);
    }
    else
    {
        VectorXd state = Vector4d::Zero();
        MatrixXd cov = Matrix4d::Zero();

        state(0) = meas.x;
        state(1) = meas.y;
        cov(0,0) = GPS_POS_VAR;
        cov(1,1) = GPS_POS_VAR;
        cov(2,2) = INIT_PSI_VAR;
        cov(3,3) = INIT_VEL_VAR;

        setState(state);
        setCovariance(cov);
    }            
}

Matrix2d KalmanFilter::getVehicleStatePositionCovariance()
{
    Matrix2d pos_cov = Matrix2d::Zero();
    MatrixXd cov = getCovariance();
    if (isInitialised() && cov.size() != 0){pos_cov << cov(0,0), cov(0,1), cov(1,0), cov(1,1);}
    return pos_cov;
}

VehicleState KalmanFilter::getVehicleState()
{
    if (isInitialised())
    {
        VectorXd state = getState(); // STATE VECTOR [X,Y,PSI,V,...]
        return VehicleState(state[0],state[1],state[2],state[3]);
    }
    return VehicleState();
}

void KalmanFilter::predictionStep(double dt){}