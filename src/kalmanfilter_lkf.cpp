// ------------------------------------------------------------------------------- //
// Advanced Kalman Filtering and Sensor Fusion Course - Linear Kalman Filter
//
// ####### STUDENT FILE #######
//
// Usage:
// -Rename this file to "kalmanfilter.cpp" if you want to use this code.

#include "kalmanfilter.h"
#include "utils.h"

// -------------------------------------------------- //
constexpr bool INIT_ON_FIRST_PREDICTION = true;
constexpr double INIT_POS_STD = 0.5;
constexpr double INIT_VEL_STD = 1.0;
constexpr double ACCEL_STD = 0.25;
constexpr double GPS_POS_STD = 3.0;

constexpr double INIT_POS_VAR = INIT_POS_STD*INIT_POS_STD;
constexpr double INIT_VEL_VAR = INIT_VEL_STD*INIT_VEL_STD;
constexpr double ACCEL_VAR = ACCEL_STD*ACCEL_STD;
constexpr double GPS_POS_VAR = GPS_POS_STD*GPS_POS_STD;
// -------------------------------------------------- //


void KalmanFilter::predictionStep(double dt)
{
    if (!isInitialised() && INIT_ON_FIRST_PREDICTION)
    {
        VectorXd state = Vector4d::Zero();
        MatrixXd cov = Matrix4d::Zero();

        state << 0.0, 0.0, 5*cos(M_PI/4), 5*sin(M_PI/4);
        cov(0,0) = INIT_POS_VAR;
        cov(1,1) = INIT_POS_VAR;
        cov(2,2) = INIT_VEL_VAR;
        cov(3,3) = INIT_VEL_VAR;

        setState(state);
        setCovariance(cov);
    }

    if (isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        MatrixXd F = Matrix4d::Identity();
        F(0, 2) += dt;
        F(1, 3) += dt;
        
        Vector2d w = ACCEL_VAR * Vector2d::Ones(); 
        Matrix2d Qa = w.asDiagonal();

        MatrixXd L(4, 2);
        L << 0.5*dt*dt, 0,
            0, 0.5*dt*dt,
            dt, 0,
            0, dt;

        MatrixXd R = L * Qa * L.transpose();

        state = F * state;
        cov   = F * cov * F.transpose() + R;

        setState(state);
        setCovariance(cov);

    }
}

void KalmanFilter::handleGPSMeasurement(GPSMeasurement meas)
{
    if(isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        MatrixXd H(2, 4);
        H << 1, 0, 0, 0,
            0, 1, 0, 0;

        Vector2d z{meas.x, meas.y};
        Vector2d innovation = z - H * state;

        Matrix2d R = Matrix2d::Zero();
        R(0, 0) = GPS_POS_VAR;
        R(1, 1) = GPS_POS_VAR;

        Matrix2d S = H * cov * H.transpose() + R; 
        MatrixXd K = cov * H.transpose() * S.inverse(); 

        state = state + K * innovation;
        cov = (Matrix4d::Identity() - K * H) * cov; 
        
        setState(state);
        setCovariance(cov);
    }
    else
    {
        Vector4d state{meas.x, meas.y, 0.0, 0.0};

        MatrixXd cov = Matrix4d::Zero();
        cov(0, 0) = GPS_POS_VAR;
        cov(1, 1) = GPS_POS_VAR;
        cov(2,2) = INIT_VEL_VAR;
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
        VectorXd state = getState(); // STATE VECTOR [X,Y,VX,VY]
        double psi = std::atan2(state[3],state[2]);
        double V = std::sqrt(state[2]*state[2] + state[3]*state[3]);
        return VehicleState(state[0],state[1],psi,V);
    }
    return VehicleState();
}

void KalmanFilter::predictionStep(GyroMeasurement gyro, double dt){predictionStep(dt);}
void KalmanFilter::handleLidarMeasurements(const std::vector<LidarMeasurement>& dataset, const BeaconMap& map){}
void KalmanFilter::handleLidarMeasurement(LidarMeasurement meas, const BeaconMap& map){}

