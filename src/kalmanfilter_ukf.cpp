// ------------------------------------------------------------------------------- //
// Advanced Kalman Filtering and Sensor Fusion Course - Unscented Kalman Filter
//
// ####### STUDENT FILE #######
//
// Usage:
// -Rename this file to "kalmanfilter.cpp" if you want to use this code.

#include "kalmanfilter.h"
#include "utils.h"

// -------------------------------------------------- //
// YOU CAN USE AND MODIFY THESE CONSTANTS HERE
constexpr double ACCEL_STD = 1.0;
constexpr double GYRO_STD = 0.01/180.0 * M_PI;
constexpr double INIT_VEL_STD = 10.0;
constexpr double INIT_PSI_STD = 45.0/180.0 * M_PI;
constexpr double GPS_POS_STD = 3.0;
constexpr double LIDAR_RANGE_STD = 3.0;
constexpr double LIDAR_THETA_STD = 0.02;
// -------------------------------------------------- //

VectorXd normaliseState(VectorXd state)
{
    state(2) = wrapAngle(state(2));
    return state;
}
VectorXd normaliseLidarMeasurement(VectorXd meas)
{
    meas(1) = wrapAngle(meas(1));
    return meas;
}
std::vector<VectorXd> generateSigmaPoints(VectorXd state, MatrixXd cov)
{
    std::vector<VectorXd> sigmaPoints;

    int n = state.rows();
    int kappa = 3 - n;
    Eigen::MatrixXd L( cov.llt().matrixL() );
    L = sqrt((kappa+n))*L;

    sigmaPoints.push_back(state);

    for(int i=0; i<n; i++)
    {   
        VectorXd del_x(n);
        del_x = L.col(i);
        sigmaPoints.push_back(state + del_x);
        sigmaPoints.push_back(state - del_x);
    }

    return sigmaPoints;
}

std::vector<double> generateSigmaWeights(unsigned int numStates)
{
    std::vector<double> weights;

    int kappa = 3 - numStates;
    int denom = kappa + numStates;

    weights.clear();
    weights.push_back(static_cast<double>(kappa) / denom);

    for (int i = 0; i < 2 * numStates; ++i) {
        weights.push_back(1.0 / (2 * denom));
    }

    return weights;
}

VectorXd lidarMeasurementModel(VectorXd aug_state, double beaconX, double beaconY)
{
    VectorXd z_hat = VectorXd::Zero(2);

    double x     = aug_state(0);
    double y     = aug_state(1);
    double psi   = aug_state(2);
    double v     = aug_state(3);
    double wr   = aug_state(4);
    double wt   = aug_state(5);

    double range = sqrt(pow((beaconX - x), 2) + pow((beaconY - y), 2)) + wr;
    double theta = atan2((beaconY - y), (beaconX - x)) - psi + wt;
    
    z_hat << range, theta;

    return z_hat;
}

VectorXd vehicleProcessModel(VectorXd aug_state, double psi_dot, double dt)
{
    VectorXd updated_aug_state = aug_state;

    double x     = aug_state(0);
    double y     = aug_state(1);
    double psi   = aug_state(2);
    double v     = aug_state(3);
    double w     = aug_state(4);
    double a     = aug_state(5);

    updated_aug_state(0) = x + dt * v * std::cos(psi);
    updated_aug_state(1) = y + dt * v * std::sin(psi);
    updated_aug_state(2) = psi + dt * (psi_dot + w);
    updated_aug_state(3) = v + dt * a;

    return updated_aug_state;
}

void KalmanFilter::handleLidarMeasurement(LidarMeasurement meas, const BeaconMap& map)
{
    if (isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        int ns = state.rows() + 2;

        BeaconData map_beacon = map.getBeaconWithId(meas.id); // Match Beacon with built in Data Association Id
        if (meas.id != -1 && map_beacon.id != -1) // Check that we have a valid beacon match
        {
            VectorXd augmented_state(6);
            augmented_state << state(0), state(1), state(2), state(3), 0, 0;

            MatrixXd augmented_cov = MatrixXd::Zero(ns, ns);
            augmented_cov.topLeftCorner(4, 4) = cov;
            augmented_cov(4, 4) = LIDAR_RANGE_STD * LIDAR_RANGE_STD;
            augmented_cov(5, 5) = LIDAR_THETA_STD * LIDAR_THETA_STD;

            std::vector<VectorXd> sigma_points = generateSigmaPoints(augmented_state, augmented_cov);
            std::vector<double> weights = generateSigmaWeights(ns);

            std::vector<Vector2d> sigma_points_projected;
            for (const auto& sigma : sigma_points) {
                sigma_points_projected.push_back(lidarMeasurementModel(sigma, map_beacon.x, map_beacon.y));
            }

            Vector2d z; 
            z << meas.range, meas.theta;

            Vector2d z_hat = Vector2d::Zero();
            VectorXd x_mean = VectorXd::Zero(6);
            for (size_t i = 0; i < sigma_points.size(); ++i) {
                z_hat += weights[i] * sigma_points_projected[i];
                x_mean += weights[i] * sigma_points[i];
            }

            Matrix2d S = Matrix2d::Zero();
            for (size_t i = 0; i < sigma_points.size(); ++i) {
                Vector2d diff = normaliseLidarMeasurement(sigma_points_projected[i] - z_hat);
                S += weights[i] * diff * diff.transpose();
            }

            MatrixXd Pxz = MatrixXd::Zero(6, 2);
            for (size_t i = 0; i < sigma_points.size(); ++i) {
                VectorXd diff_x = normaliseState(sigma_points[i] - x_mean);
                Vector2d diff_z = normaliseLidarMeasurement(sigma_points_projected[i] - z_hat);
                Pxz += weights[i] * diff_x * diff_z.transpose();
            }

            MatrixXd K = Pxz * S.inverse();

            Vector2d z_diff = normaliseLidarMeasurement(z - z_hat);
            VectorXd update = K * z_diff;

            state += update.head(4); 
            augmented_cov -= K * S * K.transpose();

            state = normaliseState(state);
            cov = augmented_cov.topLeftCorner(4, 4);
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
        
        int nx = state.rows();
        int ns = nx + 2;

        VectorXd aug_state = VectorXd::Zero(ns);
        aug_state.head(nx) = state;

        MatrixXd aug_cov = MatrixXd::Zero(ns, ns);
        aug_cov.topLeftCorner(nx, nx) = cov;
        aug_cov(nx, nx) = GYRO_STD * GYRO_STD;
        aug_cov(nx + 1, nx + 1) = ACCEL_STD * ACCEL_STD;

        std::vector<VectorXd> sigma_points = generateSigmaPoints(aug_state, aug_cov);
        std::vector<double> weights = generateSigmaWeights(ns);
        assert(weights.size() == sigma_points.size());

        std::vector<VectorXd> sigma_points_processed;
        sigma_points_processed.reserve(sigma_points.size());
        for (const auto& pt : sigma_points) {
            sigma_points_processed.push_back(vehicleProcessModel(pt, gyro.psi_dot, dt));
        }

        VectorXd y = VectorXd::Zero(ns);
        for (size_t i = 0; i < sigma_points_processed.size(); ++i) {
            y += weights[i] * sigma_points_processed[i];
        }

        MatrixXd new_cov = MatrixXd::Zero(ns, ns);
        for (size_t i = 0; i < sigma_points_processed.size(); ++i) {
            VectorXd diff = normaliseState(sigma_points_processed[i] - y);
            new_cov += weights[i] * diff * diff.transpose();
        }

        state = normaliseState(y.head(nx));
        cov = new_cov.topLeftCorner(nx, nx);

        setState(state);
        setCovariance(cov);
    } 
}

void KalmanFilter::handleGPSMeasurement(GPSMeasurement meas)
{
    // All this code is the same as the LKF as the measurement model is linear
    // so the UKF update state would just produce the same result.
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
        cov = (MatrixXd::Identity(4,4) - K*H) * cov;

        setState(state);
        setCovariance(cov);
    }
    else
    {
        VectorXd state = Vector4d::Zero();
        MatrixXd cov = Matrix4d::Zero();

        state(0) = meas.x;
        state(1) = meas.y;
        cov(0,0) = GPS_POS_STD*GPS_POS_STD;
        cov(1,1) = GPS_POS_STD*GPS_POS_STD;
        cov(2,2) = INIT_PSI_STD*INIT_PSI_STD;
        cov(3,3) = INIT_VEL_STD*INIT_VEL_STD;

        setState(state);
        setCovariance(cov);
    }             
}

void KalmanFilter::handleLidarMeasurements(const std::vector<LidarMeasurement>& dataset, const BeaconMap& map)
{
    // Assume No Correlation between the Measurements and Update Sequentially
    for(const auto& meas : dataset) {handleLidarMeasurement(meas, map);}
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
