/*
// ------------------------------------------------------------------------------- //
// Advanced Kalman Filtering and Sensor Fusion Course - Unscented Kalman Filter
//
// ####### STUDENT FILE #######
//
// Usage:
// - Rename this file to "kalmanfilter.cpp" if you want to use this code.
//
// 1. Compensate for constant sensor bias in gyroscope measurements (done)
//    - Tip: augment the state vector to include gyroscope bias.
//    - Random or stochastic biases cannot be corrected through calibration alone.
//      Instead, we use bias parameter estimation. The bias evolves as:
//         b = b + w_b, where db/dt = 0.
//    - The bias must be observable by the system to be estimated correctly.
//
// 2. Implement rejection of faulty GPS measurements (done)
//    - Use thresholding based on the innovation metric: e * e^T < threshold.
//    - Apply Fault Detection and Isolation (FDI) using methods such as:
//         * Innovation checking
//         * Normalized innovation squared (NIS)
//         * Chi-squared test (Pearson's chi-squared test)
//    - For 2D measurements, the chi-squared threshold is 5.99 at a 0.05 significance level (p = 0.05).
//
// 3. Handle missing LiDAR beacon IDs (done)
//    - Estimate and associate beacon IDs before processing the measurements.
//    - Filter out landmarks which can be in view
//    - Use current estimated position and heading to best match up possible landmarks to find data association.
//
// 4. Improve state initialization to include yaw angle (ψ) and initial velocity estimates (done)
//    - The filter should not start until a reliable estimate of the initial state and its uncertainty is available.
//    - If reasonable assumptions cannot be made, use an initialization filter or process to estimate the starting state.
//    - Use the GPS position measurements to initialize the X/Y position states
//    - Filter GPS measurements to estimate X/Y velocity or total velocity
//    - GPS Position and LIDAR measurements to landmarks to estimate heading of the vehicle.
//    - All the error covariances should be as small as possible to avoid non-linear effects.
*/

#include "kalmanfilter.h"
#include "utils.h"
#include <tuple>

constexpr double ACCEL_STD = 1.0;
constexpr double GYRO_STD = 0.01 / 180.0 * M_PI;
constexpr double INIT_VEL_STD = 10.0;
constexpr double INIT_PSI_STD = (45.0 / 180.0) * M_PI;
constexpr double GPS_POS_STD = 5.0;
constexpr double LIDAR_RANGE_STD = 3.0;
constexpr double LIDAR_THETA_STD = 0.02;

constexpr double CHI_SQUARED_TEST_THRES = 5.99;
constexpr double STATE_INIT_THRES = 10.0;
constexpr double LIDAR_RANGE = 180.0;
constexpr double INIT_BIAS_STD = 10.0 / 180.0 * M_PI;
constexpr double GATING_THRESHOLD = 9.21;

double calculate_mean(const std::vector<double>& xs) {
    int n = xs.size();
    if (n == 0) return 0.0;
    double sum = 0.0;
    for (double x : xs) sum += x;
    return sum / n;
}

VectorXd normaliseState(VectorXd state) {
    state(2) = wrapAngle(state(2));
    return state;
}

VectorXd normaliseLidarMeasurement(VectorXd meas) {
    meas(1) = wrapAngle(meas(1));
    return meas;
}

std::vector<VectorXd> generateSigmaPoints(VectorXd state, MatrixXd cov) {
    std::vector<VectorXd> sigmaPoints;
    int n = state.rows();
    int kappa = 3 - n;
    Eigen::MatrixXd L(cov.llt().matrixL());
    L = sqrt((kappa + n)) * L;
    sigmaPoints.push_back(state);
    for (int i = 0; i < n; i++) {
        VectorXd del_x = L.col(i);
        sigmaPoints.push_back(state + del_x);
        sigmaPoints.push_back(state - del_x);
    }
    return sigmaPoints;
}

std::vector<double> generateSigmaWeights(unsigned int numStates) {
    std::vector<double> weights;
    int kappa = 3 - numStates;
    int denom = kappa + numStates;
    weights.push_back(static_cast<double>(kappa) / denom);
    for (int i = 0; i < 2 * numStates; ++i) {
        weights.push_back(1.0 / (2 * denom));
    }
    return weights;
}

VectorXd lidarMeasurementModel(VectorXd aug_state, double beaconX, double beaconY) {
    VectorXd z_hat = VectorXd::Zero(2);
    double x = aug_state(0);
    double y = aug_state(1);
    double psi = aug_state(2);
    double v = aug_state(3);
    double bias = aug_state(4);
    double wr = aug_state(5);
    double wt = aug_state(6);
    double range = sqrt(pow((beaconX - x), 2) + pow((beaconY - y), 2)) + wr;
    double theta = atan2((beaconY - y), (beaconX - x)) - psi + wt;
    z_hat << range, theta;
    return z_hat;
}

VectorXd vehicleProcessModel(VectorXd aug_state, double psi_dot, double dt) {
    VectorXd updated_aug_state = aug_state;
    double x = aug_state(0);
    double y = aug_state(1);
    double psi = aug_state(2);
    double v = aug_state(3);
    double bias = aug_state(4);
    double w = aug_state(5);
    double a = aug_state(6);
    updated_aug_state(0) = x + dt * v * std::cos(psi);
    updated_aug_state(1) = y + dt * v * std::sin(psi);
    updated_aug_state(2) = psi + dt * (psi_dot - bias + w);
    updated_aug_state(3) = v + dt * a;
    return updated_aug_state;
}

void KalmanFilter::handleLidarMeasurement(LidarMeasurement meas, const BeaconMap& map) {
    if (isInitialised()) {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();
        int nx = state.rows();
        int ns = nx + 2;

        BeaconData map_beacon = map.getBeaconWithId(meas.id);
        if (map_beacon.id == -1) {
            std::vector<BeaconData> beacons = map.getBeaconsWithinRange(state(0), state(1), LIDAR_RANGE);
            if (!beacons.empty()) {
                Matrix2d R;
                R << LIDAR_RANGE_STD * LIDAR_RANGE_STD, 0,
                     0, LIDAR_THETA_STD * LIDAR_THETA_STD;

                double x = state(0), y = state(1), psi = state(2);
                Vector2d z_meas; z_meas << meas.range, meas.theta;

                double min_cost = std::numeric_limits<double>::max();
                int best_id = -1;

                for (const BeaconData& beacon : beacons) {
                    double dx = beacon.x - x;
                    double dy = beacon.y - y;
                    double predicted_range = std::sqrt(dx * dx + dy * dy);
                    double predicted_theta = std::atan2(dy, dx) - psi;
                    predicted_theta = std::atan2(std::sin(predicted_theta), std::cos(predicted_theta));

                    Vector2d z_pred; z_pred << predicted_range, predicted_theta;
                    Vector2d error = z_meas - z_pred;
                    error(1) = std::atan2(std::sin(error(1)), std::cos(error(1)));

                    double cost = error.transpose() * R.inverse() * error;
                    if (cost < min_cost && cost < GATING_THRESHOLD) {
                        min_cost = cost;
                        best_id = beacon.id;
                        map_beacon = beacon;
                    }
                }

                map_beacon.id = best_id;
                meas.id = best_id;
            }
        }

        if (meas.id != -1 && map_beacon.id != -1) {
            VectorXd augmented_state = VectorXd::Zero(ns);
            augmented_state.head(nx) = state;

            MatrixXd augmented_cov = MatrixXd::Zero(ns, ns);
            augmented_cov.topLeftCorner(nx, nx) = cov;
            augmented_cov(nx, nx) = LIDAR_RANGE_STD * LIDAR_RANGE_STD;
            augmented_cov(nx + 1, nx + 1) = LIDAR_THETA_STD * LIDAR_THETA_STD;

            std::vector<VectorXd> sigma_points = generateSigmaPoints(augmented_state, augmented_cov);
            std::vector<double> weights = generateSigmaWeights(ns);

            std::vector<Vector2d> sigma_points_projected;
            for (const auto& sigma : sigma_points) {
                sigma_points_projected.push_back(lidarMeasurementModel(sigma, map_beacon.x, map_beacon.y));
            }

            Vector2d z; z << meas.range, meas.theta;
            Vector2d z_hat = Vector2d::Zero();
            VectorXd x_mean = VectorXd::Zero(ns);

            for (size_t i = 0; i < sigma_points.size(); ++i) {
                z_hat += weights[i] * sigma_points_projected[i];
                x_mean += weights[i] * sigma_points[i];
            }

            Matrix2d S = Matrix2d::Zero();
            for (size_t i = 0; i < sigma_points.size(); ++i) {
                Vector2d diff = normaliseLidarMeasurement(sigma_points_projected[i] - z_hat);
                S += weights[i] * diff * diff.transpose();
            }

            MatrixXd Pxz = MatrixXd::Zero(ns, 2);
            for (size_t i = 0; i < sigma_points.size(); ++i) {
                VectorXd diff_x = normaliseState(sigma_points[i] - x_mean);
                Vector2d diff_z = normaliseLidarMeasurement(sigma_points_projected[i] - z_hat);
                Pxz += weights[i] * diff_x * diff_z.transpose();
            }

            MatrixXd K = Pxz * S.inverse();
            Vector2d z_diff = normaliseLidarMeasurement(z - z_hat);
            VectorXd update = K * z_diff;

            state += update.head(nx);
            augmented_cov -= K * S * K.transpose();
            state = normaliseState(state);
            cov = augmented_cov.topLeftCorner(nx, nx);
        }

        setState(state);
        setCovariance(cov);
    } else {
        std::vector<double> init_x = getInitX();
        std::vector<double> init_y = getInitY();
        std::vector<double> init_psi = getInitPsi();
        int n = init_x.size();

        if (n != 0) {
            double x = init_x[n - 1];
            double y = init_y[n - 1];
            BeaconData map_beacon = map.getBeaconWithId(meas.id);

            if (map_beacon.id == -1) {
                std::vector<BeaconData> beacons = map.getBeaconsWithinRange(x, y, LIDAR_RANGE);
                if (!beacons.empty()) {
                    double min_cost = std::numeric_limits<double>::max();
                    int best_id = -1;

                    for (const BeaconData& beacon : beacons) {
                        double dx = beacon.x - x;
                        double dy = beacon.y - y;
                        double predicted_range = std::sqrt(dx * dx + dy * dy);
                        double cost = abs(predicted_range - meas.range);

                        if (cost < min_cost) {
                            min_cost = cost;
                            best_id = beacon.id;
                            map_beacon = beacon;
                        }
                    }

                    map_beacon.id = best_id;
                    meas.id = best_id;
                }
            }

            if (meas.id != -1 && map_beacon.id != -1) {
                double psi = wrapAngle(atan2(map_beacon.y - y, map_beacon.x - x) - meas.theta);
                init_psi.push_back(psi);
                setInitPsi(init_psi);
            }
        }
    }
}

void KalmanFilter::predictionStep(GyroMeasurement gyro, double dt) {
    if (isInitialised()) {
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
    } else {
        std::vector<double> init_psi_dot = getInitPsiDot();
        init_psi_dot.push_back(gyro.psi_dot);
        setInitPsiDot(init_psi_dot);
    }
}

void KalmanFilter::handleGPSMeasurement(GPSMeasurement meas) {
    if (isInitialised()) {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();
        int nx = state.rows();

        VectorXd z(2);
        z << meas.x, meas.y;

        MatrixXd H = MatrixXd::Zero(2, nx);
        H(0, 0) = 1.0;
        H(1, 1) = 1.0;

        Matrix2d R;
        R << GPS_POS_STD * GPS_POS_STD, 0,
             0, GPS_POS_STD * GPS_POS_STD;

        VectorXd z_hat = H * state;
        VectorXd y = z - z_hat;
        MatrixXd S = H * cov * H.transpose() + R;
        MatrixXd K = cov * H.transpose() * S.inverse();

        double d_squared = y.transpose() * S.inverse() * y;

        if (d_squared < CHI_SQUARED_TEST_THRES) {
            state = state + K * y;
            cov = (MatrixXd::Identity(nx, nx) - K * H) * cov;

            setState(state);
            setCovariance(cov);
        }
    } else {
        std::vector<double> init_x = getInitX();
        std::vector<double> init_y = getInitY();

        init_x.push_back(meas.x);
        init_y.push_back(meas.y);
        setInitX(init_x);
        setInitY(init_y);

        if (init_x.size() == 3) {
            VectorXd state = VectorXd::Zero(5);
            MatrixXd cov = MatrixXd::Zero(5, 5);

            std::vector<double> init_psi = getInitPsi();
            std::vector<double> init_psi_dot = getInitPsiDot();

            state(0) = calculate_mean(init_x);
            state(1) = calculate_mean(init_y);
            state(2) = calculate_mean(init_psi);
            state(3) = 0.0;
            state(4) = calculate_mean(init_psi_dot);

            cov(0, 0) = GPS_POS_STD * GPS_POS_STD;
            cov(1, 1) = GPS_POS_STD * GPS_POS_STD;
            cov(2, 2) = INIT_PSI_STD * INIT_PSI_STD;
            cov(3, 3) = INIT_VEL_STD * INIT_VEL_STD;
            cov(4, 4) = INIT_BIAS_STD * INIT_BIAS_STD;

            init_x.clear();
            init_y.clear();
            init_psi.clear();
            init_psi_dot.clear();

            setInitX(init_x);
            setInitY(init_y);
            setInitPsi(init_psi);
            setInitPsiDot(init_psi_dot);

            setState(state);
            setCovariance(cov);
        }
    }
}

void KalmanFilter::handleLidarMeasurements(const std::vector<LidarMeasurement>& dataset, const BeaconMap& map) {
    for (const auto& meas : dataset) handleLidarMeasurement(meas, map);
}

Matrix2d KalmanFilter::getVehicleStatePositionCovariance() {
    Matrix2d pos_cov = Matrix2d::Zero();
    MatrixXd cov = getCovariance();
    if (isInitialised() && cov.size() != 0) {
        pos_cov << cov(0, 0), cov(0, 1), cov(1, 0), cov(1, 1);
    }
    return pos_cov;
}

VehicleState KalmanFilter::getVehicleState() {
    if (isInitialised()) {
        VectorXd state = getState();
        return VehicleState(state[0], state[1], state[2], state[3]);
    }
    return VehicleState();
}

void KalmanFilter::predictionStep(double dt) {}
