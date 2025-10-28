#ifndef INCLUDE_AKFSFSIM_KALMANFILTER_H
#define INCLUDE_AKFSFSIM_KALMANFILTER_H

#include <vector>
#include <Eigen/Dense>

#include "car.h"
#include "sensors.h"
#include "beacons.h"

using Eigen::VectorXd;
using Eigen::Vector2d;
using Eigen::Vector4d;

using Eigen::MatrixXd;
using Eigen::Matrix2d;
using Eigen::Matrix4d;

class KalmanFilterBase
{
    public:

        KalmanFilterBase() : m_initialised(false) {}
        virtual ~KalmanFilterBase(){}
        void reset(){m_initialised = false;}
        bool isInitialised() const {return m_initialised;}

    protected:
    
        VectorXd getState() const {return m_state;}
        MatrixXd getCovariance()const {return m_covariance;}
        void setState(const VectorXd& state ) {m_state = state; m_initialised = true;}
        void setCovariance(const MatrixXd& cov ){m_covariance = cov;}

        std::vector<double> getInitX() const { return init_x;}
        void setInitX(const std::vector<double>& initX) { init_x = initX;}
        std::vector<double> getInitY() const { return init_y;}
        void setInitY(const std::vector<double>& initY) { init_y = initY;}
        std::vector<double> getInitPsi() const { return init_psi;}
        void setInitPsi(const std::vector<double>& initPsi) { init_psi = initPsi;}
        std::vector<double> getInitPsiDot() const { return init_psi_dot;}
        void setInitPsiDot(const std::vector<double>& initPsiDot) { init_psi_dot = initPsiDot;}

    private:
        bool m_initialised;
        VectorXd m_state;
        MatrixXd m_covariance;

        std::vector<double> init_x;
        std::vector<double> init_y;
        std::vector<double> init_psi;
        std::vector<double> init_psi_dot;
};

class KalmanFilter : public KalmanFilterBase
{
    public:

        VehicleState getVehicleState();
        Matrix2d getVehicleStatePositionCovariance();

        void predictionStep(double dt);
        void predictionStep(GyroMeasurement gyro, double dt);
        void handleLidarMeasurements(const std::vector<LidarMeasurement>& meas, const BeaconMap& map);
        void handleLidarMeasurement(LidarMeasurement meas, const BeaconMap& map);
        void handleGPSMeasurement(GPSMeasurement meas);

};

#endif  // INCLUDE_AKFSFSIM_KALMANFILTER_H