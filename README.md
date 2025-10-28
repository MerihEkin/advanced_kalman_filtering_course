# Advanced Kalman Filtering and Sensor Fusion Simulation

This repository contains my solutions for the **Advanced Kalman Filtering and Sensor Fusion** course on Udemy by Dr. Steven Dumble. It includes implementations of:

- **Linear Kalman Filter (LKF):** Implemented in `kalmanfilter_lkf.cpp`  
- **Extended Kalman Filter (EKF):** Implemented in `kalmanfilter_ekf.cpp`  
- **Unscented Kalman Filter (UKF):** Implemented in `kalmanfilter_ukf.cpp`  
- **Real-world UKF Variant:** Implemented in `kalmanfilter.cpp`. This version initializes with GPS and LiDAR data, and includes handling for faulty GPS measurements, gyroscope bias, and missing LiDAR beacon IDs.


## Setup

This project uses the **Ubuntu 20.04.2 LTS 64-bit VM C++ development environment** provided in the course. Refer to the *Setting up the Development Environment* lesson for full instructions.

### 0. Install Dependencies

This project will use the Ubuntu 64 20.04.2.0 LTS VM C++ development environment that is setup for this course. (Refer to the Setting up the Development Environment Lesson) Please follow the steps below to compile the simulation.

 0. Install the dependencies
 ```
 sudo apt install libeigen3-dev libsdl2-dev libsdl2-ttf-dev
 ```
 
 1. Clone the repository
 ```
 git clone https://github.com/technitute/AKFSF-Simulation-CPP.git
 ```
 2. Setup the cmake build
 ```
 cd AKFSF-Simulation-CPP
 mkdir build
 cd build
 cmake ../
 ```

 3. Compile the code
 ```
 make
 ```
 
 3. You should now be able to and run the estimation simulator
 ```
 ./AKFSF-Simulation
 ```

## Simulation Operation ##
The simulation can be run with different motion and sensor profiles to test the different scenarios and evaluate how the filter performs with different conditions. The different profiles can be activated pressing the number keys 1-9,0, to active the corresponding profile.

### Motion Profiles ###
* 1 - Constant Velocity + GPS + GYRO + Zero Initial Conditions
* 2 - Constant Velocity + GPS + GYRO + Non-zero Initial Conditions
* 3 - Constant Speed Profile + GPS + GYRO
* 4 - Variable Speed Profile + GPS + GYRO
* 5 - Constant Velocity + GPS + GYRO + LIDAR+ Zero Initial Conditions
* 6 - Constant Velocity + GPS + GYRO + LIDAR + Non-zero Initial Conditions
* 7 - Constant Speed Profile + GPS + GYRO + LIDAR
* 8 - Variable Speed Profile + GPS + GYRO + LIDAR
* 9 - CAPSTONE
* 0 - CAPSTONE BONUS (with No Lidar Data Association)



## The Tasks ##
The tasks for this simulator can be broken down into 4 different areas:

- Linear Kalman Filter (LKF) – Profiles 1–4
- Extended Kalman Filter (EKF) – Profiles 1–8
- Unscented Kalman Filter (UKF) – Profiles 1–8
- Capstone Project – All profiles


## Authors ##

This project was developed for the Technitute Course - Advanced Kalman Filtering and Sensor Fusion. Developed and produced by Dr. Steven Dumble.