# Continuous-Discrete Multiplicative Kalman Filter (CI-MEKF) for Attitude Estimation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2023a-blue.svg)](https://www.mathworks.com/products/matlab.html)

This repository contains MATLAB models and simulation scripts for a Continuous-Discrete Multiplicative Kalman Filter (CI-MEKF), developed for robust aircraft attitude estimation.

## Project Overview

The objective of this project is to accurately estimate the orientation (attitude) of an aircraft flying a "figure-eight" trajectory using a 9-axis sensor suite (Gyroscope, Accelerometer, Magnetometer)

Estimating attitude in a global navigation frame (Fixed NED) presents significant challenges due to highly non-linear system dynamics, gyroscope bias drift, and error accumulation over time. To address this, we implemented and compared a standard Multiplicative Extended Kalman Filter (MEKF) against a Covariance-Inflated MEKF (CI-MEKF).

### System Dynamics & Sensor Fusion
The state vector is defined as $x = [\delta\tilde{q}, b_g]^T$, encompassing the quaternion error and gyroscope bias.
*   **Prediction Step:** High-frequency gyroscope data is used to propagate the system dynamics and extrapolate the error covariance.
*   **Correction Step:** Accelerometer and magnetometer measurements are fused to calculate the Kalman Gain and correct the apriori state estimate.

## Simulation & Monte Carlo Analysis
To validate the robustness of the CI-MEKF, we modeled severe magnetic disturbances (a 100 mT sudden jump) in the navigation frame to simulate real-world sensor anomalies. 

A comprehensive Monte Carlo Simulation was performed across 6 different scenarios, with 50 trials per scenario. The evaluation metric was the Root Mean Square Error (RMSE) for Roll, Pitch, Yaw, and Gyroscope Bias

### Key Findings
1.  **Disturbance Rejection:** During sudden magnetic field jumps, the standard MEKF estimates diverged significantly (e.g., ~40-degree errors in Yaw). 
2.  **Covariance Inflation Advantage:** The CI-MEKF successfully detected the anomalous measurements, dynamically reduced the Kalman gain, and relied on the system dynamics, capping the maximum error at ~15 degrees in Roll.
3.  **Overall Performance:** Across all 6 Monte Carlo scenarios, the CI-MEKF demonstrated consistently lower angular RMSE compared to the standard MEKF despite injected noise and measurement faults.

## Execution
To run the simulation and generate the comparison plots:
1. Open MATLAB.
2. Navigate to this directory.
3. Run the main script:
```matlab
run('CIMEKFSimulation.m') % Update this filename to match your actual script
```

## For More Details & Plots
- Refer to the presentation file for the project, it's prepared in Turkish, but you can easily have it translated.