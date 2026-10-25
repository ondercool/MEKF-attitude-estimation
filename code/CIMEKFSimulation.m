%%
% True Trajectory
t_end = 60;
dt = 0.025;
N = t_end / dt;
t = (0:dt:t_end);
R = 60;
H = -30;
dH = 10;
[true_pos, true_quat, true_vel, true_acc, true_angvel] = figure8Trajectory(t_end, dt, R, H, dH);

%%
% generate magnetic anomaly(disturbance)
% magnetic field disturbance will be taken about 100 u Tesla
d_N = generateMagneticDisturbance(t, round(t_end * [0.2, 0.4, 0.6, 0.8], 0), [2, 1.5, 3, 2.5] / dt , 0.4 * [0.8 1.2 1.6]);
% define bias and init angle estimates for each case
biasInitErrs = [
    0 0 0; 
    0 0 0; 
    0 0 0; 
    1.0 -1.25 1.5;
    2.0 -2.5 3.0;
    4.0 -5.0 6.0]';
angleInitErrs = [
    0 0 0; 
    10 10 10;
    20 20 20;
    10 10 10;
    20 20 20; 
    60 60 60]';
%%
% Define constants
inc_angle = deg2rad(60.0);
dec_angle = deg2rad(5.0);
% Magnetic Field is about 50 u tesla->convert to NED measurements
mag_vec = 0.5 * [cos(inc_angle)*cos(dec_angle) cos(inc_angle)*sin(dec_angle) sin(inc_angle)]';
g = 9.8065;
var_g = 0.003;
var_b = 1e-5; % d / dt (b_g) = w_b, w_b is a Gaussian Noise in system error dynamics
var_a = 0.0056;
var_m = 0.001;
var_dq_init = 0.25; %For initial error covariance matrix
var_b_init = 0.01; %For initial error covariance matrix
P0 = [var_dq_init * eye(3) zeros(3,3);
     zeros(3,3) var_b_init * eye(3)];
Q = Variances2QMat(var_g, var_b, dt);
R_m = eye(3) * var_m;
R_a = eye(3) * var_a;

%%
% Define true sensor measurements in the case of no disturbance
z_acc0 = zeros(3, N-1); % Acceleration in the body frame
z_mag0 = zeros(3, N-1); % Magnetic field vector in the body frame
z_gyro0 = zeros(3, N-1);
g_vec = [0 0 9.8065]';
for i = 1:N-1
    C = quat2rotMat(true_quat(:, i)); % nav->body frame transformation matrix
    z_acc0(:, i) = C * (-g_vec + true_acc(:, i));
    z_mag0(:, i) = C * (mag_vec + d_N(:, i));
    z_gyro0(:, i) = C * true_angvel(:, i); %true gyro bias is zero, thus not added
end

%%
% Monte Carlo Simulation
% MEKF - Allocate Memory for state predictions, error covariance diagonals
MEKFrmse_angle_errs = zeros(6, 3, 50); % #cases, #trials, #euler angles
MEKFrmse_bias_errs = zeros(6, 3, 50);
MEKFSampleTrialsAngErr = zeros(6, 3, N);
MEKFSampleTrialsBiasErr = zeros(6, 3, N);

% CI-MEKF - Allocate Memory for state predictions, error covariance diagonals
CIMEKFrmse_angle_errs = zeros(6, 3, 50); % #cases, #trials, #euler angles
CIMEKFrmse_bias_errs = zeros(6, 3, 50);
CIMEKFSampleTrialsAngErr = zeros(6, 3, N);
CIMEKFSampleTrialsBiasErr = zeros(6, 3, N);
c_acc = 3;
c_mag = 3;
true_gyro_bias = [0 0 0]';
for case_no = 1:6
    gyro_bias_err0 = biasInitErrs(:, case_no); % in degree/s
    angle_init_err0 = angleInitErrs(:, case_no); %in degrees
    for trial_no = 1:50
        gyro_bias_err = gyro_bias_err0 + randn(3,1); %add 1 degree/s random noise to the bias estimate
        angle_init_err = angle_init_err0 + sqrt(5) * randn(3,1); % add 5 degrees random noise to the orientation angle estimate
        q_est0 = angErr2InitQuatEst(deg2rad(angle_init_err), true_quat(:, 1));
        bias_est0 = deg2rad(true_gyro_bias + gyro_bias_err); % 3x1
        z_acc = z_acc0 + sqrt(var_a) * randn(3,N-1); % 3xN 
        z_mag = z_mag0 + sqrt(var_m) * randn(3,N-1); % 3xN
        z_gyro = z_gyro0 + true_gyro_bias + sqrt(var_g) * randn(3,N-1); % 3xN
        angleErr = zeros(3, N-1); %angle errors in degrees
        biasErr = zeros(3, N-1); %bias estimate errors in degrees / s

        x0 = [0 0 0 bias_est0']';
        x_ = zeros(6, 1);
        P_ = zeros(6, 6);
        x_est = x0;
        P_est = P0;
        q_est = q_est0;

        for time_step = 1:N-2
        % --- 1. PREDICT state from t to t+1 ---
            omega = z_gyro(:, time_step); 
            omega_corr = omega - x_est(4:6, :); 
            
            dq_nom = [1; 0.5 * omega_corr * dt];
            q_ = quatmultiply(q_est', dq_nom')'; 
            q_ = q_ / norm(q_);
        
            x_est(1:3, :) = zeros(3,1); 
            [x_, P_] = MEKFPredict(x_est, P_est, Q, omega_corr, dt);
            x_(1:3, :) = zeros(3,1); % Enforce zero orientation error in prediction 
            
            % --- 2. UPDATE using measurements at t+1 ---
            z_a = z_acc(:, time_step + 1);
            z_m = z_mag(:, time_step + 1);
            
            [x_est, P_est] = updateAccMEKF(x_, P_, R_a, q_, z_a);
            q_est = quatUpdateNormalize(x_est, q_)';
            x_est(1:3, :) = zeros(3,1); 
            
            [x_est, P_est] = updateMagMEKF(x_est, P_est, R_m, q_est, z_m, mag_vec);
            q_est = quatUpdateNormalize(x_est, q_est)'; 
            x_est(1:3, :) = zeros(3,1); 
            
            % --- 3. ERROR CALCULATION at t+1 ---
            q_err = quatmultiply(quatinv(q_est'), true_quat(:, time_step + 1)');
            q_err = q_err / norm(q_err);
            
            angle_err = wrapToPi(quat2eul(q_err, 'ZYX')') * 180 / pi;
            bias_err = (true_gyro_bias - x_est(4:6, :)) * 180 / pi;
            
            angleErr(:, time_step + 1) = angle_err;
            biasErr(:, time_step + 1) = bias_err;
            
            if(trial_no == 10)
                MEKFSampleTrialsAngErr(case_no, :, time_step + 1) = angle_err;
                MEKFSampleTrialsBiasErr(case_no, :, time_step + 1) = bias_err;
            end
        end
        MEKFrmse_angle_errs(case_no, :, trial_no) = errArr2RMSE(angleErr(:, 2:N-1))';
        MEKFrmse_bias_errs(case_no, :, trial_no) = errArr2RMSE(biasErr(:, 2:1:N-1))';
    end
end

%%
true_gyro_bias = [0 0 0]';
for case_no = 1:6
    gyro_bias_err0 = biasInitErrs(:, case_no); % in degree/s
    angle_init_err0 = angleInitErrs(:, case_no); %in degrees
    for trial_no = 1:50
        gyro_bias_err = gyro_bias_err0 + randn(3,1); %add 1 degree/s random noise to the bias estimate
        angle_init_err = angle_init_err0 + sqrt(5) * randn(3,1); % add 5 degrees random noise to the orientation angle estimate
        q_est0 = angErr2InitQuatEst(deg2rad(angle_init_err), true_quat(:, 1));
        bias_est0 = deg2rad(true_gyro_bias + gyro_bias_err); % 3x1
        z_acc = z_acc0 + sqrt(var_a) * randn(3,N-1); % 3xN 
        z_mag = z_mag0 + sqrt(var_m) * randn(3,N-1); % 3xN
        z_gyro = z_gyro0 + true_gyro_bias + sqrt(var_g) * randn(3,N-1); % 3xN
        angleErr = zeros(3, N-1); %angle errors in degrees
        biasErr = zeros(3, N-1); %bias estimate errors in degrees / s

        x0 = [0 0 0 bias_est0']';
        x_ = zeros(6, 1);
        P_ = zeros(6, 6);
        x_est = x0;
        P_est = P0;
        q_est = q_est0;

        for time_step = 1:N-2
            % --- 1. PREDICT state from t to t+1 ---
            omega = z_gyro(:, time_step); 
            omega_corr = omega - x_est(4:6, :); 

            dq_nom = [1; 0.5 * omega_corr * dt];
            q_ = quatmultiply(q_est', dq_nom')'; 
            q_ = q_ / norm(q_);

            x_est(1:3, :) = zeros(3,1); 
            [x_, P_] = MEKFPredict(x_est, P_est, Q, omega_corr, dt);
            x_(1:3, :) = zeros(3,1); % Enforce zero orientation error in prediction 

            % --- 2. UPDATE using measurements at t+1 ---
            z_a = z_acc(:, time_step + 1);
            z_m = z_mag(:, time_step + 1);

            [x_est, P_est] = updateAccCIMEKF(x_, P_, R_a, q_, z_a, c_acc);
            q_est = quatUpdateNormalize(x_est, q_)';
            x_est(1:3, :) = zeros(3,1); 

            [x_est, P_est] = updateMagCIMEKF(x_est, P_est, R_m, q_est, z_m, mag_vec, c_mag);
            q_est = quatUpdateNormalize(x_est, q_est)'; 
            x_est(1:3, :) = zeros(3,1); 

            % --- 3. ERROR CALCULATION at t+1 ---
            q_err = quatmultiply(quatinv(q_est'), true_quat(:, time_step + 1)');
            q_err = q_err / norm(q_err);

            angle_err = wrapToPi(quat2eul(q_err, 'ZYX')') * 180 / pi;
            bias_err = (true_gyro_bias - x_est(4:6, :)) * 180 / pi;

            angleErr(:, time_step + 1) = angle_err;
            biasErr(:, time_step + 1) = bias_err;

            if(trial_no == 10)
                CIMEKFSampleTrialsAngErr(case_no, :, time_step + 1) = angle_err;
                CIMEKFSampleTrialsBiasErr(case_no, :, time_step + 1) = bias_err;
            end
        end
        steady_state_start = round(20/dt);
        CIMEKFrmse_angle_errs(case_no, :, trial_no) = errArr2RMSE(angleErr(:, 2:N-1))';
        CIMEKFrmse_bias_errs(case_no, :, trial_no) = errArr2RMSE(biasErr(:, 2:N-1))';
    end
end

%% Graphing and Tables
% Angle error of 2 methods for specific case compared with matlab legend
% Plotting the bias errors for the first two cases
% Define the target case
target_case = 5;

% Extract and reshape the 3xN data matrices
CIMEKF_angles = squeeze(CIMEKFSampleTrialsAngErr(target_case, :, :));
MEKF_angles = squeeze(MEKFSampleTrialsAngErr(target_case, :, :));

% Generate a matching time vector based on the actual array length
N_plot = size(CIMEKF_angles, 2);
t_plot = (1:N_plot) * dt;

% --- Plot 1: CI-MEKF Errors ---
figure('Name', 'CI-MEKF Attitude Errors - Case 5', Theme="light");
subplot(2,1,1);
plot(t_plot, squeeze(CIMEKF_angles(1, :)), "b", t_plot, squeeze(CIMEKF_angles(2, :)), "g", t_plot, squeeze(CIMEKF_angles(3, :)), "r");
title('CI-MEKF Angles Errors - Case 5'); ylabel('Error (deg)'); 
legend('Yaw', 'Pitch', 'Roll');
grid on;

subplot(2,1,2);
plot(t_plot, squeeze(MEKF_angles(1, :)), "b", t_plot, squeeze(MEKF_angles(2, :)), "g", t_plot, squeeze(MEKF_angles(3, :)), "r", 'DisplayName', "Yaw", 'DisplayName', "Pitch", 'DisplayName', "Roll");
title('MEKF Attitude Errors - Case 5'); ylabel('Error (deg)'); 
legend('Yaw', 'Pitch', 'Roll');
grid on;

%%
% Angle error of 2 methods for specific case compared with matlab legend
% Plotting the bias errors for the first two cases
% Define the target case
target_case = 5;

% Extract and reshape the 3xN data matrices
CIMEKF_biases = squeeze(CIMEKFSampleTrialsBiasErr(target_case, :, :));
MEKF_biases = squeeze(MEKFSampleTrialsBiasErr(target_case, :, :));

% Generate a matching time vector based on the actual array length
N_plot = size(CIMEKF_biases, 2);
t_plot = (1:N_plot) * dt;

% --- Plot 1: CI-MEKF Errors ---
figure('Name', 'CI-MEKF Gyro Bias Errors - Case 5', Theme="light");
subplot(2,1,1);
plot(t_plot, squeeze(CIMEKF_biases(1, :)), "b", t_plot, squeeze(CIMEKF_biases(2, :)), "g", t_plot, squeeze(CIMEKF_biases(3, :)), "r");
title('CI-MEKF Bias Errors - Case 5'); ylabel('Error (deg/s)'); 
legend('Yaw', 'Pitch', 'Roll');
grid on;

subplot(2,1,2);
plot(t_plot, squeeze(MEKF_biases(1, :)), "b", t_plot, squeeze(MEKF_biases(2, :)), "g", t_plot, squeeze(MEKF_biases(3, :)), "r", 'DisplayName', "Yaw", 'DisplayName', "Pitch", 'DisplayName', "Roll");
title('MEKF Gyro Bias Errors - Case 5'); ylabel('Error (deg/s)'); 
legend('Yaw', 'Pitch', 'Roll');
grid on;

%%
% Plot Magnetic Disturbance, acceleration and angular speed
N_plot = size(CIMEKF_biases, 2);
t_plot = (1:N_plot-1) * dt;
figure('Name', 'Magnetic Disturbance', Theme="light");
plot(t_plot, squeeze(z_mag(1,:)), 'b', t_plot, squeeze(z_mag(2,:)), 'g', t_plot, squeeze(z_mag(3,:)), 'r');
title('Magnetic Disturbance - Nav Frame'); ylabel('dN (G)'); xlabel('time (s)');
legend('X', 'Y', 'Z'); grid on;

figure('Name', 'Acceleration', Theme="light");
plot(t_plot, squeeze(z_acc(1,:)), 'b', t_plot, squeeze(z_acc(2,:)), 'g', t_plot, squeeze(z_acc(3,:)), 'r');
title('Ground Truth Acceleration - Body Frame'); ylabel('aB (m/s^2)'); xlabel('time (s)');
legend('X', 'Y', 'Z'); grid on;

figure('Name', 'Angular Velocity', Theme="light");
plot(t_plot, squeeze(z_acc(1,:)), 'b', t_plot, squeeze(z_acc(2,:)), 'g', t_plot, squeeze(z_acc(3,:)), 'r');
title('Ground Truth Angular Velocity - Body Frame'); ylabel('Omega(rad/s)'); xlabel('time (s)');
legend('Yaw', 'Pitch', 'Roll'); grid on;

%%
% Plot the trajectory of the robot
plot3Position(true_pos);
%%
% Table of RMSE errors for each case
mu = mean(CIMEKFrmse_angle_errs, 3);
sigma = std(CIMEKFrmse_angle_errs, 0, 3);

figure('Name', 'CI-MEKF Angle RMSE', Theme="light");
b = bar(mu, 'grouped');
hold on;
[ngroups, nbars] = size(mu);
x = nan(ngroups, nbars);
for i = 1:nbars
    x(:, i) = b(i).XEndPoints;
end
errorbar(x, mu, sigma, 'k', 'linestyle', 'none', 'LineWidth', 1);
set(gca, 'XTickLabel', "Case " + (1:ngroups));
legend({'Yaw', 'Pitch', 'Roll'}, 'Location', 'best');
ylabel('RMSE (degrees)');
hold off;

tableData = strings(ngroups, nbars);
for i = 1:ngroups
    for j = 1:nbars
        tableData(i, j) = sprintf('%.2f ± %.2f', mu(i, j), sigma(i, j));
    end
end

Cases = "Case " + (1:ngroups)';
Yaw = tableData(:, 1);
Pitch = tableData(:, 2);
Roll = tableData(:, 3);
ResultTable = table(Cases, Yaw, Pitch, Roll);
disp(ResultTable);