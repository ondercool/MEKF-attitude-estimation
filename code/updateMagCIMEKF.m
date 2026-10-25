function [x_, P_] = updateMagCIMEKF(x, P, R, q_hat, p_bm, p_N, c)
    % EKF Measurement Update
    % Actual mn value in the navigation frame is calculated beforehand with
    % the magnetic field model
    p_B_hat = quat2rotMat(q_hat) * p_N;
    p_B_hat_skew = skewSymmetric(p_B_hat);
    H = [2*p_B_hat_skew zeros(3,3)];
    y = p_bm - p_B_hat;
    y_res = y - H*x;
    S = H*P*H' + R + c*diag(y_res.^2);
    K = P * H' / S; % Compute the Kalman gain
    x_ = x + K * y_res; % Update the state estimate
    P_ = (eye(size(P)) - K * H) * P * (eye(size(P)) - K * H)' + K*R*K'; % For numerical stability
end