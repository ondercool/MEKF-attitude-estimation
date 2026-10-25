function [x_, P_] = updateAccMEKF(x, P, R, q_hat, p_bm)
    % EKF Measurement Update
    p_N = [0 0 -9.8065]';
    p_B_hat = quat2rotMat(q_hat) * p_N;
    p_B_hat_skew = skewSymmetric(p_B_hat);
    H = [2*p_B_hat_skew zeros(3,3)];
    y = p_bm - p_B_hat;
    y_res = y - H*x;
    S = H*P*H' + R;
    K = P * H' / S; % Compute the Kalman gain
    x_ = x + K * y_res; % Update the state estimate
    P_ = (eye(size(P)) - K * H) * P * (eye(size(P)) - K * H)' + K*R*K'; % For numerical stability
end
