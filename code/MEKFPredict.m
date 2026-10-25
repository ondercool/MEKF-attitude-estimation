function [x_,P_] = MEKFPredict(x, P, Q, w_hat, dt)
    x_ = [0 0 0 x(4:6, :)']';
    F = cont2DistTfMat(w_hat, dt);
    P_ = F*P*F' + Q;
end