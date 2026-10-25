function [A] = Omega2AMat(w_hat)
    w_hat_skew = skewSymmetric(w_hat);
    A = [-w_hat_skew eye(3); zeros(3) zeros(3)];
end