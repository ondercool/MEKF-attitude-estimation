function [Q] = Variances2QMat(var_g, var_b, T)
    Phi1 = 0.25*(var_g + 1/3 * var_b) * T;
    Phi2 = 0.25*var_b*T^2;
    Phi3 = var_b * T;
    Q = [Phi1*eye(3) Phi2*eye(3); Phi2*eye(3) Phi3*eye(3)];
end