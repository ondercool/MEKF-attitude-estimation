function [exp_mat] = skewExponential(A, t)
    % Extract the magnitude of the angular velocity directly from the skew matrix
    w_norm = sqrt(A(3,2)^2 + A(1,3)^2 + A(2,1)^2);
    
    if w_norm > 1e-6
        % General Rodrigues' rotation formula
        exp_mat = eye(3) + (sin(w_norm * t) / w_norm) * A + ((1 - cos(w_norm * t)) / w_norm^2) * (A * A);
    else
        % First-order approximation to prevent divide-by-zero for near-zero rotations
        exp_mat = eye(3) + A * t;
    end
end