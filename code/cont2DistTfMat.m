function [F] = cont2DistTfMat(w_hat, T)
    w_skew = skewSymmetric(w_hat);
    w_norm = norm(w_hat);
    
    % Prevent divide-by-zero during low angular velocities
    if w_norm > 1e-6
        % Correct Rodrigues formula for exp(-[w x] T)
        exp_w_skew = eye(3) - (sin(w_norm*T)/w_norm)*w_skew + ((1-cos(w_norm*T))/(w_norm^2))*(w_skew*w_skew);
        
        % Correct Integral term (Phi_2)
        term1 = eye(3) * T;
        term2 = -((1-cos(w_norm*T))/(w_norm^2)) * w_skew; % Fixed sign
        term3 = ((w_norm*T - sin(w_norm*T))/(w_norm^3)) * (w_skew*w_skew); % Fixed parentheses
        
        exp_int_term = -0.5 * (term1 + term2 + term3);
    else
        % Limit approximation for near-zero rotation
        exp_w_skew = eye(3) - w_skew * T;
        exp_int_term = -0.5 * eye(3) * T;
    end
    
    % Assemble F matrix using eye(3) for the bias block
    F = [exp_w_skew, exp_int_term; 
         zeros(3,3), eye(3)]; 
end