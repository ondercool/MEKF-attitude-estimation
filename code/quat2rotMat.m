function C = quat2rotMat(q)
    % Extracts elements assuming scalar-first convention: q = [q0; q1; q2; q3]
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);
    
    % Constructs the rotation matrix C(q)
    C = 2 * [ (q0^2 + q1^2 - 0.5), (q1*q2 + q0*q3),     (q1*q3 - q0*q2);
        (q1*q2 - q0*q3),     (q0^2 + q2^2 - 0.5), (q2*q3 + q0*q1);
        (q1*q3 + q0*q2),     (q2*q3 - q0*q1),     (q0^2 + q3^2 - 0.5) ];
end