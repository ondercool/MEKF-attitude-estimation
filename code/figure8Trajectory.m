function [pos, q, vel, acc, angvel] = figure8Trajectory(t_end, dt, R, H, dH)
    % 1. Simulation Parameters
    t = (0:dt:t_end)';
    omega = 2*pi / t_end; 
    N = length(t);
    
    % 2. Position
    North = R * sin(omega * t);
    East  = R * sin(omega * t) .* cos(omega * t);
    Down  = H - dH * sin(2 * omega * t);
    pos = [North, East, Down]'; % 3xN
    
    % 3. Velocity
    v_North = R * omega * cos(omega * t);
    v_East  = R * omega * (cos(omega * t).^2 - sin(omega * t).^2);
    v_Down  = -2 * dH * omega * cos(2 * omega * t);
    vel = [v_North, v_East, v_Down]'; % 3xN
    
    % 4. Analytical Acceleration (Bypassing noisy splines entirely)
    a_North = -R * omega^2 * sin(omega * t);
    a_East  = -2 * R * omega^2 * sin(2 * omega * t);
    a_Down  = 4 * dH * omega^2 * sin(2 * omega * t);
    acc = [a_North, a_East, a_Down]'; % 3xN
    
    % 5. Orientation
    horizontal_speed = sqrt(v_North.^2 + v_East.^2);
    Yaw   = atan2(v_East, v_North);
    Pitch = atan2(-v_Down, horizontal_speed); 
    
    % Fix radians: 34 degrees is exactly ~0.593 radians
    Roll  = 0.2 * sin(omega * t); 
    
    true_orientations = quaternion([Yaw, Pitch, Roll], 'euler', 'ZYX', 'frame');
    q = compact(true_orientations)'; % 4xN
    
    % 6. Stable Angular Velocity (Avoids waypointTrajectory noise)
    % Compute smooth Euler rates
    % Pre-calculate squared terms and horizontal magnitude
    V_horizontal_sq = v_North.^2 + v_East.^2;
    V_horizontal = sqrt(V_horizontal_sq);
    V_total_sq = V_horizontal_sq + v_Down.^2;

    % Derivative of V_horizontal (using chain rule)
    V_horizontal_dot = (v_North .* a_North + v_East .* a_East) ./ V_horizontal;

    % 1. Exact Yaw Rate (derivative of atan2(v_East, v_North))
    Yaw_dot = (v_North .* a_East - v_East .* a_North) ./ V_horizontal_sq;
    % 2. Exact Pitch Rate (derivative of atan2(-v_Down, V_horizontal))
    Pitch_dot = (V_horizontal .* (-a_Down) - (-v_Down) .* V_horizontal_dot) ./ V_total_sq;
    % 3. Exact Roll Rate (derivative of 0.593 * sin(omega * t))
    Roll_dot = 0.2 * omega * cos(omega * t);
    angvel = zeros(3, N);
    for i = 1:N
        phi = Roll(i);
        theta = Pitch(i);
        
        % Map Euler rates to Body frame (Standard ZYX mapping)
        W = [1,  0,        -sin(theta);
             0,  cos(phi),  sin(phi)*cos(theta);
             0, -sin(phi),  cos(phi)*cos(theta)];
             
        euler_rates = [Roll_dot(i); Pitch_dot(i); Yaw_dot(i)];
        angvel_body = W * euler_rates;
        
        % Rotate back to Navigation frame so your CIMEKFSimulation 
        % (z_gyro0 = C * true_angvel) processes it correctly
        C_nav2body = quat2rotMat(q(:, i));
        angvel(:, i) = C_nav2body' * angvel_body;
    end
end