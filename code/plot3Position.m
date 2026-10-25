function h = plot3Position(pos)
    % plot3Position  Plot 3-D trajectory from pos (3xN)
    %   h = plot3Position(pos) plots the columns of pos as [x;y;z] points.
    %   Returns a handle to the main plot line.
    %
    %   Example:
    %     % pos is 3xN
    %     plot3Position(pos);
    
    validateattributes(pos, {'numeric'}, {'2d','nrows',3,'nonempty'}, mfilename, 'pos', 1);
    
    x = pos(1, :);
    y = pos(2, :);
    z = pos(3, :);
    
    figure('Name', 'Flight Trajectory', Theme="light");
    h = plot3(x, y, z, '-b', 'LineWidth', 1.5);
    hold on;
    plot3(x(1), y(1), z(1), 'go', 'MarkerFaceColor', 'g'); % start
    plot3(x(end), y(end), z(end), 'ro', 'MarkerFaceColor', 'r'); % end
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('3-D Trajectory');
    grid on;
    axis equal;
    view(3);
    hold off;
end