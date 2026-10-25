function d_N = generateMagneticDisturbance(t, stepTimes, durations, magnitudes)
    N = length(t);
    d_N = zeros(3, N); % Initialize the disturbance array
    for i = 1:length(stepTimes)
        startIdx = find(t >= stepTimes(i), 1);
        endIdx = startIdx + durations(i) - 1;
        if endIdx <= N
            d_N(1, startIdx:endIdx) = d_N(1, startIdx:endIdx) + magnitudes(1);
            d_N(2, startIdx:endIdx) = d_N(2, startIdx:endIdx) + magnitudes(2);
            d_N(3, startIdx:endIdx) = d_N(3, startIdx:endIdx) + magnitudes(3);
        end
    end
end