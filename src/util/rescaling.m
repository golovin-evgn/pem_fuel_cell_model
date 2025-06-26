function y = rescaling(x, x_min, x_max, y_min, y_max)
%% rescaling(in, x_min, x_max, y_min, y_max) scales input vector in from
% range [x_min, x_max] to output range [y_min, y_max]

y = y_min + (y_max-y_min) .* (x - x_min) ./ (x_max - x_min);