function [breg y_noseas y_fit] = fitseasoncycle(y, time);
% [breg y_noseas] = fitseasoncycle(y, time);
%
% Function to estimate the seasonal cycle from a time series.
%
% Inputs:
% 
% + y: data vector
% + time: time vector (must be in years)
%
% Outputs:
%
% + breg: least-square coefficients [mean trend annual_1 annual_2
% semiannual_1 semiannual_2]
% + y_noseas: data vector minus the seasonal cycle (annual + semiannual).
% + y_fit: modeled data
% 
% Manuel O. Gutierrez-Villanueva          2022/07/06

if length(y)~=length(time)
    error('Data and time vector are not the same size')
end

y = y(:);
time = time(:);

yy = y(~isnan(y)); % remove nans
xx = time(~isnan(yy));

% Ordinary least-squares
phianual = 2*pi*xx(:);
phisemi = 2*phianual;

Xfull = [ones(length(xx), 1) xx(:) cos(phianual(:)) sin(phianual(:)) ...
    cos(phisemi(:)) sin(phisemi(:))];
breg = regress(yy(:), Xfull);

% Removes seasonal cycle
y_noseas = y;
y_noseas(~isnan(y)) = yy(:) - Xfull(:, 3:end)*breg(3:end);

y_fit = Xfull*breg;