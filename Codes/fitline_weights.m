function [breg y_fit] = fitline_weights(y, w, time);
% [breg y_fit] = fitline_weights(y, time);
%
% Function to fit a line using weighted least-squares.
%
% Inputs:
% 
% + y: data vector
% + w: weights vector
% + time: time vector (must be in years)
%
% Outputs:
%
% + breg: least-square coefficients [mean trend annual_1 annual_2
% semiannual_1 semiannual_2]
% + y_fit: modeled data
% 
% Manuel O. Gutierrez-Villanueva          2022/07/06

if length(y)~=length(time)
    error('Data and time vector are not the same size')
elseif length(y)~=length(w)
    error('Data and weights vector are not the same size')
elseif length(w)~=length(time)
    error('Time and weights vector are not the same size')
end

y = y(:);
time = time(:);
w = w(:);

% remove nans
yy = y(~isnan(y)); 
xx = time(~isnan(y));
ww = w(~isnan(y));

Xmodel = [ones(size(yy)) xx];

% Weights
I = eye(length(ww));
W = I; W(boolean(eye(sum(~isnan(yy))))) = 1./ww.^2;


% Coefficients
breg = inv(Xmodel'*W*Xmodel)*Xmodel'*W*yy;
y_fit = Xmodel*breg;