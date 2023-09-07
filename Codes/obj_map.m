function [datamap Emap] = obj_map(xx1, yy1, data, decorr_x, decorr_y, err, xg, yg);
% [datamap Emap] = obj_map(xx1, yy1, data, decorr_x, decorr_y, err, xg, yg)
%
% Function for objectively map data
%
% Manuel O. Gutierrez Villanueva
% 2021/08/24

if length(xx1) ~= length(yy1);
    error('xx1 is not the same length as yy1');
elseif length(xx1) ~= length(data) || length(yy1) ~= length(data);
    error('data is not the same length as xx1 and yy1');
end

if length(decorr_x) ~= 1 || length(decorr_y) ~= 1 || length(err) ~= 1
    error('Decorrlation scale and error must be a constant');
end

if err > 1 || err <= 0;
    error('Mapping error must be 0 < Error <= 1');
end

%Calculate the signal covariance matrix E plus noise.
[Xx1 Xx2] = meshgrid(xx1);
[Yy1 Yy2] = meshgrid(yy1);
Dx = Xx1 - Xx2;
Dy = Yy1 - Yy2;
% dxr=dx*cosphi-dy*sinphi;
% dyr=dx*sinphi+dy*cosphi;
Ee = exp(-Dx.^2/(decorr_x)^2 - Dy.^2/(decorr_y)^2);
Ee = Ee + err*eye(length(xx1));%covariance matrix plus error

% Calculate my covariance matrix
Cc = zeros(length(xx1),length(xg(:)));
% [Yy,Xx] = meshgrid(yg(:)',xg(:)');

clear n

for n = 1:length(xx1);
    Dx2(n,:) = xg(:)'- xx1(n);
    Dy2(n,:) = yg(:)'- yy1(n);
%    dxr = dx*cosphi - dy*sinphi;
%    dyr=dx*sinphi+dy*cosphi;
    Cc(n,:) = exp(-Dx2(n,:).^2/(decorr_x)^2 - Dy2(n,:).^2/(decorr_y)^2);
end

EC1 = Ee\Cc;%inverse of data-signal covariance matrix
% tprime = temp_prime(r_nonan1);


datamap = EC1'*data;%objective mapping
%removes mean a bit dirty-naughty


% Now calculates the normalized error
Emap = 1 - diag(Cc'*EC1); % psi

% Error1 = diag(1 - Cc'*EC1);%/nanvar(temp_gd));