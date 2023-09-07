% [x,y] = ll2xy(lat,lon,lat0,lon0, pressure);
% convert from lat-lon coordinates to xy coordinates, with respect to some origin point
%     lat      = decimal degrees (+ve N, -ve S) [- 90.. +90]
%     lon      = decimal degrees (+ve E, -ve W) [-180..+180]
% OPTIONAL:
%     pressure     =  sea pressure ( default is 0 )               [ dbar ]
%            ( i.e. absolute pressure - 10.1325 dbar )
% Modified by MOGV (08/20/2019): Units is removed and now the last input is
% pressure.
function[x,y] = ll2xy2(lat,lon,lat0,lon0,press);

lon(lon>180) = lon(lon>180)-360;
if intersect(size(lat), size(lon)) == 1 & length(lat)~=length(lon);
    [lon,lat] = meshgrid(lon,lat);
end


xlon = lon0*ones(length(lon(:))*2,1);
xlon(2:2:end) = lon(:); % intersperse longitudinal series with reference longitude
xlat = repmat(lat(:)',2,1)'; % set latitude pairs to match data-reference pairs for longitude. 

xlon = reshape(xlon,2,length(lon))';
x = gsw_distance(xlon,xlat,press); % calculate cartesian x co-ordinates
% x = x(1:2:end); % select pairs with correct latitudes
ind_xng = find(lon < lon0);
x(ind_xng) = -x(ind_xng); % check signs.

ylat = lat0*ones(length(lat(:))*2,1); 
ylat(2:2:end) = lat(:); %intersperse lat data with reference latitude.
ylat = reshape(ylat,2,length(lon))';

ylon = repmat(lon,1,2);

y = gsw_distance(ylon,ylat,press); % calculate cartesian y co-ordinates, independant of longitude
% y = y(1:2:end); % select only neccesary data
ind_yneg = find(lat<lat0);
y(ind_yneg) = -y(ind_yneg);%check signs

% x = reshape(x, size(lat,1), size(lat,2));
% y = reshape(y, size(lat,1), size(lat,2));
