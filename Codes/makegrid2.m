function [gpos,xr,yr] = makegrid2(lat, lon, binlength,press,theta,lat0,lon0);
% [gpos,xr,yr] = makegrid2(lat, lon, binlength, press, theta, lat0, lon0) 
% uses the original lat, lon defined by vectors [lat, lon] to compute 
% grid positions separated by 'binlength' (km) in both the projected 
% cartesian coordinates and lat, lon.
%
% grid point positions are in gpos, while xr and yr are original 
% positions in passage coordinate system
% 
% theta is optional and defaults to YDL's from DPshort
% 
% Modified by MOGV (08/19/2019): Allows to specifies pressure, theta, lat0
% and lon0 for the along-across DP grids.


%transform lat, lons into cartesian co-ordinates
[x,y] = ll2xy2(lat, lon,lat0,lon0, press); 

%this is the angle to subtract--rotate clockwise by
%   (YLF NOTES) so the angle from zonal to the new x
%theta= 16;
% if length(varargin)==0
%    load (['./data/all_drake/DPshort' 'theta']);
     %(YLF NOTES) x_pc is 23 degrees ccw from zonal
% else
%    theta = varargin{1};
% end


if theta ~= 0;
    xr = real(complex(x,y)*exp(-i*theta/180*pi));
    yr = imag(complex(x,y)*exp(-i*theta/180*pi));
else
    xr = real(complex(x,y));
    yr = imag(complex(x,y));
end
%define x and y grid points
gpos.x = floor(min(xr)/binlength)*binlength:binlength:...
    floor(max(xr)/binlength)*binlength; 
gpos.y = ceil(min(yr)/binlength)*binlength:binlength:...
    floor(max(yr)/binlength)*binlength;

[xg, yg] = meshgrid(gpos.x, gpos.y);

if theta ~= 0
    x0 = real(complex(xg,yg)*exp(i*theta/180*pi));
    y0 = imag(complex(xg,yg)*exp(i*theta/180*pi));
else
    x0 = real(complex(xg,yg));
    y0 = imag(complex(xg,yg));

end

[gpos.lat, gpos.lon] = xy2ll2(x0,y0,lat0,lon0,press);
