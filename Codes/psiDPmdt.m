% Program to extract MDT from Maximenko et al. (2012). and interpolate it
% to the along-across Drake Passage grid. Grid is a 25 km x 25 km
% horizontal grid used to estimate mean temperature and objectively mapped
% velocities.
%
% Inputs:
%
% + dir_folder: Specificy directory where MDT is saved.
% + lat,lon: Latitude, longitude used to interpolate the MDT.
% + mean_psi: Objectively mapped mean geostrophic streamfunction in
%             meters. Dimmensions of lat x lon.
%
% Created by MOGV 09/06/2019


function psi_mean = psiDPmdt(dir_folder,lat,lon,mean_psi);

if size(lat,1) ~= size(mean_psi,1);
    error(' Mean streamfunction and lat,lon do not have the same dimensions');
end

if size(lon,2) ~= size(mean_psi,2);
    error(' Mean streamfunction and lat,lon do not have the same dimensions');
end

psi = mean_psi;

% Extracts Mean Dynamic Topography (MDT) and reduces to DP area
mdt_niimax = ncread([dir_folder 'niiler_maximenko_1992_2012_mdt.nc'],...
    'MDOT')*1e-2;
lon_niimax = ncread([dir_folder 'niiler_maximenko_1992_2012_mdt.nc'],...
    'LONN319_N207');
lat_niimax = ncread([dir_folder 'niiler_maximenko_1992_2012_mdt.nc'],...
    'LAT37_109');

% Drake Passage area
wlon = -90; elon = -35; 
ind_lon = find(lon_niimax > wlon & lon_niimax < elon);

slat = -80; nlat = -40; 
ind_lat = find(lat_niimax > slat & lat_niimax < nlat);

lon_niimax = lon_niimax(ind_lon);
lat_niimax = lat_niimax(ind_lat);
mdt_niimax = mdt_niimax(ind_lon,ind_lat);
mdt_niimax = permute(mdt_niimax,[2,1]);

% Interpolates to sADCP grid
mdt_niimax_new = interp2(lon_niimax,lat_niimax,mdt_niimax,...
    lon,lat);


% Now interpolates the maps to the location of the objectively mapped
% streamfunction
ind_geos_nan = find(isnan(psi) == 1);
psi(ind_geos_nan) = mdt_niimax_new(ind_geos_nan) - 0.6;

psi_mean.psi = psi;
psi_mean.lon = lon;
psi_mean.lat = lat;
psi_mean.doc = char('Combines MDT of Maximenko et al. 2014 with the ',...
    'objectively mapped mean streafunction. The MDT is used where the ',...
    'mean OM mean streamfunction is not defined (i.e. away from the area',...
    'enclosed by the most repeated transects.',...
    '.psi is the combined streamfunction interpolated to lon,lat');


