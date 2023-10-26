% Code to estimate streamwise binned and averaged eddy momentum fluxes from
% the OS38 ADCP.
%
% MOGV 01/29/2023

clear all;
close all;

pycmd = '/Users/manuelgutierrez/anaconda3/bin/python';

% Parameters
gsize = 25; %km
cont = [-1.0:0.1:1.0];%vel contours
contmaperr = [0:0.1:0.5];%mapping error contours
decorr_x = 50000; %meters distance
error = 0.2;

% loads colormaps
load('./ekeemf_cmap.mat')


adcp = 'os38nb'; % os38nb or nb150
if strcmp(adcp,'nb150') == 1
    maxz = 210; %maximum depth
    decorr_y = 100; %meters depth
else
    maxz = 780;
    decorr_y = 300; %meters depth
end
varia = 1;% Variance 
error = 0.2;% noise-to-signal ratio
op_plot_m = 0;% Option to plot geostrophic streamfunctions
boxlim = [0    1000   -maxz   0]; %Axis limits
Umean = 91.21e6;%Mean transport to Dec2019

arch = ['../Datasets/lmgvel' adcp '_' num2str(decorr_x/1000) '_' num2str(decorr_y)];


% loads transects
load(['../BigDatasets/lmg_' adcp '_rotated_25pg.mat']);


% Loads bathymetry to determine which cruises are good for transport
b = load('../Datasets/drakembf_avgd'); 
b.lat = [b.lats; b.lat; b.latn]; 
b.lon = [b.lonw b.lon]; 
b.bathyf = [b.bathyfw [b.bathyfs; b.bathyf; b.bathyfn]];

% Bathymetry boundaries
dpn = [-55, -54.70];
dps1 = [-0.00954, 5.9987, -1002.8768];
dps2 = [-0.00988, 6.2136, -1036.9253];

% loads good transects
good = load(['../Datasets/all_transects_' adcp '.txt']);
flag = good(good(:, 1)<3 | good(:, 1)==6, 1);
good = datenum(good(good(:, 1)<3 | good(:, 1)==6, 2:end));
[good, iip] = sort(good, 'ascend');
flag = flag(iip);


lat0 = -55;
lon0 = -65;
omega = 2*pi/86400;
overlapping = 0.0;% Overlapping 0.5 = 50%
depths_psi = [0 250];% Averages mean geostrophic velocities and psi over these two depths


% Splits into transects
I = interval(all_lmg.t, 2, 'n');
fprintf('\n\nNumber of transects available: %3.0f \n\n', length(I));
fprintf('Starting date: %s \n', datestr(all_lmg.t(1)));
fprintf('Ending date: %s \n\n', datestr(all_lmg.t(end)));

for i = 1:length(I);
    meantime(i, 1) = nanmean(all_lmg.t(I{i}));
end

[x y] = meshgrid(meantime, good);
[t1, ~] = meshgrid(1:length(meantime), 1:length(good));

% Finds only good transects
inds = abs(x-y)<1;
inds = t1(inds);


clear x y t1

% Removes tides
all_lmg.uv = all_lmg.uv - repmat(complex(tide.u, tide.v)', ...
    size(all_lmg.uv,1), 1);

% Removes data below specific depth
all_lmg.uv = all_lmg.uv(all_lmg.z<=maxz,:);
all_lmg.z = all_lmg.z(all_lmg.z<=maxz);
all_lmg.lon = all_lmg.lon - 360;
dz = nanmean(diff(all_lmg.z));
nz = length(all_lmg.z);

% Rotates into x,y (distance) grid
% load gridded_mean_uv_upd0;
load(['../Datasets/DPshort.mat'],'theta');% Angle to rotate

[all_x, all_y] = ll2xy2(all_lmg.lat, all_lmg.lon, lat0, lon0, 0);

all_xy = complex(all_x, all_y).*exp(-sqrt(-1)*theta*pi/180);
all_x = real(all_xy); all_y = imag(all_xy);


% Loads OM mean currents (2004 - 2016)
load(['../Datasets/omgeosvel_' adcp '.mat']);
om.uvmean = smean_om3.uvmean(smean_om3.z>22 & smean_om3.z<=maxz, :, :); % depth x lat x lon
om.uvmean = om.uvmean.*exp(sqrt(-1)*theta*pi/180);% rotates the mean vectors
om.uvmean_err = smean_om3.err(smean_om3.z>22 & smean_om3.z<=maxz, :, :);
om.lon_mean = smean_om3.lon;
om.lat_mean = smean_om3.lat;
om.z_mean = smean_om3.z(smean_om3.z>22 & smean_om3.z<=maxz);
om.uvmean = permute(om.uvmean, [2, 3, 1]); % permutes to lat x lon x depth
om.uvmean_err = permute(om.uvmean_err, [2, 3, 1]); % permutes to lat x lon x depth
om.x = smean_om3.x;
om.y = smean_om3.y;

clear smean_om3 sm3_lmg gadcp_lmg

all_lmg.uvprime = nan(size(all_lmg.uv));

% Removes mean from total velocity vector
for i = 1:size(om.lon_mean, 1) % along latitude
    for j = 1:size(om.lon_mean, 2) % along longitude

        if ~isnan(squeeze(abs(om.uvmean(i, j, 1))));
            indx = find((abs(om.x(i, j) - all_x) < gsize*1e3/2) & ...
                (abs(om.y(i, j) - all_y) < gsize*1e3/2));

            if ~isempty(indx) % Removes mean
                all_lmg.uvprime(:, indx) = all_lmg.uv(:, indx) - ...
                    squeeze(om.uvmean(i, j, :));
               
                
            end
            clear indx
        end
    end
end

% Depth integrated eke and u'v'
eke = 0.5*(real(all_lmg.uvprime).^2 + imag(all_lmg.uvprime).^2);
all_lmg.eke = sum(eke, 1)*mean(diff(all_lmg.z)) + ...
    eke(1, :)*[all_lmg.z(1)-mean(diff(all_lmg.z))/2];
uv_eddy = real(all_lmg.uvprime).*imag(all_lmg.uvprime);
all_lmg.uv_prime = sum(uv_eddy)*mean(diff(all_lmg.z)) + ...
    uv_eddy(1, :)*[all_lmg.z(1)-mean(diff(all_lmg.z))/2];
all_lmg.doc = 'eke and uv_eddy are depth integrated over the upper 780 m';


% Preallocates data
distlmg = nan(length([12.5:12:1200.5]), length(inds));
latlmg = distlmg;
lonlmg = distlmg;
ekelmg = nan(length([12.5:12:1200.5]), length(inds));
uv_primelmg = ekelmg;
stheta = []; cs = []; Us = []; uoffs = []; nb = []; 

r = 0;
s = 0;

for i = 1:length(inds);
    
    clear data
    
    % Data for each transect
    data.uv_prime = all_lmg.uv_prime(I{inds(i)});
    data.eke = all_lmg.eke(I{inds(i)});
    data.time = all_lmg.t(I{inds(i)});
    data.lon = all_lmg.lon(I{inds(i)});
    data.lat = all_lmg.lat(I{inds(i)});
    data.bi = interp2(b.lon, b.lat, b.bathyf, data.lon-360, data.lat);
   
    % Calculates maximum distance (i.e. between minimum and maximum
    % transect)
    maxdist = sw_dist([nanmax(data.lat) nanmin(data.lat)], ...
        [data.lon(data.lat == nanmax(data.lat)) ...
        data.lon(data.lat == nanmin(data.lat))], 'km');

    % Sorts bottom data
    [latsort ii] = sort(data.lat);
    bisort = data.bi(ii);
    
%     % Checks if data is within some minimum bottom depth and if the
%     % transect length is equal or larger than 700 km
%     if (bisort(1) >= -1000) & (bisort(end) >= -1000) & (maxdist >= 700);
        
%     % Checks if data has a gap larger than 150 km.
%     data.dist = sw_dist(latsort, data.lon(ii), 'km');
%     dxdist = diff(data.dist);
        
%         if isempty(find(abs(dxdist) >= 75))
            
    s = s + 1;

    % Removes data away from boundaries
    dix = (data.lat < dpn(1)) & (data.lat > polyval(dps2, data.lon))...
        & (data.lon > 293.5-360);
    data.lat = data.lat(dix);
    data.lon = data.lon(dix);
    data.uv_prime = data.uv_prime(dix);
    data.eke = data.eke(dix);
    data.time = data.time(dix);
    data.bi = data.bi(dix);
    data.z = all_lmg.z;

        % If gaps are larger than some distance, make some fake data
%         data = fillgapslmg(data, 0.15, 'latitude');
    data.bi = interp2(b.lon, b.lat, b.bathyf, data.lon, data.lat);
%     end

    
    % If there gaps between 25 km < gap < 75 km,
    % fill with NaNs and linearly-interpolated latiudes/longitudes
    % for the objective mapping
            
    dist = [sw_dist(data.lat', data.lon', 'km')];
    indgap = find(dist>25);
       
    % Calculates again alongtrack distance
    [data.dist data.ang] = sw_dist(data.lat, data.lon, 'km');
    

    % Checks if it transect is south or northbound, and corrects
    % for the reversing orientation
    if (data.lat(1) > -60);% & (sum(diff(data.lat)>0)>0); %southbound
        if sum(diff(data.lat)>0)>0;
            data.dist(diff(data.lat) > 0) = ...
                data.dist(diff(data.lat)>0)*-1;
        end
        northbound(s, 1) = 0;
        nanbath = 'end';

        % Removes data at the southern end that crosses the 500 m
        % isobath
        indbi = data.bi>=-500 & data.lat < -60;
        indbi = find(indbi==1);

        if ~isempty(indbi);
            data.eke(indbi(1):end) = [];
            data.uv_prime(indbi(1):end) = [];
            data.lon(indbi(1):end) = [];
            data.lat(indbi(1):end) = [];
            data.bi(indbi(1):end) = [];
            data.time(indbi(1):end) = [];
%             data.uv(:, indbi(1):end) = [];
%             data.dist(indbi(1):end) = [];
%             data.uvmean(:, indbi(1):end) = [];
        end
%                 
     elseif (data.lat(1) < -60);% & (sum(diff(data.lat)<0)>0); %northbound
         if sum(diff(data.lat)<0)>0;
             data.dist(diff(data.lat) < 0) = ...
                 data.dist(diff(data.lat)<0)*-1;
         end
                
%          % Positive towards the east*
%          data.ur = data.ur*-1;
%          data.ur_mean = data.ur_mean*-1;
         northbound(s, 1) = 1;
         nanbath = '1';

         % Removes data at the southern end that crosses the 500 m
         % isobath
         indbi = data.bi>=-500 & data.lat < -60;
         indbi = find(indbi==1);

         if ~isempty(indbi);
             data.eke(1:indbi(end)) = [];
             data.uv_prime(1:indbi(end)) = [];
             data.lon(1:indbi(end)) = [];
             data.lat(1:indbi(end)) = [];
             data.bi(1:indbi(end)) = [];
             data.time(1:indbi(end)) = [];
%              data.uv(:, 1:indbi(end)) = [];
%              data.dist(1:indbi(end)) = [];
%              data.uvmean(:, 1:indbi(end)) = [];
         end
    end

    % Cumulative distance
    [data.dist ~] = sw_dist(data.lat, data.lon, 'km');
    data.dist = [0; cumsum(data.dist)];

    % Bin averages every certain distance
    [bins, lataver] = line_bin(data.dist, data.lat', 12);
    [bins, lonaver] = line_bin(data.dist, data.lon', 12);
    [bins, uv_prime] = line_bin(data.dist, data.uv_prime, 12);
    [bins, eke] = line_bin(data.dist, data.eke, 12);
%     [bins, biaver] = line_bin(data.dist, data.bi', gsize);
    lataver = lataver(1,:);
    lonaver = lonaver(1,:);
            
    
    distlmg(1:length(bins), s) = bins;
    lonlmg(1:length(bins), s) = lonaver;
    latlmg(1:length(bins), s) = lataver;
    lmgdata.time(1, s) = nanmean(data.time);
    lmgdata.flag(1, s) = flag(i);
    ekelmg(1:length(bins), s) = eke;
    uv_primelmg(1:length(bins), s) = uv_prime;

    
end
                
% Saves data 
lmgdata.eke = ekelmg;
lmgdata.uv_prime = uv_primelmg;
lmgdata.dist = distlmg;
lmgdata.lon = lonlmg;
lmgdata.lat = latlmg;
lmgdata.northbound = northbound;
lmgdata.adcp = adcp;
lmgdata.doc =  {'eke is depth-integrated EKE (0-780 m)';
                'uv_prime is u_prime*v_prime integrated (0-780 m)';
                'dist is along track distance [km]';
                'lon and lat are coordinates averaged for each 25 km bin';
                'northbound = 1 is a transect going from south to north';
                'NOTE: northbound transects are already flipped'};

% Remove data at the end of transects
ind = find(floor(lmgdata.time) == datenum(2014, 11, 20));
indr = find(lmgdata.lat(:, ind)<-62);

lmgdata.lat(indr, ind) = NaN;
lmgdata.lon(indr, ind) = NaN;
lmgdata.dist(indr, ind) = NaN;
lmgdata.uv_prime(indr, ind) = NaN;
lmgdata.eke(indr, ind) = NaN;

% Makes the transect that start from south to north to now be north to
% south
indsouth = find(lmgdata.northbound == 1);
indnorth = find(lmgdata.northbound == 0);


for i = 1:length(indsouth);
    [lmgdata.dist(:, indsouth(i)) ii] = ...
        sort(abs(nanmax(lmgdata.dist(:, indsouth(i))) - ...
        lmgdata.dist(:, indsouth(i))) + 12.5);
    lmgdata.eke(:, indsouth(i)) = lmgdata.eke(ii, indsouth(i));
    lmgdata.uv_prime(:, indsouth(i)) = lmgdata.uv_prime(ii, indsouth(i));
    lmgdata.lon(:, indsouth(i)) = lmgdata.lon(ii, indsouth(i));
    lmgdata.lat(:, indsouth(i)) = lmgdata.lat(ii, indsouth(i));
end


            
% Saves data
save('eke_uvprime_all_os38.mat', 'lmgdata');

