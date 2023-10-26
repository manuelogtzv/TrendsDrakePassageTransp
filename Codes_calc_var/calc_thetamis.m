% Calculates transduce misalignment angle following Firing et al. (2011).
% Program to load different ADCP transects from the LMG starting from 2004
% till 2019. Program seperates data into transects and separates them into 
% southbound and northbound transects
% 
% Fills gaps <150 km in the alongtrack distance with an objective mapping.
% ONLY the across-transect velocity is employed.
%
% Manuel O. Gutierrez V.
% 2022/04/20 - Calculates misalignment angle following Firing et al 2011.

clear all; 
close all;

% Parameters
gsize = 25; %km
cont = [-1.0:0.1:1.0];%vel contours
contmaperr = [0:0.1:0.5];%mapping error contours
decorr_x = 50000; %meters distance
error = 0.2;
cmap = getPyPlot_cMap('seismic', length(cont)-1);
cmaperr = flipud(getPyPlot_cMap('pink', length(contmaperr)-1));
% cmaperrvar = flipud(getPyPlot_cMap('YlGnBu', length(contmaperr)-1));

adcp = 'os38nb'; % os38nb or nb150
if strcmp(adcp,'nb150') == 1
    maxz = 210; %maximum depth
    decorr_y = 100; %meters depth
else
    maxz = 970;
    decorr_y = 300; %meters depth
end
varia = 1;% Variance 
error = 0.2;% noise-to-signal ratio
op_plot_m = 0;% Option to plot geostrophic streamfunctions
boxlim = [0    1000   -maxz   0]; %Axis limits
Umean = 92.80e6;%Mean transport to Dec2019
dir1 = ['./Figures/OMTransects_flagged/' adcp];

% Dir for the ADCP data
path1 = ['./drake_mean_2021/data/all_drake/'];

% load gridded_mean_uv_upd0;
load([path1 'DPshort.mat'],'theta');% Angle to rotate

% loads transects
load(['./drake_mean_2021/data/all_drake/lmg_' adcp '.mat']);

% loads sla data
load([path1 'drake_sla_1993_2019.mat']);

% loads good transects
good = load(['./all_transects_' adcp '.txt']);
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
dz = nanmean(diff(all_lmg.z));
nz = length(all_lmg.z);

% Loads OM mean currents (2004 - 2016)
load(['./omgeosvel_' adcp '.mat']);
om.uvmean = smean_om3.uvmean(smean_om3.z<=maxz, :, :); % lat x lon x depth
om.uvmean = om.uvmean.*exp(sqrt(-1)*theta*pi/180);% rotates the mean vectors
om.uvmean_err = smean_om3.err(smean_om3.z<=maxz, :, :);
om.lon_mean = smean_om3.lon;
om.lat_mean = smean_om3.lat;
om.z_mean = smean_om3.z(smean_om3.z<=maxz); 
om.uvmean = permute(om.uvmean, [2, 3, 1]);
om.uvmean_err = permute(om.uvmean_err, [2, 3, 1]);


clear smean_om3 sm3_lmg gadcp_lmg

% Loads bathymetry to determine which cruises are good for transport
b = load('./drake_mean_2021/data/all_drake/drakembf_avgd'); 
b.lat = [b.lats; b.lat; b.latn]; 
b.lon = [b.lonw b.lon]; 
b.bathyf = [b.bathyfw [b.bathyfs; b.bathyf; b.bathyfn]];

% Bathymetry boundaries
dpn = [-55, -54.70];
dps1 = [-0.00954, 5.9987, -1002.8768];
dps2 = [-0.00988, 6.2136, -1036.9253];

r = 0;
s = 0;



% Preallocates data
distlmg = nan(length([12.5:25:1200.5]), length(inds));
latlmg = distlmg;
lonlmg = distlmg;
ulmg = nan(length(all_lmg.z) + 1, length([12.5:25:1200.5]), length(inds));
ulmg_om = ulmg;
ulmg_om_mean = ulmg;
ulmg_raw = ulmg;
uprimeom = nan(length(26:8:210)+1, length([12.5:25:1200.5]), length(inds));
lmgdata.uprime_om_nb150 = uprimeom;
stheta = []; cs = []; Us = []; uoffs = []; nb = []; 

for i = 1:length(inds);
    
    clear data
    
    % Data for each transect
    data.uv = all_lmg.uv(:, I{inds(i)});
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
        & (data.lon > 293.5);
    data.lat = data.lat(dix);
    data.lon = data.lon(dix);
    data.uv = data.uv(:, dix);
    data.time = data.time(dix);
    data.bi = data.bi(dix);
    data.z = all_lmg.z;

    % Checks if data deviates from the main ship's direction
%     if flag(i) == 6;
    data = transectaverlmg(data);
    data.z = all_lmg.z;

        % If gaps are larger than some distance, make some fake data
%         data = fillgapslmg(data, 0.15, 'latitude');
    data.bi = interp2(b.lon, b.lat, b.bathyf, data.lon-360, data.lat);
%     end

    
    % If there gaps between 25 km < gap < 75 km,
    % fill with NaNs and linearly-interpolated latiudes/longitudes
    % for the objective mapping
            
    dist = [sw_dist(data.lat', data.lon', 'km')];
    indgap = find(dist>25);
            
    if ~isempty(indgap)
        data = fillgapslmg(data, gsize, 'distance');
        data.bi = interp2(b.lon, b.lat, b.bathyf, data.lon-360, data.lat);
    end
                
                
    % Interpolates mean to each depth bin and position
    data.uvmean = nan(size(data.uv)) + sqrt(-1)*nan(size(data.uv));           
            
    for zz  = 1:length(data.z);
        uv1 = squeeze(om.uvmean(:, :, zz));
        u = real(uv1);
        v = imag(uv1);
                
        umean = griddata(...
            om.lon_mean(~isnan(u)),...
            om.lat_mean(~isnan(u)),...
            u(~isnan(u)), data.lon(:)-360, data.lat(:), 'nearest');
        vmean = griddata(...
            om.lon_mean(~isnan(v)),...
            om.lat_mean(~isnan(v)),...
            v(~isnan(v)), data.lon(:)-360, data.lat(:), ...
            'nearest');
        data.uvmean(zz, :) = complex(umean, vmean);
    end
            
    % Calculates again alongtrack distance
    [data.dist data.ang] = sw_dist(data.lat, data.lon, 'km');
    data.ang = mod(data.ang, 360);
    data.ang = [data.ang(1); data.ang]';
    data.ang = nanmedian(data.ang);
            
    % Rotates vector to along/acrosstrack velocity components
    data.ur = imag(data.uv .*exp(-sqrt(-1)*...
        repmat(data.ang, size(data.z, 1), 1)*pi/180));
    data.ur_mean = imag(data.uvmean .*exp(-sqrt(-1)*...
        repmat(data.ang, size(data.z, 1), 1)*pi/180));

    % Ships velocity m/s
    data.uship = nanmstd(data.dist./diff(data.time))/86.4;

    % Calculates absolute geostrophic velocities from aviso sla + mdt
    if flag(i) == 2;
        sla = squeeze(nanmean(upd.sla(:, :, ...
            datenum(upd.dates)>= floor(data.time(1)) & ...
            datenum(upd.dates)<= ceil(data.time(end))), 3));
        uvadt = calcuv_sat(upd.lon-360, upd.lat, sla);
        [latsla lonsla] = meshgrid(upd.lat, upd.lon-360);
    end


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
            data.ur(:, indbi(1):end) = [];
            data.ur_mean(:, indbi(1):end) = [];
            data.lon(indbi(1):end) = [];
            data.lat(indbi(1):end) = [];
            data.bi(indbi(1):end) = [];
            data.time(indbi(1):end) = [];
            data.uv(:, indbi(1):end) = [];
%             data.dist(indbi(1):end) = [];
            data.uvmean(:, indbi(1):end) = [];
        end
%                 
     elseif (data.lat(1) < -60);% & (sum(diff(data.lat)<0)>0); %northbound
         if sum(diff(data.lat)<0)>0;
             data.dist(diff(data.lat) < 0) = ...
                 data.dist(diff(data.lat)<0)*-1;
         end
                
         % Positive towards the east*
         data.ur = data.ur*-1;
         data.ur_mean = data.ur_mean*-1;
         northbound(s, 1) = 1;
         nanbath = '1';

         % Removes data at the southern end that crosses the 500 m
         % isobath
         indbi = data.bi>=-500 & data.lat < -60;
         indbi = find(indbi==1);

         if ~isempty(indbi);
             data.ur(:, 1:indbi(end)) = [];
             data.ur_mean(:, 1:indbi(end)) = [];
             data.lon(1:indbi(end)) = [];
             data.lat(1:indbi(end)) = [];
             data.bi(1:indbi(end)) = [];
             data.time(1:indbi(end)) = [];
             data.uv(:, 1:indbi(end)) = [];
%              data.dist(1:indbi(end)) = [];
             data.uvmean(:, 1:indbi(end)) = [];
         end
    end

    % Cumulative distance
    [data.dist ~] = sw_dist(data.lat, data.lon, 'km');
    data.dist = [0; cumsum(data.dist)];

    % Bin averages every certain distance
    [bins, lataver] = line_bin(data.dist, ...
        repmat(data.lat', size(data.uv, 1), 1), gsize);
    [bins, lonaver] = line_bin(data.dist, ...
        repmat(data.lon', size(data.uv, 1), 1), gsize);
    [bins, uvaver] = line_bin(data.dist, data.ur, gsize);
    [bins, uvaver_mean] = line_bin(data.dist, data.ur_mean, gsize);
    [bins, biaver, bistd] = line_bin(data.dist, ...
        repmat(data.bi', size(data.uv, 1), 1), gsize, 1);%             [bins, distaver] = line_bin(data.dist, data.diust, gsize);

    lataver = lataver(1,:);
    lonaver = lonaver(1,:);
            
    % Removes mean
    urprime = uvaver - uvaver_mean;
            
    % Grid (depth x dist)
    [xx1, yy1] = meshgrid(bins(1):gsize:bins(end), data.z);
    indgd = ~isnan(urprime);
            
    % Objective mapping
    [urprimeOM Emap] = obj_map...
        (xx1(indgd)*1e3, yy1(indgd), urprime(indgd), ...
        decorr_x, decorr_y, error, xx1(:)*1e3, yy1(:));

    urprimeOM = reshape(urprimeOM, size(uvaver_mean));
            
    % Applies bathimetry mask where bottom depth + std dev is larger than
    % adcp depth bins.
    maskbat = floor(biaver + bistd) + repmat(data.z, 1, length(bins));
    if sum(maskbat(:)>0)>0
        urprimeOM(maskbat>0) = NaN;
    end

%     eval(['mask = uvaver(:, ' nanbath ')./uvaver(:, ' nanbath ');']);
%     eval(['urprimeOM(:, ' nanbath ') = urprimeOM(:, ' nanbath ').*mask;']);
    
    % Adds back the mean
    urfill = urprimeOM + uvaver_mean; 
    errmap  = reshape(Emap, size(uvaver_mean));
%     urfill = urprimeOM + uvaver_mean; %adds the mean back

    % Fills averaged data with objective mapping
    uvaver_fill = uvaver;
    uvaver_fill(isnan(uvaver_fill)) = urfill(isnan(uvaver_fill));

    % Applies mask
    if sum(maskbat(:)>0)>0
        uvaver_fill(maskbat>0) = NaN;
    end
            
    % Calculates error = 1 - sum([map - obs]^2)/var(obs)           
    varerr(1, s) = sqrt(nansum((urfill(:) - uvaver(:)).^2)/...
        sum(~isnan(uvaver(:))))/nanstd(uvaver(:));
        
    % Slab layer and saves data
%     ulmg(:, 1:length(bins), s) = [uvaver_fill(1, :); uvaver_fill];
%     ulmg_om(:, 1:length(bins), s) = [urfill(1, :); urfill];
%     ulmg_raw(:, 1:length(bins), s) = [uvaver(1, :); uvaver];
%     ulmg_om_mean(:, 1:length(bins), s) = [uvaver_mean(1, :); uvaver_mean];
%     distlmg(1:length(bins), s) = bins;
%     lonlmg(1:length(bins), s) = lonaver;
%     latlmg(1:length(bins), s) = lataver;
    lmgdata.time(1, s) = nanmean(data.time);
    lmgdata.flag(1, s) = flag(i);



    % Velocity offset (Firing et al. 2011) for os38nb
    if strcmp(adcp,'os38nb') == 1
        % Velocity offset (Firing et al. 2011)
        u = [(all_lmg.z(1)-dz/2)*uvaver_fill(1, :); ...
            uvaver_fill*dz];
        U = nansum(nansum(u)*gsize*1e3, 2);
	    zsp = ~isnan(u); 
        zsp = double(zsp).*([(all_lmg.z(1)-dz/2); ...
            dz*ones(nz,1)]*ones(1,size(u,2))); 
    	sp = sum(zsp(:)*gsize*1e3);
	            
        % Velocity and angle offset 
        uoffs = (Umean - U)/sp; 
% 	      uoffs(1, s) = uoff; 
  	    thetamis = asin(uoffs/data.uship);

        lmgdata.uoffs(1, s) = uoffs;
        lmgdata.thetaoff(1, s) = thetamis;

      
    end
end
               
% Saves data on a .mat             
% save(['lmgvel' adcp '_' num2str(decorr_x/1000) '_' num2str(decorr_y)], ...
%     'lmgdata');
    
% Calculates mean by leaving outliers out
aa = find(abs(abs(lmgdata.thetaoff)-nanmean(abs(lmgdata.thetaoff)))<...
    3*nanstd(abs(lmgdata.thetaoff)));

% Theta missalignment
thetamiss = nanmean(abs(lmgdata.thetaoff(aa)));

% Uoffset
% southbound transects
bb = find(abs(lmgdata.uoffs(northbound == 0) - ...
    nanmean(lmgdata.uoffs(northbound == 0))) < ...
    3*nanstd(lmgdata.uoffs(northbound == 0)));
uoffs_south = lmgdata.uoffs(northbound == 0);
uoffs_south = nanmean(uoffs_south(bb));

% northbound transects
dd = find(abs(lmgdata.uoffs(northbound == 1) - ...
    nanmean(lmgdata.uoffs(northbound == 1))) < ...
    3*nanstd(lmgdata.uoffs(northbound == 1)));
uoffs_north = lmgdata.uoffs(northbound == 1);
uoffs_north = nanmean(uoffs_north(dd));


fprintf('\nTransducer misalignment angle = %0.4f ^o\n', thetamiss*180/pi)
save('trandsmiss_theta.mat', 'thetamiss', 'uoffs_south', 'uoffs_north');
fprintf('Number of transects used: %s\n', num2str(s));


figure('color','w');
ax = tight_subplot(1, 1, [0.01 0.01], [0.10 0.01], [0.15 0.02]);
axes(ax(1));
t1 = plot([lmgdata.time-datenum(1996,1,1)]/365.24, ...
    abs(lmgdata.thetaoff*180/pi), '.-b', 'linewidth',1.5);
hold on
plot([lmgdata.time(aa)-datenum(1996,1,1)]/365.24, ...
    abs(lmgdata.thetaoff(aa)*180/pi), 'sr', 'linewidth',1.5);
ylabel('\theta_m_i_s [^o]');
ax(1).FontSize = 12;
ax(1).YTick = [0:0.1:1]';
ylim([-0.05 1])
xlabel('Years since Jan 1996')
text(16, 0.8, ['\theta_m_i_s_s = ' num2str(thetamiss*180/pi, '%0.4f') ...
    '^o'])'

figure('color', 'w');
[~, ~, s1] = histf((lmgdata.uoffs(northbound==0)), [-0.1:0.005:0.1]);
hold on
[~, ~, s2] = histf((lmgdata.uoffs(northbound==1)), [-0.1:0.005:0.1]);
s1.FaceAlpha = 0.5;
s2.FaceAlpha = 0.5;
set(gca, 'TickLength', [0.03 0.03]);
legend([s1; s2], 'Southbound', 'Northbound', 'Location', 'NorthEast');
text(-0.08, 27, ...
    ['Southbound Mean |\Deltau| = ' num2str(uoffs_south, '%0.4f') ' m s^-^1']);
text(-0.08, 25, ...
    ['Northbound Mean |\Deltau| = ' num2str(uoffs_north, '%0.4f') ' m s^-^1']);
xlabel('|\Deltau| [m s^-^1]')
ylabel('Frequency');
ylim([0 30]);


