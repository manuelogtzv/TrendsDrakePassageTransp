% Program to estimate time series of transport per pair of streamlines. For
% each unique pair of streamlines (Gutierrez-Villanueva et al. 2020),
% velocity estimates are binned and integrated in depth and alongtrack
% distance to estimate the transport [m^3 s^{-1}] per pair of streamlines
% per transect. The total velocity profiles are interpolated from 25 km to
% 11 km alongtrack resolution to match that of the geostrophic and
% reference velocity.
%
% 2023/06/28 - Fixes an issue with selection of the transects.
% More consistent naming of variables.


clear all; 
% close all;


lat0 = -55;
lon0 = -65;
omega = 2*pi/86400;
% gsize = 25000;% Grid-box size in meters
overlapping = 0.0;% Overlapping 0.5 = 50%
depths_psi = [100 200];% Averages mean geostrophic velocities and psi over these two depths
psi_saf = -0.52;
psi_pf = -1.25;
psi_saccf = -1.6;
nday = 20;
gsize = 25;
ngsize = 12;
alpha = 0.05;
op_season = 1; % removes seasonal cycle
op_vid = 0;
op_mostreptran = 1;
dPsi = 0.10;
adcp = 'geos';

% Loads aviso SSH daily maps for Drake Passage
load('../BigDatasets/SeaLevelAnomalies.mat');
ssh = sla;
clear sla;


% Sampling period
din2 = datenum(2005, 10, 1);
dfin1 = datenum(2019, 4, 30);



% Mean streamwise position of the ACC fronts
saf_psi = [-0.40 -0.60];%Subantarctic Front
pf_psi = [-1.00 -1.30];%Polar Front
saccf_psi = [-1.55 -1.65];%Southern ACC Front




if strcmp(adcp, 'os38nb');
    % Loads data
    load(['./lmgvelos38nb_50_300.mat']);
    % Makes the transect that start from south to north to now be north to
    % south
    indsouth = find(lmgdata.northbound == 1);
    indnorth = find(lmgdata.northbound == 0);

    % Adds the uoffs
    load('trandsmiss_theta.mat');
    uoffs = mean(abs([uoffs_north, uoffs_south]));
    lmgdata.u(:, :, indsouth) = lmgdata.u(:, :, indsouth) - uoffs;
    lmgdata.u(:, :, indnorth) = lmgdata.u(:, :, indnorth) + uoffs;
    nomz = 760;
    lmgdata.u = lmgdata.u(lmgdata.z<=nomz, :, :);
    lmgdata.u_om = lmgdata.u_om(lmgdata.z<=nomz, :, :);
    lmgdata.z = lmgdata.z(lmgdata.z<=nomz);

    % Makes all transect to start from north
    for i = 1:length(indsouth);
        [lmgdata.dist(:, indsouth(i)) ii] = ...
            sort(abs(nanmax(lmgdata.dist(:, indsouth(i))) - ...
            lmgdata.dist(:, indsouth(i))) + 12.5);
        lmgdata.u(:, :, indsouth(i)) = lmgdata.u(:, ii, indsouth(i));
%         lmgdata.u_om(:, :, indsouth(i)) = lmgdata.u_om(:, ii, indsouth(i));
        lmgdata.lon(:, indsouth(i)) = lmgdata.lon(ii, indsouth(i));
        lmgdata.lat(:, indsouth(i)) = lmgdata.lat(ii, indsouth(i));
    end

%     dist = [lmgdata.dist(1):gsize:size(lmgdata.dist,1)*gsize ];
    
    % Renames data
    data.u = lmgdata.u(:, :, lmgdata.time>=din2 & lmgdata.time<=dfin1);
    data.lon = lmgdata.lon(:, lmgdata.time>=din2 & lmgdata.time<=dfin1) - 360;
    data.lat = lmgdata.lat(:, lmgdata.time>=din2 & lmgdata.time<=dfin1);
    data.dist = lmgdata.dist(:, lmgdata.time>=din2 & lmgdata.time<=dfin1);
    data.z = lmgdata.z;
    data.time = lmgdata.time(1, lmgdata.time>=din2 & lmgdata.time<=dfin1);;
 %     data.distdiff = ones(size(data.dist))*25;
   
    clear lmgdata uoffs* indsouth indnorth 

    S_unique = -[0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90...
        1.0 1.15 1.25 1.35 1.45 1.50 1.55 1.60 1.65 1.70];
    
    clear totaltransp;


elseif strcmp(adcp, 'geos');
    load(['geostrophic_noaver.mat']);

%     data = geosbaroc;
    data.time = geosbaroc.time(geosbaroc.time>=din2 & geosbaroc.time<=dfin1);
    data.u = geosbaroc.u(:, :, geosbaroc.time>=din2 & geosbaroc.time<=dfin1);
    data.lon = geosbaroc.lon_u(:, geosbaroc.time>=din2 & geosbaroc.time<=dfin1);
    data.lat = geosbaroc.lat_u(:, geosbaroc.time>=din2 & geosbaroc.time<=dfin1);
    data.dist = geosbaroc.dist_u(:, geosbaroc.time>=din2 & geosbaroc.time<=dfin1);
    data.distdx = diff(geosbaroc.dist(:, geosbaroc.time>=din2 & geosbaroc.time<=dfin1));
    data.z = geosbaroc.z;
    % geosbci.distu = distubci; 
    clear geosbaroc 
    S_unique = -[0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90...
        1.0 1.15 1.25 1.35 1.45 1.50 1.55 1.60 1.65 1.70];
else
    load(['reference_os38nb_760.mat' ]);

    data.time = geosref.time(geosref.time>=din2 & geosref.time<=dfin1)
    data.u = geosref.u(:, :, geosref.time>=din2 & geosref.time<=dfin1);
    data.lon = geosref.lon(:, geosref.time>=din2 & geosref.time<=dfin1);
    data.lat = geosref.lat(:, geosref.time>=din2 & geosref.time<=dfin1);
    data.dist = geosref.dist(:, geosref.time>=din2 & geosref.time<=dfin1);
    data.distdx = [zeros(size(data.dist(1, :))); data.dist];
    data.distdx = diff(data.distdx);
    data.ek = geosref.ektransp(:, geosref.time>=din2 & geosref.time<=dfin1);

    S_unique = -[0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90...
        1.0 1.15 1.25 1.35 1.45 1.50 1.55 1.60 1.65 1.70];
end


% Option to make per frontal region
% S_unique = [-0.30 saf_psi pf_psi saccf_psi];
S_m = 0.5*(S_unique(1:end-1) + S_unique(2:end));

% % Finds coincident xbt and adcp transects
% % os38nb 
% [x y] = meshgrid(timegeo_760, timeos38_970);
% [t1 t2] = meshgrid(1:length(timegeo_760), 1:length(timeos38_970));
% ind_os38_a = t2(abs(x-y)<=4);


% Option for using the most repeated transects
if op_mostreptran == 1;
    fprintf('\n\n total number of transects: %s\n',...
        num2str(length(data.time)));

    box1 = [ -64.606712283408200 -56.324379209501110;
             -64.534099360692622 -56.568973205588243;
             -64.373586584163448 -56.541796094911895;
             -64.438556041330017 -56.351556320177458];

    box2 = [-63.009254066180596 -61.670918367346943;
            -62.908300616937744 -61.915816326530610;
            -62.504486819966345 -61.752551020408163;
            -62.625630959057766 -61.426020408163268];

    ii1 = inpolygon(data.lon, data.lat, box1(:, 1), ...
        box1(:, 2));

    [~, yj] = find(ii1 == 1);
    yj = unique(yj);

    mm1 = inpolygon(data.lon(:, yj), data.lat(:, yj), ...
        box2(:, 1), box2(:, 2));

    [~, zj] = find(mm1 == 1);
    zj = unique(zj);

    % HOw many transects along the most repeated line
    fprintf('\n\n total number of transects (most repeated): %s\n',...
        num2str(length(data.time(yj(zj)))));
end

% dz = diff(data.z(2:3));

% Averages time-mean streamfunction, velocity over specific depth ranges
load('omgeosvel_nb150.mat'); %mean geostrophic vel
ind_depths_psi = find(smean_om3.z >= depths_psi(1) & ...
    smean_om3.z <= depths_psi(2));

psi_mean = ...
    squeeze(nanmstd(squeeze(smean_om3.psi(ind_depths_psi, :, :))));%Average

% Combines MDT with \Psi
psiMean = psiDPmdt('../Datasets/', smean_om3.lat, smean_om3.lon, psi_mean);


clear r data_xbt* data_lmg*

fig22 = figure('color','w');


r = 0;

for i = 1:length(data.time);% Loops on aviso time
    ind = find(abs(floor(data.time(i)) - ssh.time) <= 0);
    
    if ~isempty(ind);% finds coincident xbt and aviso times

        r = i;

        % Gradients of SSH + \Psi = \Psi^\ast
        sla = nanmstd(ssh.sla(:, :, ind-nday:ind+nday), 3);
        sla = interp2(double(ssh.lon)', double(ssh.lat)', sla',...
            psiMean.lon, psiMean.lat);

        maps(:, :, r) = sla + psiMean.psi;
%             maps(:, :, r) = psiMean.psi;
                    
        
        %%%%%%%%%%%%% FIX CODE HERE FOR TOTAL TRANSPORT ONLY %%%%%%%%%%%%%
        if strcmp(adcp, 'os38nb');
            % Interpolates total velocity to a 11 km distance grid
            ndist = [12.5:ngsize:nanmax(data.dist(~isnan(data.lat(:, i)), i))+ngsize/2];
            
            % Interpolate coordinates
            nlat = interp1(data.dist(~isnan(data.lat(:, i)), i), ...
                data.lat(~isnan(data.lat(:, i)), i), ndist, 'linear');
            nlon = interp1(data.lat(~isnan(data.lat(:, i)), i), ...
                data.lon(~isnan(data.lat(:, i)), i), nlat, 'linear');

            [z dist] = meshgrid(data.z, data.dist(~isnan(data.lat(:, i)), i));
            [nz ndist] = meshgrid(data.z, ndist);

            % Interpolates data to finer horizontal grid
            nu = griddata(dist, z, ...
                squeeze(data.u(:, ~isnan(data.lat(:, i)), i))', ...
                ndist, nz, 'linear')';
            ndx = ndist(~isnan(nlat));

        else
            % For geostrophic and reference, just rename and remove NaNs
            nlat = data.lat(~isnan(data.lat(:, i)), i);
            nlon = data.lon(~isnan(data.lon(:, i)), i);
            nu = squeeze(data.u(:, ~isnan(data.lat(:, i)), i));
            ndx = data.distdx(~isnan(data.lat(:, i)), i);
        end


               
        
        figure(fig22);
        clf;
        [cc hh] = contour(smean_om3.lon, smean_om3.lat, ...
            squeeze(maps(:, :, r)), S_unique, 'linewidth',2);
        clabel(cc, hh, 'labelspacing', 300, 'fontsize', 16);
        hold on; 
        ax = axis; caxis([-2.0 -0.3])
        colormap(viridis(50));
        plot(data.lon(~isnan(data.lon(:, i)), i), ...
            data.lat(~isnan(data.lon(:, i)), i), '.k');
        plot(nlon, nlat, 'xr');
        set(gca, 'fontsize', 16, 'xtick', [-3:0]*1e5, 'xticklabel',...
            num2str([-3:0]'*1e2), 'ytick', [-8:0]*1e5, ...
            'yticklabel', num2str([-8:0]'*1e2));
%         xlabel(' Longitude ', 'interpreter', 'latex', 'fontsize', 18);
%         ylabel(' km ', 'interpreter', 'latex', 'fontsize', 18);
        title(datestr(double(ssh.time(ind))), 'fontsize', 18, ...
            'Interpreter',  'latex');
        axis equal;
            
        % Gets position of each contour level
        S = contourdata(cc);
        cont_le = length(S);%count contour levels

        % Saves contour x,y and level in other variables
        m = 0;
        data_nan = ones(size(nlon));

        % Removes contours that are closed
        for h = 1:cont_le;
            if S(h).isopen == 1;
                m = m + 1;
                S_level(r, m) = S(h).level;
                S_x{r, m} = S(h).xdata;
                S_y{r, m} = S(h).ydata;
                S_curv(r, m) = sqrt((S_x{r, m}(1) - S_x{r, m}(end))^2 + ...
                    (S_y{r, m}(1) - S_y{r, m}(end))^2);
                
                else S(h).isopen == 0;% If closed contour is found
                
                % Finds data within closed contours
                idata = inpolygon(nlon(:), nlat(:),...
                    S(h).xdata, S(h).ydata);
                
                if ~isempty(idata);
                    % Assigns NaN if data is found within closed
                    % contours
                    data_nan(idata) = NaN;
                end
            end
        end
        
        for j = 1:length(S_unique);%number of unique contours
            ind_coinc = find(S_level(r,:) == S_unique(j));
            
            for m = 1:length(ind_coinc);%length of each contour
                len(m,1) = length(S_x{r,ind_coinc(m)});
            end
            
            % Find max length
            ind_max = find(len == nanmax(len));
            ind_coinc = ind_coinc(ind_max);%keeps the longest contour
            
            % Saves contours
            cont_x{r,j} = S_x{r,ind_coinc};
            cont_y{r,j} = S_y{r,ind_coinc};
            
            clear ind_max ind_coinc len
        end
           
        clear m j 
        
        % Now uses polygons to find data inside them
        for b = 1:length(S_m);
            pol_x = [cont_x{r,b} ; flipud(cont_x{r,b+1})];
            pol_y = [cont_y{r,b} ; flipud(cont_y{r,b+1})];
            
            in_lmg = inpolygon(nlon(:), nlat(:), pol_x, pol_y);
            
            % Velocity data
            u_lmg{i, b} = nu(:, in_lmg);
            lon_lmg{i, b} = nlon(in_lmg);
            lat_lmg{i, b} = nlat(in_lmg);
            dx_lmg{i, b} = ndx(in_lmg);
            num(i, b) = sum(in_lmg);
            
            if strcmp(adcp, 'os38nb');
                dz = 24; % meters
                % transport per streamline
                ulmg_sum = nansum(u_lmg{i, b}(2:end, :), 1)*dz*ngsize*1e3;

                % Slab layer
                ulmg_sum = ulmg_sum + u_lmg{i, b}(1, :)*ngsize*(data.z(2)-dz/2)*1e3; 

                transp_lmg(i, b) = sum(ulmg_sum(:)); % sum

            elseif strcmp(adcp, 'geos');
                dz = 10;
                % transport per streamline
%                 ulmg_sum = nansum(u_lmg{i, b}*dz, 1).*dx_lmg{i, b}*1e3;
%                 transp_lmg(i, b) = sum(ulmg_sum(:));
                ulmg_sum = nansum(u_lmg{i, b}, 1)*dz*dx_lmg{i, b}*1e3;

                transp_lmg(i, b) = sum(ulmg_sum(:)); % sum

            else
                dldp = 648;
                
                % transport per streamline
                ulmg_sum = u_lmg{i, b}*dldp*dx_lmg{i, b}*1e3;

                transp_lmg(i, b) = sum(ulmg_sum(:)); % sum

                transp_ek(i, b) = sum(data.ek(in_lmg, i));
            end
        end
    end
end

% Nans zeros
transp_lmg(transp_lmg == 0) = NaN;

[nPsi ntime] = meshgrid(S_m, data.time);
anom = transp_lmg;%-nanmean(transp_lmg); % Calculates anomalies

% Saves index for transects along the most repeated line
indx = yj(zj);

% Saves data
stream.transp = transp_lmg;
stream.lon = lon_lmg;
stream.lat = lat_lmg;
stream.num = num;
stream.time = data.time;
stream.psi = S_unique;
stream.psic = S_m;
stream.adcp = adcp;
stream.indx_rep = yj(zj);
stream.dx = ngsize;
stream.nday = nday;
if strcmp(adcp, 'ref'); % Includes ekman transport if reference is used
    stream.ek = transp_ek;
end
stream.doc = {'STREAMWISE BINNING OF VELOCITY AND TRASNPORT';...
              'adcp is the adcp used';...
              'transp is the transport per pair of streamlines [m^3 s^{-1}]';...
              'lon and lat are the longitude and latitude found per pair of streamlines per transect';...
              'num is the number of profiles found per pair streamlines per transect';...
              'psi is the streamfunction pairs [m]';...
              'psic is the centered streamfunction values';...
              'time is the time vector [days]';...
              'indx_rep lists the indexes of the transects along the most repeated line';...
              'dx is the along-track resolution of the data';...
              'nday is the half the number of days that time series of sea level anomalies were averaged prior the binning'};
save(['Streamwise_binning_' adcp '.mat'], 'stream');





