% Program to load different ADCP transects from the LMG starting from 2004
% till 2019. Program seperates data into transects and checks
% for the conditions below:
%
% a) Separates into southbound and northbound transects
% b) Checks that each transects is completed in <4 days
% c) Checks that transects are reach the 500 m isobath at both ends (south
% america and antarctica)
% d) Removes transects with gaps larger than 150 km
% 
% Fills gaps <150 km in the alongtrack distance with an objective mapping.
% ONLY the across-transect velocity is employed.
%
% Manuel O. Gutierrez V.
% 2021/08/23
% 2021/08/31 - Includes a list of manually neglected transects.
% 2021/11/29 - Includes adcp name
% 2023/06/12 - Uses 5 min ping data


clear all; 
close all;

pycmd = '/Users/manuelgutierrez/anaconda3/bin/python';

gsize = 25; %km
cont = [-1.0:0.1:1.0];%vel contours
contmaperr = [0.01 0.05 0.10 0.20 0.5 0.80 1.0];%mapping error contours
decorr_x = 50000; %meters distance
decorr_y = 300; %meters depth
error = 0.2;
cmap = getPyPlot_cMap('seismic', length(cont)-1, [], pycmd);
cmaperr = flipud(getPyPlot_cMap('pink', length(contmaperr), [], pycmd));
cmaperrvar = flipud(getPyPlot_cMap('YlGnBu', length(contmaperr), [], pycmd));
boxlim = [0.0125    0.8625   -1.0300   -0.0460]*1e3;
adcp = 'os38nb'; % os38nb or nb150
plot_fig = 0;

% loads transects
load(['./drake_mean_2021/data/all_drake/lmg_' adcp '_5min_25pg.mat']);

% loads bad transects
bad = datenum(load('./drake_mean_2021/badtransects'));


% Splits into transects
I = interval(all_lmg.t, 2, 'n');
fprintf('\n\nNumber of transects available: %3.0f \n\n', length(I));
fprintf('Starting date: %s \n', datestr(all_lmg.t(1)));
fprintf('Ending date: %s \n\n', datestr(all_lmg.t(end)));


% Removes tides
all_lmg.uv = all_lmg.uv - repmat(complex(tide.u, tide.v)', ...
    size(all_lmg.uv,1), 1);

% Removes last 8 bins of data
all_lmg.uv = all_lmg.uv(1:end-8,:);
all_lmg.z = all_lmg.z(1:end-8);

% Loads OM mean currents (2004 - 2016)
load('../mean_eddy_sADCP_om3.mat');
om.uvmean = smean_om3.uvmean; % lat x lon x depth
om.uvmean_err = smean_om3.err;
om.lon_mean = smean_om3.lon;
om.lat_mean = smean_om3.lat;
om.z_mean = smean_om3.z; 

clear smean_om3

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

fig1 = figure('color','w');
ha = tight_subplot(2, 3, [0.2 0.11],[0.20 0.05], [0.20 0.05]);

% Preallocates data
distlmg = nan(length([12.5:25:1200.5]), length(I));
latlmg = distlmg;
lonlmg = distlmg;
ulmg = nan(length(all_lmg.z)+1, length([12.5:25:1200.5]), length(I));
ulmg_om = ulmg;


for i = 1:length(I);
    
    clear data
    
    % Data for each transect
    data.uv = all_lmg.uv(:, I{i});
    data.time = all_lmg.t(I{i});
    data.lon = all_lmg.lon(I{i});
    data.lat = all_lmg.lat(I{i});
    data.bi = interp2(b.lon, b.lat, b.bathyf, data.lon-360, data.lat);
    
    % If transect is not a bad one
    if isempty(find(abs(bad - nanmean(data.time)) < 1));
    
    % Calculates maximum distance (i.e. between minimum and maximum
    % transect)
    maxdist = sw_dist([nanmax(data.lat) nanmin(data.lat)], ...
        [data.lon(data.lat == nanmax(data.lat)) ...
        data.lon(data.lat == nanmin(data.lat))], 'km');

    % Sorts bottom data
    [latsort ii] = sort(data.lat);
    bisort = data.bi(ii);
    
    % Checks if data is within some minimum bottom depth and if the
    % transect length is equal or larger than 700 km
    if (bisort(1) >= -1000) & (bisort(end) >= -1000) & (maxdist >= 700);
        
        % Checks if data has a gap larger than 150 km.
        data.dist = sw_dist(latsort, data.lon(ii), 'km');
        dxdist = diff(data.dist);
        
        if isempty(find(abs(dxdist) >= 75))
            
            s = s + 1;

            % Removes data away from boundaries
            dix = (data.lat < dpn(1)) & ...
                (data.lat > polyval(dps2, data.lon))...
                & (data.lon > 293.5);
            data.lat = data.lat(dix);
            data.lon = data.lon(dix);
            data.uv = data.uv(:, dix);
            data.time = data.time(dix);
            data.z = all_lmg.z;
            
            if floor(nanmean(data.time)) == datenum(2007,11,06);
                data.lat(47) = [];
                data.lon(47) = [];
                data.uv(:, 47) = [];
                data.time(47) = [];
            end
                      
            % If there gaps between 25 km < gap < 75 km,
            % fill with NaNs and linearly-interpolated latiudes/longitudes
            % for the objective mapping
            
            dist = [sw_dist(data.lat', data.lon', 'km')];
            indgap = find(dist>25);
            
            if ~isempty(indgap)
                
                if ~isempty(diff(data.lat)==0)
                    data.lat(diff(data.lat) == 0) = ...
                        data.lat(diff(data.lat) == 0) + 0.0001;
                end

                [v, w] = unique(data.lat, 'stable' );
                dupind = setdiff( 1:numel(data.lat), w);
                if ~isempty(dupind);
                    data.lat(dupind(1)) = data.lat(dupind(1)) + 0.001;
                end

                data = fillgapslmg(data, indgap);
                data.bi = interp2(b.lon, b.lat, b.bathyf,...
                    data.lon-360, data.lat);
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
                    u(~isnan(u)), data.lon(:)-360, data.lat(:), ...
                    'nearest');
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
            
            % Rotates vector to along/acrosstrack velocity components
            data.ur = imag(data.uv .*exp(-sqrt(-1)*...
                repmat(data.ang, size(data.z, 1), 1)*pi/180));
            data.ur_mean = imag(data.uvmean .*exp(-sqrt(-1)*...
                repmat(data.ang, size(data.z, 1), 1)*pi/180));
            
            
            % Checks if it transect is south or northbound, and corrects
            % for the reversing orientation
            if (data.lat(1) > -60);% & (sum(diff(data.lat)>0)>0); %southbound
                if sum(diff(data.lat)>0)>0;
                    data.dist(diff(data.lat) > 0) = ...
                        data.dist(diff(data.lat)>0)*-1;
                end
                northbound(s, 1) = 0;
                nanbath = 'end';
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
            end
            
            
            % Cumulative distance
            data.dist = [0; cumsum(data.dist)];
            
            % Bin averages every certain distance
            [bins, lataver] = line_bin(data.dist, ...
                repmat(data.lat', size(data.uv, 1), 1), gsize);
            [bins, lonaver] = line_bin(data.dist, ...
                repmat(data.lon', size(data.uv, 1), 1), gsize);
            [bins, uvaver] = line_bin(data.dist, data.ur, gsize);
            [bins, uvaver_mean] = line_bin(data.dist, data.ur_mean, gsize);
            [bins, biaver] = line_bin(data.dist, ...
                repmat(data.bi', size(data.uv, 1), 1), gsize);%             [bins, distaver] = line_bin(data.dist, data.diust, gsize);

            lataver = lataver(1,:);
            lonaver = lonaver(1,:);
            
            % Removes mean
            urprime = uvaver - uvaver_mean;
            
             % New grid (depth x dist)
            [xx1, yy1] = meshgrid(bins(1):gsize:bins(end), data.z);

            indgd = ~isnan(urprime);
            
            % Objective mapping
            [urprimeOM Emap] = obj_map...
                    (xx1(indgd)*1e3, yy1(indgd), urprime(indgd), ...
                    decorr_x, decorr_y, error, xx1(:)*1e3, yy1(:));

            % Adds back the mean
            urprimeOM = reshape(urprimeOM, size(uvaver_mean));
            
            % Applies bathimetry mask
            eval(['mask = uvaver(:, ' nanbath ')./uvaver(:, ' nanbath ');']);
            eval(['urprimeOM(:, ' nanbath ') = urprimeOM(:, ' nanbath ').*mask;']);

            
            urfill = urprimeOM + uvaver_mean; 
            errmap  = reshape(Emap, size(uvaver_mean));
%             urfill = urprimeOM + uvaver_mean; %adds the mean back

            % Fills averaged data with objective mapping
            uvaver_fill = uvaver;
            uvaver_fill(isnan(uvaver_fill)) = urfill(isnan(uvaver_fill));
            
            
            
            % Calculates error = 1 - sum([map - obs]^2)/var(obs)           
            varerr(1, s) = sqrt(nansum((urfill(:) - uvaver(:)).^2)/...
                sum(~isnan(uvaver(:))))/nanstd(uvaver(:));

            % Plotting
            if plot_fig == 1;
            if floor(nanmean(data.time)) > datenum(2007,03,24) & ...
                    floor(nanmean(data.time)) < datenum(2007,04,05) || ...
                    floor(nanmean(data.time)) == datenum(2006,03,25) || ...
                    floor(nanmean(data.time)) == datenum(2007,09,27) || ...
                    floor(nanmean(data.time)) == datenum(2009,09,21) || ...
                    floor(nanmean(data.time)) == datenum(2011,05,07) || ...
                    floor(nanmean(data.time)) == datenum(2010,11,15) || ...
                    floor(nanmean(data.time)) == datenum(2015,06,28);
            figure(fig1);
            
            axes(ha(1));
            as = quiver(data.lon-360, data.lat, real(data.uv(2, :))'*1.15, ...
                imag(data.uv(2, :))'*1.15, 0, 'color', 'k');
            hold on
            quiver(om.lon_mean, om.lat_mean, ...
                squeeze(real(om.uvmean(:, :, 2)))*0.85,...
                squeeze(imag(om.uvmean(:, :, 2)))*0.85, 0, 'color', 'r');
            d1 = plot(data.lon(1)-360, data.lat(1), '^k', 'markersize', 10,...
                'markerfacecolor', 'k');
            d2 = plot(data.lon(end)-360, data.lat(end), 'sr', ...
                'markersize', 10, 'markerfacecolor', 'k');
            axis equal; axis([-67 -56 -65 -53]);
            title(datestr(nanmean(data.time), 2));
            legend([d1; d2], 'ti', 'tf');
            
            
            axes(ha(2));
            cla;
            contourf(bins, -data.z, uvaver, cont, 'edgecolor', 'none');
            hold on
            contour(bins, -data.z, uvaver, [0 0], 'color', 'k', ...
                'linewidth', 1.5);
            bat = plot(bins, biaver, 'k');
            caxis(ha(2), [cont(1) cont(end)]);
            colormap(ha(2), cmap);
            caxis(ha(2), [cont(1) cont(end)]);
            colormap(ha(2), cmap);
            title(ha(2), 'U_{r} alongtrack-averaged 25 km');
            ylabel(ha(2), 'Depth [m]');
            xlabel(ha(2), 'Alongtrack distance [km]');
            axis(ha(2), boxlim);
            
            
            axes(ha(3));
            cla;
            contourf(bins, -data.z, urfill, cont, 'edgecolor', 'none');
            hold on
            contour(bins, -data.z, urfill, [0 0], 'color', 'k', ...
                'linewidth', 1.5);
            caxis(ha(3), [cont(1) cont(end)]);
            colormap(ha(3), cmap);
            caxis(ha(3), [cont(1) cont(end)]);
            colormap(ha(3), cmap);
            col3 = colorbar(ha(3), 'Location', 'SouthOutside');
            set(col3, 'Position', [ha(3).Position(1) 0.54 ...
                ha(3).Position(3) 0.02], 'Ticks', cont(1:2:end),...
                'TickLabels', num2str(cont(1:2:end)'));
            title(col3, '[m s^-^1]');
            title(ha(3), 'Obj. Mapped U_{r}');
            ylabel(ha(3), 'Depth [m]');
            xlabel(ha(3), 'Alongtrack distance [km]');
            axis(ha(3), boxlim);

                        
            axes(ha(4));
            cla;
            contourf(bins, -data.z, log10(errmap), log10(contmaperr),...
                'edgecolor', 'none');
            hold on
            colormap(ha(4), cmaperr);
            caxis(ha(4), log10([contmaperr(1) contmaperr(end)]));
            colormap(ha(4), cmaperr);
            col4 = colorbar(ha(4), 'Location', 'South');
            set(col4, 'Position', [ha(4).Position(1) 0.07 ...
                ha(4).Position(3) 0.02], 'Ticks', log10(contmaperr),...
                'TickLabels', num2str(contmaperr'));
            title(col4, 'Normalized mapping error');
            xlabel(ha(4), 'Alongtrack distance [km]');
            ylabel(ha(4), 'Depth [m]');
            axis(ha(4), boxlim);

            axes(ha(5));
            cla;
            contourf(bins, -data.z, uvaver_fill, cont,...
                'edgecolor', 'none');
            hold on
            contour(bins, -data.z, uvaver_fill, [0 0], 'color', 'k', ...
                'linewidth', 1.5);
            caxis(ha(5), [cont(1) cont(end)]);
            colormap(ha(5), cmap);
            caxis(ha(5), [cont(1) cont(end)]);
            colormap(ha(5), cmap);
            title(ha(5), 'U_r filled with Obj. Mapped U_r ');
            xlabel(ha(5), 'Alongtrack distance [km]');
            ylabel(ha(5), 'Depth [m]');
            axis(ha(5), boxlim);
            
            axes(ha(6));
            cla;
            contourf(bins, -data.z, uvaver - urfill, cont,...
                'edgecolor', 'none');
            hold on
            caxis(ha(6), [cont(1) cont(end)]);
            colormap(ha(6), cmap);
            caxis(ha(6), [cont(1) cont(end)]);
            colormap(ha(6), cmap);
            title(ha(6), 'U_r - Obj. Mapped U_r ');
            xlabel(ha(6), 'Alongtrack distance [km]');
            ylabel(ha(6), 'Depth [m]');
            axis(ha(6), boxlim);
            
%             pause;%(0.1)
            print('-depsc', ['./Figures/transect' adcp '_' ...
                datestr(floor(nanmean(data.time)), 'yyyymmdd')]);
            delete([as; d1; d2; col3; col4; bat])
            end
            end
            
            
            % Slab layer and saves data
            ulmg(:, 1:length(bins), s) = [uvaver_fill(1, :); uvaver_fill];
            ulmg_om(:, 1:length(bins), s) = [urfill(1, :); urfill];
            distlmg(1:length(bins), s) = bins;
            lonlmg(1:length(bins), s) = lonaver;
            latlmg(1:length(bins), s) = lataver;
            lmgdata.time(1, s) = nanmean(data.time);
            
            
            
        else
            fprintf('Transect discarded: %s \n', ...
                datestr(nanmean(data.time)));
            r = r + 1;
        end
        
    else
        fprintf('Transect discarded: %s \n', datestr(nanmean(data.time)));
        r = r + 1;
    end
    
    else
        fprintf('Transect discarded: %s \n', datestr(nanmean(data.time)));
        r = r + 1;
    end
end





lmgdata.u = ulmg(:, :, 1:s);
lmgdata.u_om = ulmg_om(:, :, 1:s);
lmgdata.dist = distlmg(:, 1:s);
lmgdata.lon = lonlmg(:, 1:s);;
lmgdata.lat = latlmg(:, 1:s);;
lmgdata.northbound = northbound;
lmgdata.z = [0; data.z];
lmgdata.rms_map = varerr;
lmgdata.decorrx = decorr_x;
lmgdata.decorrz = decorr_y;
lmgdata.res = gsize;
lmgdata.n2s = error;
lmgdata.adcp = adcp;
lmgdata.doc = {'u is across-transect velocity [m/s] (positive to east); gaps filled with OM';
                'u_om is the objectively mapped velocity [m/s]';
                'dist is along track distance [km]';
                'lon and lat are coordinates averaged for each 25 km bin';
                'northbound = 1 is a transect going from south to north';
                'z is depth [m]';
                'rms map is the normalized RMS between obj. map. and data';
                'decorrx and decorrz are the decorrlation scales [m]';
                'n2s is the mapping signal-to-noise ratio';
                'res is the alongtrack resolution [km]'};

% Remove data at the end of transects
ind = find(floor(lmgdata.time) == datenum(2014, 11, 20));

lmgdata.lat(lmgdata.lat(:, ind)<-61.99, ind) = NaN;
lmgdata.lon(lmgdata.lat(:, ind)<-61.99, ind) = NaN;
lmgdata.dist(lmgdata.lat(:, ind)<-61.99, ind) = NaN;
% lmgdata.uprime_om_nb150(:, lmgdata.lat<-61.99, ind) = NaN;
lmgdata.u(:, lmgdata.lat(:, ind)<-61.99, ind) = NaN;
lmgdata.u_om(:, lmgdata.lat(:, ind)<-61.99, ind) = NaN;
% lmgdata.u_aver(:, lmgdata.lat(:, ind)<-61.99, ind) = NaN;
% lmgdata.u_mean(:, lmgdata.lat(:, ind)<-61.99, ind) = NaN;

            
% save(['lmgvel' adcp '_' num2str(decorr_x/1000) '_' num2str(decorr_y) ...
%     '_' num2str(gsize, '%3.0f') 'km'], 'lmgdata');

save(['lmgvel' adcp '_' num2str(decorr_x/1000) '_' num2str(decorr_y) ], 'lmgdata');
    

fprintf('Total number of transects discarded: %s\n', num2str(r));
fprintf('Number of transects used: %s\n', num2str(length(I) - r));

figure('color','w');
plot(varerr,' o-k', 'linewidth',1.5);
hold on
ylim([0 0.40]);
title(['Decorr. scale(x, z) = ' num2str(decorr_x*1e-3) ' km, ' num2str(decorr_y) ...
    ' m,     Error = ' num2str(error)]);
ylabel('normalized Root Mean Square error');
xlabel(' Transect number');

