% Program to load different XBT/XCTD transects from the LMG starting from 
% 1996 till 2019 (XBT) or 2016 (XCTD). Program seperates data into 
% transects and checks for the conditions below:
%
% a) Separates into southbound and northbound transects
% b) Checks that each transects is completed in <4 days
% c) Checks that transects are reach the 500 m isobath at both ends (south
% america and antarctica)
% d) Removes transects with gaps larger than 150 km
% 
% Unlike the ADCP data, the XBT/XCTD data does not require objective
% mapping to fill the gaps in. Transects have been interpolated and filled
% already.
%
% Some transects have been manually discarded since large gaps were found
% or the transects orientation changed several times. 
%
% Manuel O. Gutierrez V.
% 2021/11/22
%
% 2021/11/29 - Allows for calculating geostrophic quantities with raw data
% 2022/02/25 - Applies mask to longitude and latitude
% 2022/05/19 - Saves temperature and salinity data, and velocity and
% transport in the same .mat file.



clear all; 
close all;

gsize = 25; %km
cont_sal = [33.0:0.1:35.0];%salinity contours
cont_temp = [-3:0.5:8];%temperature contours
cont_vel = [-0.3:0.05:0.3];%velocity contours
cont_ster = [0 0.025 0.050 0.075 0.1:0.1:0.8];%steric height
boxlim = [0    1.2   -0.8   -0.0460]*1e3;
binaver = 0;

% loads colormaps
load('./geos_cmap.mat')


if binaver == 0
    arch = 'noaver';
else
    arch = 'aver';
end
    

% loads transects
load('../Datasets/ax22om_1996_2019.mat');

% Flips data so transects start in the north
ax22om.lat = flip(ax22om.lat);
ax22om.lon = flip(ax22om.lon, 1);
ax22om.temp = flip(ax22om.temp, 1);
ax22om.sal = flip(ax22om.sal, 1);


fprintf('\n\nNumber of transects available: %3.0f \n\n', length(ax22om.time));
fprintf('Starting date: %s \n', datestr(ax22om.time(1)));
fprintf('Ending date: %s \n\n', datestr(ax22om.time(end)));


% Removes data away of specific boundaries
indgd = ax22om.lat <= -54.9 & ax22om.lat >= -64.0;
ax22om.lat = ax22om.lat(indgd);
ax22om.lon = ax22om.lon(indgd, :);
ax22om.temp = ax22om.temp(indgd, :, :);
ax22om.sal = ax22om.sal(indgd, :, :);

% Loads bathymetry to determine which cruises are good for transport
b = load('../Datasets/drakembf_avgd'); 
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
figure_width = 36;
figure_height = 12;
figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
set(fig1,'Visible', figuresVisible);
set(fig1, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height]);
set(fig1, 'PaperPositionMode', 'auto');
ha = tight_subplot(2, 3, [0.08 0.08],[0.20 0.05], [0.20 0.05]);


% Preallocates data
distgeo = nan(length(ax22om.lat), size(ax22om.sal, 3));
latgeo = distgeo;
longeo = distgeo;
ugeo = nan(length(ax22om.depth_sal), size(distgeo, 1)-1, size(ax22om.sal, 3));
geoanom = nan(length(ax22om.depth_sal), size(distgeo, 1), size(ax22om.sal, 3));
sals = nan(length(ax22om.depth_sal), size(distgeo, 1), size(ax22om.sal, 3));
temps = nan(length(ax22om.depth_sal), size(distgeo, 1), size(ax22om.sal, 3));


for i = 1:size(ax22om.sal, 3);
    
    clear data
    
    % Data for each transect
    data.sal = squeeze(ax22om.sal(:, :, i));
    data.temp = squeeze(ax22om.temp(:, 1:length(ax22om.depth_sal), i));
    data.time = ax22om.time(i);
    data.lon = ax22om.lon(:, i);
    data.lat = ax22om.lat(:);
    data.bi = interp2(b.lon, b.lat, b.bathyf, data.lon, data.lat);
    data.z = ax22om.depth_sal;
    
      
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
                & (data.lon+360 > 293.5);
            data.lat = data.lat(dix);
            data.lon = data.lon(dix);
            data.sal = data.sal(dix, :);
            data.temp = data.temp(dix, :);
            data.z = ax22om.depth_sal;
            data.bi = data.bi(dix);
                                 
            
            % Calculates again alongtrack distance
            [data.dist data.ang] = sw_dist(data.lat, data.lon, 'km');
            data.ang = mod(data.ang, 360);
            data.ang = [data.ang(1); data.ang]';
          
            
            % Checks if it transect is south or northbound, and corrects
            % for the reversing orientation
            if (data.lat(1) > -60);% & (sum(diff(data.lat)>0)>0); %southbound
                if sum(diff(data.lat)>0)>0;
                    data.dist(diff(data.lat) > 0) = ...
                        data.dist(diff(data.lat)>0)*-1;
                end
%                 
            elseif (data.lat(1) < -60);% & (sum(diff(data.lat)<0)>0); %northbound
                if sum(diff(data.lat)<0)>0;
                    data.dist(diff(data.lat) < 0) = ...
                        data.dist(diff(data.lat)<0)*-1;
                end
            end
            
            
            % Cumulative distance
            data.dist = [0; cumsum(data.dist)];
            
            if binaver == 1;
                % Bin averages every certain distance
                [bins, lataver] = line_bin(data.dist, ...
                    repmat(data.lat', size(data.sal, 1), 1), gsize);
                [bins, lonaver] = line_bin(data.dist, ...
                    repmat(data.lon', size(data.sal, 1), 1), gsize);
                [bins, sal] = line_bin(data.dist, data.sal', gsize);
                [bins, temp] = line_bin(data.dist, data.temp', gsize);
                [bins, biaver] = line_bin(data.dist, ...
                    repmat(data.bi', size(data.sal, 1), 1), gsize);
%                 [bins, distaver] = line_bin(data.dist, data.diust, gsize);

                lataver = lataver(1,:);
                lonaver = lonaver(1,:);
                bins = bins';
                
                
            else
                
                bins = data.dist;
                lataver = data.lat';
                lonaver = data.lon';
                sal = data.sal';
                temp = data.temp';
                biaver = data.bi;
            end
                

                
            
            % Geopotential anomaly m2 s-2
            gpanaver = sw_gpan(sal, temp, ...
                repmat(data.z, 1, size(sal, 2)));
            gpanaver = gpanaver(end, :) - gpanaver; %Referenced to 760
            
            f = sw_f(0.5*(lataver(1:end-1) + lataver(2:end)));
            ugeos = gpanaver(:, 1:end-1) - gpanaver(:, 2:end);
            ugeos = -ugeos./repmat(f, size(ugeos, 1), 1);
            ugeos = ugeos./diff(bins*1e3)'; %m/s
            g = sw_g(repmat(lataver, size(data.z, 1), 1), ...
                repmat(data.z, 1, size(lataver, 2)));
            
                       
             % New grid (depth x dist)
            [xx1, yy1] = meshgrid(bins(1):gsize:bins(end), data.z);

            indgd = ~isnan(sal.*temp);
            
            
            % Plotting
            
            if mod(s, 5) == 1;
                figure(fig1);
                axes(ha(1));
                if s == 1;
                    contour(b.lon, b.lat, b.bathyf, [0 0], ...
                        'linewidth', 1.5, 'color', 'k');
                    hold on;
                    contour(b.lon, b.lat, b.bathyf, -[500 1000], ...
                        'linewidth', 1, 'color', [0.5 0.5 0.5]);
                    axis(ha(1), [-69 -55 -67 -54]); axis equal;
                    ha(1).TickLength = [0.02 0.02];
                end
                p1 = plot(data.lon, data.lat, '*r');
                p2 = plot(lonaver, lataver, '.k');
                title(datestr(data.time));
            
            
                axes(ha(2));
                cla;
                [~, c1] = contourf(bins, -data.z, sal, cont_sal, ...
                    'edgecolor', 'none');
                hold on
                bat = plot(bins, biaver, 'k');
                caxis(ha(2), [cont_sal(1) cont_sal(end)]);
                colormap(ha(2), cmap_sal);
                caxis(ha(2), [cont_sal(1) cont_sal(end)]);
                colormap(ha(2), cmap_sal);
                if s == 1;
                    title(ha(2), 'Salinity ');
                    ylabel(ha(2), 'Depth [m]');
%                     xlabel(ha(2), 'Alongtrack distance [km]');
                    ha(2).XTickLabel = [];
                    axis(ha(2), boxlim);
                    colsal = colorbar('Position', [0.64 0.71 0.0040 0.09]);
                    ylabel(colsal, '');
                    ha(2).TickLength = [0.02 0.02];
                end
            
                axes(ha(3));
                cla;
                [~, c2] = contourf(bins, -data.z, temp, cont_temp, ...
                    'edgecolor', 'none');
                hold on
                [h1, c6] = contour(bins, -data.z, temp, [2 2], ...
                    'edgecolor', 'k');
                clabel(h1, c6, 'fontsize', 10);
                bat = plot(bins, biaver, 'k');
                caxis(ha(3), [cont_temp(1) cont_temp(end)]);
                colormap(ha(3), cmap_temp);
                caxis(ha(3), [cont_temp(1) cont_temp(end)]);
                colormap(ha(3), cmap_temp);
                if s == 1;
                    title(ha(3), 'Temperature ');
%                     ylabel(ha(3), 'Depth [m]');
                    xlabel(ha(3), 'Alongtrack distance [km]');
                    axis(ha(3), boxlim);
                    coltmp = colorbar('Position', ...
                        [0.9183 0.8243 0.0040 0.1100]);
                    ylabel(coltmp, '^oC');
                    ha(3).TickLength = [0.02 0.02];
                end
                        
                axes(ha(4));
                cla;
                [~, c3] = contourf(bins, -data.z, gpanaver./g, 10,...
                    'edgecolor', 'none');
                hold on
                bat = plot(bins, biaver, 'k');
                caxis(ha(4), [cont_ster(1) cont_ster(end)]);
                colormap(ha(4), cmap_sterh);
                caxis(ha(4), [cont_ster(1) cont_ster(end)]);
                colormap(ha(4), cmap_sterh);
                if s == 1;
                    title(ha(4), 'Steric Height ');
                    ylabel(ha(4), 'Depth [m]');
                    xlabel(ha(4), 'Alongtrack distance [km]');
                    axis(ha(4), boxlim);
                    colster = colorbar('Position', ...
                        [0.3621 0.3629 0.0049 0.1037]);
                    ylabel(colster, '[m]');
                    ha(4).TickLength = [0.02 0.02];
                end
            
                axes(ha(5));
                cla;
                [~, c4] = contourf(0.5*(bins(1:end-1)+bins(2:end)),...
                    -data.z, ugeos, cont_vel, 'edgecolor', 'none');
                hold on
                [~, c5] = contour(0.5*(bins(1:end-1)+bins(2:end)),...
                    -data.z, ugeos, [0 0], 'color', 'k', 'linewidth', 1.5);
                bat = plot(bins, biaver, 'k');            
                caxis(ha(5), [cont_vel(1) cont_vel(end)]);
                colormap(ha(5), cmap_vel);
                caxis(ha(5), [cont_vel(1) cont_vel(end)]);
                colormap(ha(5), cmap_vel);
                if s == 1;
                    title(ha(5),...
                        ['u_g_e_o  referenced at ' num2str(data.z(end))...
                        ' m']);
                    xlabel(ha(5), 'Alongtrack distance [km]');
%                     ylabel(ha(5), 'Depth [m]');
                    axis(ha(5), boxlim);
                    colvel = colorbar('Position', ...
                        [0.6227 0.3876 0.0040 0.0775]);
                    ylabel(colvel, '[m/s]');
                    ha(5).TickLength = [0.02 0.02];
                end
            
                drawnow
                delete(ha(6));
                pause(0.1)
%                 print('-depsc', ['./Figures/geosbaroc_' ...
%                     datestr(floor(nanmean(data.time)), 'yyyymmdd')...
%                     '_' arch]);
                delete([p1; p2; c1; c2; c3; c4; c5; c6])
            end
            
            
            % Saves data
            temps(:, 1:length(bins), s) = temp;
            sals(:, 1:length(bins), s) = sal;
            ugeo(:, 1:length(bins)-1, s) = ugeos;
            geoanom(:, 1:length(bins), s) = gpanaver;
            distgeo(1:length(bins), s) = bins;
            longeo(1:length(bins), s) = lonaver;
            latgeo(1:length(bins), s) = lataver;
            geosbaroc.time(1, s) = nanmean(data.time);
                       
        else
            fprintf('Transect discarded: %s \n', ...
                datestr(nanmean(data.time)));
            r = r + 1;
        end
        
    else
        fprintf('Transect discarded: %s \n', datestr(nanmean(data.time)));
        r = r + 1;
    end
end

% Applies mask to lon and lat
mask = squeeze(sum(geoanom, 1));
mask = mask./mask;
longeo = longeo.*mask;
latgeo = latgeo.*mask;
masku = ugeo(:, :, 1:length(geosbaroc.time))./ugeo(:, :, 1:length(geosbaroc.time));
masku = squeeze(masku(1, :, :));

% Remove anomalous high westward velocities
% ugeo(:, 63, 106) = 0.5*(ugeo(:, 62, 106) + ugeo(:, 64, 106));
% ugeo(data.z>694, 1, :) = NaN;

% Calculates transport due to the geostrophic component
diffdist = diff(distgeo(:, 1:length(geosbaroc.time)));
transp = squeeze(trapz(data.z, ugeo(:, :, 1:length(geosbaroc.time)), 1))...
    .*diffdist*1e3;
cumtransp = cumsum(transp, 1, 'omitnan').*masku;
tottransp = nansum(transp, 1);

geosbaroc.temp = temps(:, :, 1:length(geosbaroc.time));
geosbaroc.sal = sals(:, :, 1:length(geosbaroc.time));
geosbaroc.u = ugeo(:, :, 1:length(geosbaroc.time));
geosbaroc.gpot = geoanom(:, :, 1:length(geosbaroc.time));
geosbaroc.dist = distgeo(:, 1:length(geosbaroc.time));
geosbaroc.dist_u = 0.5*(distgeo(1:end-1, 1:length(geosbaroc.time)) + ...
    distgeo(2:end, 1:length(geosbaroc.time)));
geosbaroc.lon = longeo(:, 1:length(geosbaroc.time));
geosbaroc.lon_u = 0.5*(geosbaroc.lon(1:end-1, :) + geosbaroc.lon(2:end, :));
geosbaroc.lat = latgeo(:, 1:length(geosbaroc.time));
geosbaroc.lat_u = 0.5*(geosbaroc.lat(1:end-1, :) + geosbaroc.lat(2:end, :));
geosbaroc.transp = transp;
geosbaroc.cumtransp = cumtransp;
geosbaroc.nettransp = tottransp;
geosbaroc.z = data.z;
geosbaroc.binaver = binaver;
geosbaroc.doc = {'temp is in-situ temperature [oC]';
                 'sal is salinity [ups]';
                 'u is across-transect velocity [m s^{-1}] (positive to east)';
                 'gpot is geopotential anomaly [m^2 s^{-2}] referenced to 760 m';
                 'dist is along track distance [km] for temp, sal and gpot';
                 'lon and lat are coordinates averaged for each 25 km bin';
                 'z is depth [m]'
                 'binaver = 1 indicates that alongtrack 25-km averages';
                 '_u is longitude, latitude and distance for velocity';
                 'transp is the transport [m^3 s^{-1}] per distance box';
                 'cumtransp is the cumulative transport [m^3 s^{-1}]';
                 'nettransp is the net Drake Passage transport [m^3 s^{-1}]'};
                
%saves data
% save(['geostrophic_' arch '.mat'], 'geosbaroc');
    

fprintf('Total number of transects discarded: %s\n', num2str(r));
fprintf('Number of transects used: %s\n', num2str(s));


