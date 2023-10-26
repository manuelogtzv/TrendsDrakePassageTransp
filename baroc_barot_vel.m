% Code for estimating the difference between total transport and baroclinic
% transport referenced at 760 m. Barotropic component is then calculated as
% the depth-averaged residual velocity (total - baroclinic) at each
% distance position and transect, multiplied by water depth and distance.
% 
% Manuel O. Gutierrez Villanueva
% 
% 2021/11/29
%
% 2021/12/14 - Includes nb150 option
% 2023/06/26 - No longer includes the deep estimate, changed the name of
% the output .mat file, saves velocity and transport in a single .mat file.
% 2023/06/27 - Changed the interpolation method. Now uses interp1 instead
% of griddata.

clear all;
close all;
pycmd = '/Users/manuelgutierrez/anaconda3/bin/python';

adcp = 'os38nb';
cont_bar = [-10:2:10]; 
cont_bci = [-4:0.5:4];
cont_tot = [-30:5:30];
cmap_bar = getPyPlot_cMap('seismic', length(cont_bar)-1, [], pycmd);
cmap_bci = getPyPlot_cMap('seismic', length(cont_bci)-1, [], pycmd);
cmap_tot = getPyPlot_cMap('seismic', length(cont_tot)-1, [], pycmd);
datein = datenum(1996, 1, 1);
ptch = [0 0 25 25 0; 1000 -10 -10 1000 1000];
bxlim = [0 25];
map = getPyPlot_cMap('bwr', 2, [], pycmd);



% Loads baroclinic velocities
load(['geostrophic_noaver.mat']);
% geosbaroc.u = geosbaroc.u(:, 1:end, :);
% distubci = 0.5*(geosbaroc.dist(1:end-1, :) + ...
%             geosbaroc.dist(2:end, :));
% distubci = distubci - distubci(1, :);



% Velocities
if strcmp(adcp,'nb150') == 1
    nomz = 210; %maximum depth
    load(['lmgvel' adcp '_50_100.mat']);

    % loads total trasnport data from adcp
    load(['TransportDP_adcp_' adcp 'maxz210.mat']);
    lmgtransp = transpDP; clear transpDP;
else
    nomz = 760;% os38nb max z: 780 or 210
    load(['lmgvel' adcp '_50_300.mat']);
    
    % loads total trasnport data from adcp
    load(['TotalTransport_' adcp 'maxz_' num2str(nomz) '_2.mat']);
    lmgtransp = totaltransp; clear totaltransp;
end

% Keeps nominal depth
lmgdata.u = lmgdata.u(lmgdata.z<=760, :, :);
lmgdata.u_om = lmgdata.u_om(lmgdata.z<=760, :, :);
lmgdata.z = lmgdata.z(lmgdata.z<=760);

% lmgdata.dist = lmgdata.dist(:, :);
% lmgdata.lon = lmgdata.lon(:, 1:length(lmgdata.time));
% lmgdata.lat = lmgdata.lat(:, 1:length(lmgdata.time));
% lmgdata.u = lmgdata.u(1:3:end, :, :);
% lmgdata.z = lmgdata.z(1:3:end);

if strcmp(adcp, 'os38nb');

    % Load velocity offset
    load('trandsmiss_theta.mat')
    
%     % Velocity offset for os38nb
%     lmgdata.u(:, :, lmgdata.northbound==1) = ...
%         lmgdata.u(:, :, lmgdata.northbound==1) - uoffs_north;
%     lmgdata.u(:, :, lmgdata.northbound==0) = ...
%         lmgdata.u(:, :, lmgdata.northbound==0) - uoffs_south;
%     % Finds south and northbound transects
%     indsouth = find(lmgdata.northbound == 1);
%     indnorth = find(lmgdata.northbound == 0);
% 
%     % Adds the uoffs
%     load('trandsmiss_theta.mat');
%     uoffs = mean(abs([uoffs_north, uoffs_south]));
%     lmgdata.u(:, :, indsouth) = lmgdata.u(:, :, indsouth) - uoffs;
%     lmgdata.u(:, :, indnorth) = lmgdata.u(:, :, indnorth) + uoffs;
    % Finds south and northbound transects
    indsouth = find(lmgdata.northbound == 1);
    indnorth = find(lmgdata.northbound == 0);

    % Adds the uoffs
    load('trandsmiss_theta.mat');
    lmgdata.u(:, :, indsouth) = lmgdata.u(:, :, indsouth) - uoffs_south;
    lmgdata.u(:, :, indnorth) = lmgdata.u(:, :, indnorth) - uoffs_north;

    nomz = 760;
    lmgdata.u = lmgdata.u(lmgdata.z<=nomz, :, :);
    lmgdata.u_om = lmgdata.u_om(lmgdata.z<=nomz, :, :);
    lmgdata.z = lmgdata.z(lmgdata.z<=nomz);


    % Makes the transect that start from south to north to now be north to
    % south
    indsouth = find(lmgdata.northbound == 1);
    indnorth = find(lmgdata.northbound == 0);
    for i = 1:length(indsouth);
        [lmgdata.dist(:, indsouth(i)) ii] = ...
            sort(abs(nanmax(lmgdata.dist(:, indsouth(i))) - ...
            lmgdata.dist(:, indsouth(i))) + 12.5);
        lmgdata.u(:, :, indsouth(i)) = lmgdata.u(:, ii, indsouth(i));
        lmgdata.u_om(:, :, indsouth(i)) = lmgdata.u_om(:, ii, indsouth(i));
        lmgdata.lon(:, indsouth(i)) = lmgdata.lon(ii, indsouth(i));
        lmgdata.lat(:, indsouth(i)) = lmgdata.lat(ii, indsouth(i));
    end


end

% % Finds south and northbound transects
%     indsouth = find(lmgdata.northbound == 1);
%     indnorth = find(lmgdata.northbound == 0);
% 
% for i = 1:length(indsouth);
%     [lmgdata.dist(:, indsouth(i)) ii] = ...
%         sort(abs(nanmax(lmgdata.dist(:, indsouth(i))) - ...
%         lmgdata.dist(:, indsouth(i))) + 12.5);
%     lmgdata.u(:, :, indsouth(i)) = lmgdata.u(:, ii, indsouth(i));
%     lmgdata.u_om(:, :, indsouth(i)) = lmgdata.u_om(:, ii, indsouth(i));
%     lmgdata.lon(:, indsouth(i)) = lmgdata.lon(ii, indsouth(i));
%     lmgdata.lat(:, indsouth(i)) = lmgdata.lat(ii, indsouth(i));
% end

% % loads baroclinic transport estimates
% load('TransportDP_baroc_noaverg.mat');

r = 0;

ubc = nan(length(lmgdata.z), ...
    size(geosbaroc.u, 2), size(geosbaroc.u, 3));

% Interpolates baroclinic velocities profiles to adcp depths
for i = 1:size(geosbaroc.u, 2);
    for j = 1:size(geosbaroc.u, 3);
        
        if ~isnan(nanmean(squeeze(geosbaroc.u(:, i, j))));
            ubc(:, i, j) = interp1(geosbaroc.z(:), ...
                squeeze(geosbaroc.u(:, i, j)), lmgdata.z, 'pchip', ...
                'extrap');
        end
    end
end


% Finds coincident xbt and adcp transects
[x y] = meshgrid(geosbaroc.time, lmgdata.time);
[t1 t2] = meshgrid(1:length(geosbaroc.time), 1:length(lmgdata.time));

ind_lmg = t2(abs(x-y)<=5);
ind_geos = t1(abs(x-y)<=5);

% % Only coincident transects
distubci = geosbaroc.dist(:, ind_geos);% Cumulative along-track distance
diffdist = diff(geosbaroc.dist(:, ind_geos));
distubci = geosbaroc.dist_u(:, ind_geos);


clear x y t1 t2

figure('color', 'w');
h3 = tight_subplot(5, 3, [0.01 0.01], [0.10 0.07], [0.15 0.03]);
% ubar1 = nan(length(lmgdata.z), size(lmgdata.lon, 1), ...
%     length(ind_geos));

r = 0;

fig2 = figure;

% Substract baroclinic velocity from total velocity
for i = 1:length(lmgdata.z);
    for j = 1:length(ind_lmg);

        utot = squeeze(lmgdata.u(i, :, ind_lmg(j)))';
        ubci = squeeze(ubc(i, :, ind_geos(j)));

        lontot = lmgdata.lon(:, ind_lmg(j)) - 360;
        lattot = lmgdata.lat(:, ind_lmg(j));
        

        % Applies mask to data
        lonbci = geosbaroc.lon_u(:, ind_geos(j));
        latbci = geosbaroc.lat_u(:, ind_geos(j));
        mask = geosbaroc.lon_u(:, ind_geos(j))./...
            geosbaroc.lon_u(:, ind_geos(j));


        vtot = nan(size(lonbci));

        % Interpolates total velocity to baroclinic velocity grid
%         vtot(~isnan(latbci)) = griddata(lontot(~isnan(lontot)), ...
%             lattot(~isnan(lontot)), utot(~isnan(lontot)), ...
%             lonbci(~isnan(latbci)), latbci(~isnan(latbci)), 'linear');
        vtot(~isnan(latbci)) = interp1(lattot(~isnan(lontot)), utot(~isnan(lontot)), ...
            latbci(~isnan(latbci)), 'linear');
        vtot = vtot;

%         figure(fig2)
%         clf
%         plot(vtot, '-r');
%         hold on
%         plot(ubci, '-b');
%         axis equal; axis([-66   -58   -64   -55])
%         drawnow;
%         pause;


%         vbci = nan(size(lontot));
% 
%         % Interpolates total velocity to baroclinic velocity grid
%         vbci(~isnan(lontot)) = griddata(lonbci(~isnan(lonbci)), ...
%             latbci(~isnan(lonbci)), ubci(~isnan(lonbci)), ...
%             lontot(~isnan(lontot)), lattot(~isnan(lontot)), 'cubic');


        

        % Calculates depth-independent value
%         ubar1(i, :, j) = utot - vbci;
%         ubci1(i, :, j) = vbci;
                ubci1(i, :, j) = ubci;
        utot1(i, :, j) = vtot;
        ubar1(i, :, j) = vtot(:) - ubci(:);

        
%         if sum(~isnan(vtot)) ~= sum(~isnan(ubci))
%             fprintf('\n Depth %3.0f  Time %3.0f\n', num2str(i), num2str(j));
%             clf;
%             plot(vtot, '-r');
%             hold on
%             plot(ubci, '-b');
%             plot(squeeze(ubar1(i, :, j)), '--k');
%             drawnow;
%             pause
% 
%         end

    end

    vbar1 = squeeze(ubar1(i, :, :));
    vbci1 = squeeze(ubci1(i, :, :));
    vtot1 = squeeze(utot1(i, :, :));
    disttot = lmgdata.dist(:, ind_lmg);
    
    
    if mod(i, 7) == 1 & r<5;
        
        r = r + 1;
        
        figure(gcf);
        axes(h3((r*3)-2));
        scatter(reshape(repmat([lmgdata.time(ind_lmg) - ...
            datenum(1996, 1, 1)]/365.25, ...
            size(distubci, 1), 1), size(vbar1(:))), ...
            reshape(distubci, size(vbar1(:))), 8,...
            vbar1(:), 'o', 'filled');
        colormap(h3((r*3)-2), cmap_tot);
        caxis(h3((r*3)-2), [-1.0 1.0]);
        text(h3((r*3)-2), 12, 900, [num2str(lmgdata.z(i)) ' m']);
        box on;
        
        axes(h3((r*3)-1));
        scatter(reshape(repmat([geosbaroc.time(ind_geos) - ...
            datenum(1996, 1, 1)]/365.25, ...
            size(distubci, 1), 1), size(distubci(:))), distubci(:), 8,...
            vbci1(:), 'o', 'filled');
        colormap(h3((r*3)-1), cmap_tot);
        caxis(h3((r*3)-1), [-1.0 1.0]);
        text(h3((r*3)-1), 12, 900, [num2str(lmgdata.z(i)) ' m']);
        box on;

        axes(h3((r*3)));
        scatter(reshape(repmat([lmgdata.time(ind_lmg) - ...
            datenum(1996, 1, 1)]/365.25, ...
            size(distubci, 1), 1), size(distubci(:))), ...
            reshape(distubci, size(vbar1(:))), 8,...
            vtot1(:), 'o', 'filled');
        colormap(h3((r*3)), cmap_tot);
        caxis(h3((r*3)), [-1.0 1.0]);
        text(h3((r*3)), 12, 900, [num2str(lmgdata.z(i)) ' m']);
        box on;
    end
    
    
end

% Colors for different profiles
cmapdate = flipud(getPyPlot_cMap('gnuplot', size(ubar1, 3), [], pycmd));
cmappos = flipud(getPyPlot_cMap('RdYlGn', size(ubar1, 2), [], pycmd));


title(h3(1), 'u_r_e_f [m/s]');
title(h3(2), 'u_g_e_o_s [m/s]');
title(h3(3), 'u_t_o_t [m/s]');
ylabel(h3(10), 'Distance from north [km]');
xlabel([h3(end-2)],...
    ['Years since ' datestr(datenum(1996, 1, 1))]);
xlabel([h3(end-1)],...
    ['Years since ' datestr(datenum(1996, 1, 1))]);
xlabel([h3(end)],...
    ['Years since ' datestr(datenum(1996, 1, 1))]);



xx = squeeze(lmgdata.u(:, :, ind_lmg));
x1 = squeeze(lmgdata.u(:, :, :));

% dist_bar = lmgdata.dist(:, ind_lmg);

% Calculates reference velocity
uref1 = nanmean(ubar1(lmgdata.z>=100 & lmgdata.z<=nomz, :, :), 1);
dprng = lmgdata.z(lmgdata.z>=100 & lmgdata.z<=nomz);
dz = abs(diff(lmgdata.z(2:3)));
dldp = (dprng(end) + dz/2) - (dprng(1) - dz/2);

mask =ubar1(1, :, :)./ubar1(1, :, :);

% Calculates transport due to the barotropic component
transp_bar = squeeze(uref1).*diffdist*dldp*1e3;
cumtransp_bar = cumsum(transp_bar, 1, 'omitnan').*squeeze(mask);
tottransp_bar = nansum(transp_bar, 1);


dzs = lmgdata.z(2)-diff(abs(lmgdata.z(2:3)))/2;
dz = abs(diff(lmgdata.z(2:3)));

% Calculates Ekman transport
uEks = squeeze(ubar1(1, :, :)).*diffdist*1e3*dzs;
uEk = squeeze(nansum(ubar1(lmgdata.z>0 & lmgdata.z<100, :, :), 1)).*diffdist*1e3*dz;
svEk = uEks + uEk;
Ektransp = nansum(svEk, 1);


% transp_bar2 = squeeze(ubar1(end, :, :))...
%     .*diffdist*760*1e3;
% cumtransp_bar2 = cumsum(transp_bar2, 1, 'omitnan').*mask;
% tottransp_bar2 = nansum(transp_bar2, 1);
% tottransp_bar = nanmax(cumtransp_bar);

% Saves data
% geosbar.transp_2 = tottransp_bar2;
% geosbar.transpdist_2 = transp_bar2;
% geosbar.cumtransp_2 = cumtransp_bar2;
% geosref.z = lmgdata.z;
% geosbar.dist = lmgdata.dist;

geosref.adcp = adcp;
geosref.nettransp = tottransp_bar;
geosref.transp = transp_bar;
geosref.cumtransp = cumtransp_bar;
geosref.u = uref1;
geosref.ektransp = svEk;
geosref.eknettransp = Ektransp;
geosref.lon = geosbaroc.lon_u(:, ind_geos);
geosref.lat = geosbaroc.lat_u(:, ind_geos);
geosref.dist = geosbaroc.dist_u(:, ind_geos);
geosref.time = geosbaroc.time(ind_geos);
geosref.doc = {'Reference velocity and transport';...
               'adcp is the adcp used';...
               'transp is the transport [m^3 s^{-1}] per distance bin';...
               'cumtransp is the cumulative transport [m^3 s^{-1}]';...
               'nettransp is the net Drake Passage transport [m^3 s^{-1}]';...
               'ektransp is the Ekman transport [m^3 s^{-1}] per distance bin';...
               'eknettransp is the net Drake Passage Ekman transport [m^3 s^{-1}]';...
               'u is the barotropic velocity component [m s^{-1}]';...
               'dist is the alongtrack distance from the north [km]';...
               'time is the time';...
               'lon and lat are the longitude and latitude'};
           
save(['reference_' adcp '_760.mat'], 'geosref');


figure('color', 'w');
ha1 = tight_subplot(2, 2, [0.10 0.01], [0.13 0.08], [0.15 0.10]);

axes(ha1(1));
p1 = plot(squeeze(nanmean(ubar1, 2)), -lmgdata.z);
hold on;
set(p1, {'Color'}, num2cell(cmapdate, 2));
ha1(1).TickLength = [0.03 0.03];
ylabel('Depth [m]');
xlim([0.01 0.25]);
title('Mean profiles by transect');


axes(ha1(2));
p2 = plot(squeeze(nanmedian(ubar1, 2)), -lmgdata.z);
hold on;
set(p2, {'Color'}, num2cell(cmapdate, 2));
ha1(2).TickLength = [0.03 0.03];
ha1(2).YTickLabel = [];
xlim([0.01 0.25]);
title('Median profiles by transect');


axes(ha1(3));
p3 = plot(squeeze(nanmean(ubar1, 3)), -lmgdata.z);
hold on;
set(p3, {'Color'}, num2cell(cmappos, 2));
ha1(3).TickLength = [0.03 0.03];
xlabel('[m/s]');
ylabel('Depth [m]');
xlim([-0.31 0.63]);
title('Mean profiles by position');

axes(ha1(4));
p4 = plot(squeeze(nanmedian(ubar1, 3)), -lmgdata.z);
hold on;
set(p4, {'Color'}, num2cell(cmappos, 2));
ha1(4).TickLength = [0.03 0.03];
ha1(4).YTickLabel = [];
xlabel('[m/s]');
xlim([-0.31 0.63]);
title('Median profiles by position');


figure('color', 'w');
ha2 = tight_subplot(1, 1, [0.01 0.10], [0.13 0.08], [0.15 0.03]);

axes(ha2(1));
[~, ~, h2] = histf(x1(:), [-1.6:0.1:1.6], ...
    'facecolor', map(2, :));
hold on;
[~, ~, h1] = histf(xx(:), [-0.6:0.05:0.6], ...
    'facecolor', map(1, :));
h1.FaceAlpha = 0.3; h1.EdgeColor = [0.5 0.5 0.5];
h2.FaceAlpha = 0.3; h2.EdgeColor = [0.5 0.5 0.5];
xlabel(['[m/s]']);
ylabel('Frequency');
title(['Cross-transect velocity at ' num2str(lmgdata.z(end))...
    ' m']);
text(-1, 5e4, ['Mean = ' num2str(nanmstd(xx(:)), '%1.2f') ' m/s'], ...
    'color', map(1, :));
text(-1, 4.5e4, ['Median = ' num2str(nanmedian(xx(:)), '%1.2f') ' m/s'], ...
    'color', map(1, :));
text(-1, 2e4, ['Mean = ' num2str(nanmstd(x1(:)), '%1.2f') ' m/s'], ...
    'color', map(2, :));
text(-1, 1.5e4, ['Median = ' num2str(nanmedian(x1(:)), '%1.2f') ' m/s'], ...
    'color', map(2, :));
box on


figure('color','w');
fig1 = gcf;
figure_width = 16.5;
figure_height = 13;
figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
set(fig1,'Visible', figuresVisible)
set(fig1, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
set(fig1, 'PaperPositionMode', 'auto');

h3 = tight_subplot(3, 1, [0.01 0.01], [0.12 0.03], [0.19 0.16]);

axes(h3(1));
p1 = patch(ptch(1, :), ptch(2, :), [0.85 0.85 0.85]);
hold on; 
scatter([reshape(repmat(geosbaroc.time, size(geosbaroc.dist_u, 1), 1), ...
    length(geosbaroc.dist_u(:)), 1)-datein]/365.25, ...
    geosbaroc.dist_u(:), 8, geosbaroc.transp(:)*1e-6, 'o', 'filled');
colormap(h3(1), cmap_bci); caxis(h3(1), [min(cont_bci) max(cont_bci)]); 
col1 = colorbar; 
% dTk = diff(col.Limits)/(2*length(cmap));
% set(col1, 'Ticks', [col.Limits(1)+dTk:2*dTk:col.Limits(2)-dTk], ...
%     'TickLen', 0.02);
ylabel(col1, 'Transport [Sv]');
% ylabel('Distance from North [km]');
axis ij;
axis(h3(1), [bxlim -10 990]); 
set(h3(1), 'YTick', [0:200:1000], 'YTickLabel', ...
    num2str([0:200:1000]'), 'TickLen', [0.01 0.01], 'TickDir', 'out');
text(0, 900, 'a) Geostrophic (760 m)');
% xlabel(['Year since ' datestr(datein)]);
box on;
set(col1, 'TickLen', 0.02, 'Position', [sum(h3(1).Position([1, 3]))+0.1 ...
    h3(1).Position(2) 0.0089 h3(1).Position(4)]);

axes(h3(2));
p2 = patch(ptch(1, :), ptch(2, :), [0.85 0.85 0.85]);
hold on; 
scatter([reshape(repmat(geosbaroc.time(ind_geos), size(distubci, 1), 1), ...
    length(distubci(:)), 1)-datein]/365.25, ...
    distubci(:), 8, transp_bar(:)*1e-6, 'o', 'filled');
colormap(h3(2), cmap_bar); caxis(h3(2), [min(cont_bar) max(cont_bar)]); 
col2 = colorbar; 
% dTk = diff(col.Limits)/(2*length(cmap));
% set(col,'Ticks',[col.Limits(1)+dTk:2*dTk:col.Limits(2)-dTk], ...
%     'TickLen', 0.02);
col2.Ticks = cont_bar(2:end-1)';
ylabel(col2, 'Transport [Sv]');
ylabel('Distance from North [km]');
axis ij;
axis(h3(2), [bxlim -10 990]); 
set(h3(2), 'YTick', [0:200:1000], 'YTickLabel', ...
    num2str([0:200:1000]'), 'TickLen', [0.01 0.01], 'TickDir', 'out');
text(0, 900, 'b) Reference (90-760 m)');
% xlabel(['Year since ' datestr(datein)]);
box on;
set(col2, 'TickLen', 0.02, 'Position', [sum(h3(2).Position([1, 3]))+0.1 ...
    h3(2).Position(2) 0.0089 h3(2).Position(4)]);

axes(h3(3));
p3 = plot([geosbaroc.time-datein]/365.25, geosbaroc.nettransp*1e-6, 'o-k', ...
    'linewidth', 1.5);
hold on;
[hh(1), p4(1)] = addaxis([geosbaroc.time(ind_geos)-datein]/365.25, ...
    geosref.nettransp*1e-6, 'o-r', 'linewidth', 1.5);
[hh(2), p4(2)] = addaxis([geosbaroc.time(ind_geos)-datein]/365.25, ...
    (geosref.nettransp+geosbaroc.nettransp(ind_geos))*1e-6, '^-b', 'linewidth', 1.5);
% [hh(3), p4(3)] = addaxis([lmgtransp.time-datein]/365.25, lmgtransp.transp*1e-6, ...
%     's', 'color', [0.5 0.5 0.5]);
xlabel(['Year since ' datestr(datein)]);
addaxislabel(1, 'Geostrophic transport [SV]');
addaxislabel(2, 'Reference transport [Sv]');
addaxislabel(3, 'Total transport [Sv]');
xlim(h3(3), bxlim);
p4(1).Color(4) = 0.5;
p4(2).Color(4) = 0.5;
p4(3).Color(4) = 0.5;
h3(3).Box = 'on';
h3(3).Position(1) = h3(2).Position(1); 
h3(3).Position(3) = h3(2).Position(3);
hh(1).Position(1) = sum(h3(2).Position([1, 3]))+0.005;
hh(2).Position(1) = h3(2).Position(1) - 0.10;
hh(1).TickLength = [0.03 0.03]; hh(1).XColor = 'r';
hh(2).TickLength = [0.03 0.03]; hh(2).XColor = 'b';

