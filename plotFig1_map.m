% Code for plotting Figure 1 - Area of study, plus mean temperature,
% salinity and velocity fields.
%
% Oct 4/2022
% March 3/2022 - Plots entire streamfunction instead of bathymetry

clear all;
close all;


% Loads palette for bathymetry
load('../Turbulent_mixing/gmt_palette.mat');
adtcolor = map1(1:end-1,:);
pycmd = '/Users/manuelgutierrez/anaconda3/bin/python';
S_unique = [-1.75:0.05:-0.30];


% Contours and colormaps
cont_sal = [33.6:0.05:34.8];
cont_pt = [-2:0.5:8];
cont_dens = [26:0.1:30];
cont_vel = [-50:5:50];
cmap_sal = flipud(getPyPlot_cMap('PiYG', length(cont_sal)-1, [], pycmd));
cmap_pt = flipud(getPyPlot_cMap('RdYlBu', length(cont_pt)-1, [], pycmd));
cmap_vel = getPyPlot_cMap('seismic', length(cont_vel)-1, [], pycmd);
cont_sal_std = [0:0.01:0.11];
cmap_sal_std = viridis(length(cont_sal_std)-1);
cont_pt_std = [0:0.2:1.7];
cmap_pt_std = plasma(length(cont_pt_std)-1);
cmap_ssh = flipud(cmocean('deep', length(S_unique)-1));
%cmap_ssh = (getPyPlot_cMap('ocean', length(S_unique)-1, [], pycmd));%
%cmap_ssh = flipud(viridis(length(S_unique)-1));

% Other variables
datein = datenum(1996, 01, 01, 0, 0 ,0);
freqan = 1/365.25*pi;
freqseman = freqan*2;
ti1 = [datenum(2007, 1 , 1) - datein]/365.25;
ti2 = [datenum(2015, 1, 1) - datein]/365.25; 
alpha = 0.05;
coldis = distinguishable_colors(5);
gamma = 1;
season_opt = 1;
depths_psi = [100 200];% Averages mean geostrophic velocities and psi over these two depths


% Loads data
load('lmgvelos38nb_50_300.mat');

% Makes the transect that start from south to north to now be north to
% south
indsouth = find(lmgdata.northbound == 1);
indnorth = find(lmgdata.northbound == 0);

% Adds the uoffs
load('trandsmiss_theta.mat');
uoffs = mean(abs([uoffs_north, uoffs_south]));
lmgdata.u(:, :, indsouth) = lmgdata.u(:, :, indsouth) - uoffs;
lmgdata.u(:, :, indnorth) = lmgdata.u(:, :, indnorth) + uoffs;
nomz = 970;
lmgdata.u = lmgdata.u(lmgdata.z<=nomz, :, :)*100;
lmgdata.u_om = lmgdata.u_om(lmgdata.z<=nomz, :, :);
lmgdata.z = lmgdata.z(lmgdata.z<=nomz);

% Makes all transect to start from north
for i = 1:length(indsouth);
   [lmgdata.dist(:, indsouth(i)) ii] = ...
      sort(abs(nanmax(lmgdata.dist(:, indsouth(i))) - ...
      lmgdata.dist(:, indsouth(i))) + 12.5);
   lmgdata.u(:, :, indsouth(i)) = lmgdata.u(:, ii, indsouth(i));
   lmgdata.u_om(:, :, indsouth(i)) = lmgdata.u_om(:, ii, indsouth(i));
   lmgdata.lon(:, indsouth(i)) = lmgdata.lon(ii, indsouth(i));
   lmgdata.lat(:, indsouth(i)) = lmgdata.lat(ii, indsouth(i));
end
os38 = lmgdata;

load('geosbarocDP_noaver.mat'); %xbt lines
load('omgeosvel_nb150.mat'); %mean geostrophic vel

% Renames T, S, gamma and other variables
sal = geosbaroc.sal; 
dist = geosbaroc.dist(1:end-1, :);
gamma_surfaces = [26:0.05:29];
gamma_surfaces = gamma_surfaces(1:28);
gpot = geosbaroc.gpot;
ubci = geosbaroc.u;

% Pre-allocates memory
pt0 = nan(size(sal));
dens0 = pt0;
gamma_n = pt0;
sal_gamma = nan(length(gamma_surfaces), size(sal, 2), size(sal, 3));
temp_gamma = sal_gamma;
p_gamma = sal_gamma;

% Calculates potential temperature, density and neutral density
for i = 1:length(geosbaroc.time)
    pt0(:, :, i) = gsw_pt0_from_t(squeeze(sal(:, :, i)), ...
        squeeze(geosbaroc.temp(:, :, i)), geosbaroc.z(:));

    if gamma ~= 1; %calculates potential density
        dens0(:, :, i) = sw_pden(squeeze(sal(:, :, i)), ...
            squeeze(geosbaroc.temp(:, :, i)), geosbaroc.z(:), 0) - 1000;
    else % calculates gamma
        gamma_n(:, :, i) = eos80_legacy_gamma_n(squeeze(sal(:, :, i)), ...
            squeeze(geosbaroc.temp(:, :, i)), repmat(geosbaroc.z, 1, size(sal, 2)),...
            repmat(geosbaroc.lon(:, i)', size(geosbaroc.z, 1), 1),...
            repmat(geosbaroc.lat(:, i)', size(geosbaroc.z, 1), 1));
    end

end

% Calculate gamma (neutral density) or potential density
if gamma == 1;
    dens = gamma_n;
    coltext = '\gamma_n [kg m^-^3]';
else
    dens = dens0;
    coltext = '\sigma_{\theta} - 1000 [kg m^-^3]';
end

pt0(:, 77:end, :) = NaN;
sal(:, 77:end, :) = NaN;
dens(:, 77:end, :) = NaN;

% Mean streamwise position of the ACC fronts
saf_psi = [-0.40 -0.60];%Subantarctic Front
pf_psi = [-1.00  -1.30];%Polar Front
saccf_psi = [-1.60 -1.65];%Southern ACC Front

mrl_coor = [-64.9422  -55.0384; -62.4422  -62.3884];
dir_folder = './aviso/';

% Loads MDT
% Averages time-mean streamfunction, velocity over specific depth ranges
% Extracts Mean Dynamic Topography (MDT) and reduces to DP area
ssh_mdt = ncread([dir_folder 'niiler_maximenko_1992_2012_mdt.nc'],...
    'MDOT')*1e-2 - 0.6;
lon_mdt = ncread([dir_folder 'niiler_maximenko_1992_2012_mdt.nc'],...
    'LONN319_N207');
lat_mdt = ncread([dir_folder 'niiler_maximenko_1992_2012_mdt.nc'],...
    'LAT37_109');

% Drake Passage area
wlon = -68; elon = -56; 
ind_lon = find(lon_mdt > wlon & lon_mdt < elon);

slat = -65; nlat = -54.25; 
ind_lat = find(lat_mdt > slat & lat_mdt < nlat);

lon_mdt = lon_mdt(ind_lon);
lat_mdt = lat_mdt(ind_lat);
ssh_mdt = ssh_mdt(ind_lon,ind_lat);
ssh_mdt = permute(ssh_mdt,[2,1]);

% Plot levels

% figure
fig1 = figure('color', 'w');
figure_width = 18;
figure_height = 12;
figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
set(fig1,'Visible', figuresVisible)
set(fig1, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
set(fig1, 'PaperPositionMode', 'auto');

ax1 = subplot(3, 4, [1 2 5 6 9 10], 'align');
ax1.Position(1) = 0.06; ax1.FontName = 'SansSerif';
m_proj('lambert', 'lon', [-67.90 -56.90], 'lat', [-65 -54.4]); hold on;
colormap(ax1, cmap_ssh);
caxis(ax1, [min(S_unique) max(S_unique)]);            
%m_etopo2('shadedrelief', 'gradient', 10);
m_gshhs_i('patch', [.4 .4 .4], 'edgecolor', 'none');
m_grid('linestyle', 'none', 'tickdir', 'out', 'linewidth', 1.5,...
    'xtick', [-65:5:-55], 'fontsize', 7);
ax = m_contfbar(.87, [.35 .75], [S_unique([1, end])], S_unique, ...
    'edgecolor', 'none', 'endpiece', 'none');
ax.YTick = [S_unique(2:4:end)];
ylabel(ax, 'SSH [m]', 'fontsize', 7);
ax.FontSize = 7;
ax.FontName = 'SansSerif';
ax.TickLength = [0.03 0.03];
ax.TickDir = 'out';
ax.Position = [0.3900    0.1200    0.0109    0.2500];
[om h1] = m_contourf(lon_mdt, lat_mdt, ...
    ssh_mdt, S_unique, 'EdgeColor', 'none');
[CS,CH] = m_etopo2('contour',[-1000 -200], 'color','black', ...
    'linewidth', 0.5);
clabel(CS, CH, 'fontsize', 5, 'fontname', 'SansSerif', 'labelspacing', 500);
[~, ~] = m_contour(lon_mdt, lat_mdt, ...
    ssh_mdt, [saf_psi(2):0.1:saf_psi(1)], 'Color', rgb('Goldenrod'), ...
    'linewidth', 1);
[~, ~] = m_contour(lon_mdt, lat_mdt, ...
    ssh_mdt, [pf_psi(2):0.1:pf_psi(1)], 'Color', rgb('DarkGray'), ...
    'linewidth', 1);
[~, ~] = m_contour(lon_mdt, lat_mdt, ...
    ssh_mdt, [saccf_psi(2):0.05:saccf_psi(1)], 'Color', 'm', ...
    'linewidth', 1);
c1 = m_plot(-180, -180, '-', 'Color', rgb('DimGray'), 'linewidth', 2);
lmg = m_plot(os38.lon(:)-360, os38.lat(:), '-', 'linewidth', 1.5, ...
    'Color', rgb('Black'));
xbt = m_plot(geosbaroc.lon(:), geosbaroc.lat(:), '-', 'linewidth', 0.5, ...
    'Color', rgb('Red'));
leg = legend([lmg; xbt], 'sADCP transects', 'XBT/XCTD transects');
set(leg, 'location', 'NorthEast', 'fontsize', 7, 'EdgeColor', 'none', ...
    'Position', [0.1512    0.1272    0.1765    0.0647]);
leg.ItemTokenSize = [10; 10];
m_text(-64, -55.5, 'SAF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'color', 'k');
m_text(-61.1, -58.2, 'PF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'color', 'k');
m_text(-58.9, -60.5, 'SACCF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'color', [0.6 0.6 0.6]);

ax2 = axes;
ax2.Position = [0.22 0.73 0.18 0.18];
ax2.FontName = 'SansSerif';
m_proj('azimuthal equal-area', 'radius', 55, 'lat', -70, ...
    'long', -65, 'rot', 0);
colormap(ax2, adtcolor);
caxis(ax2, [-8000 000]);       
m_etopo2('shadedrelief', 'gradient', 10);
m_coast('patch', [.4 .4 .4], 'edgecolor', 'none');
m_grid('xticklabel', [], 'yticklabel', [], 'linestyle', 'none',...
    'ytick', [-70:30:90], 'linewidht', 0.5, 'xtick', [-160:60:180]);
hold on;
p1 = m_plot([-68 -68 -56 -56 -68], [-68 -54 -54 -68 -68], '-r', ...
    'linewidth', 1);

ax3 = subplot(3, 4, [3, 4]);
ax3Pos = ax3.Position;
ax3.FontName = 'SansSerif';
contourf(nanmean(os38.dist(1:39, :)'), -os38.z, ...
    squeeze(nanmean(os38.u(:, 1:39, :), 3)), cont_vel, 'edgecolor', 'none');
hold on;
axis(ax3, [0 840 -760 0]);
p1 = plot([820 745], [0 0], '-', 'linewidth', 4, 'color', 'm');
p2 = plot([300 500], [0 0], '-', 'linewidth', 4, 'color', rgb('DarkGray'));
p3 = plot([35 135], [0 0], '-', 'linewidth', 4, 'color', rgb('Goldenrod'));
[cc hh] = contour(nanmean(dist'), -geosbaroc.z, ...
    squeeze(nanmean(dens(:, 1:end-1, :), 3)), cont_dens, ...
    'color', [0.5 0.5 0.5]);
clabel(cc, hh, 'labelspacing', 140, 'fontsize', 7, 'fontname', 'SansSerif');
caxis(ax3, [cont_vel(1)  cont_vel(end)]);
colormap(ax3, cmap_vel);
set(ax3, 'TickLength', [0.01 0.01], 'TickDir', 'out', 'XTickLabel', [], ...
    'fontsize', 7);
ylabel('Depth [m]', 'fontsize', 7, 'fontname', 'SansSerif');
col3 = colorbar(ax3);
set(col3, 'fontsize', 7, 'TickLength', 0.03, 'Ticks', cont_vel(1:2:end)',...
    'TickLabels', num2str(cont_vel(1:2:end)'), 'TickDir', 'out');
ylabel(col3, 'u_{tot} [cm s^{-1}]', 'fontsize', 7, 'fontname', 'SansSerif');
text(300, -700, ['February 2005 - December 2019'], 'fontsize', 7, ...
    'fontname', 'SansSerif');


ax4 = subplot(3, 4, [7, 8]);
ax4Pos = ax4.Position;
ax4.FontName = 'SansSerif';
contourf(nanmax(dist'), -geosbaroc.z, ...
    squeeze(nanmean(pt0(:, 1:end-1, :), 3)), cont_pt, 'edgecolor', 'none');
hold on;
axis(ax4, [0 840 -760 0]);
p1 = plot([820 745], [0 0], '-', 'linewidth', 4, 'color', 'm');
p2 = plot([300 500], [0 0], '-', 'linewidth', 4, 'color', rgb('DarkGray'));
p3 = plot([35 135], [0 0], '-', 'linewidth', 4, 'color', rgb('Goldenrod'));
[cc hh] = contour(nanmean(dist'), -geosbaroc.z, ...
    squeeze(nanmean(dens(:, 1:end-1, :), 3)), cont_dens, ...
    'color', [0.5 0.5 0.5]);
clabel(cc, hh, 'labelspacing', 140, 'fontsize', 7, 'fontname', 'SansSerif');
caxis(ax4, [cont_pt(1)  cont_pt(end)]);
colormap(ax4, cmap_pt);
set(ax4, 'TickLength', [0.01 0.01], 'TickDir', 'out', 'XTickLabel', [], ...
    'fontsize', 7);
ylabel('Depth [m]', 'fontsize', 7, 'fontname', 'SansSerif');
col4 = colorbar(ax4);
set(col4, 'fontsize', 7, 'TickLength', 0.03, 'Ticks', cont_pt(1:2:end)',...
    'TickLabels', num2str(cont_pt(1:2:end)'), 'TickDir', 'out');
ylabel(col4, '\theta [^oC]', ...
    'fontsize', 7, 'fontname', 'SansSerif');
text(300, -700, ['February 1997 - October 2019'], 'fontsize', 7, ...
    'fontname', 'SansSerif');

ax5 = subplot(3, 4, [11, 12]);
ax5Pos = ax5.Position;
ax5.FontName = 'SansSerif';
contourf(nanmean(dist'), -geosbaroc.z, ...
    squeeze(nanmean(sal(:, 1:end-1, :), 3)), cont_sal, 'edgecolor', 'none');
hold on;
axis(ax5, [0 840 -760 0]);
p1 = plot([820 745], [0 0], '-', 'linewidth', 4, 'color', 'm');
p2 = plot([300 500], [0 0], '-', 'linewidth', 4, 'color', rgb('DarkGray'));
p3 = plot([35 135], [0 0], '-', 'linewidth', 4, 'color', rgb('Goldenrod'));
[cc hh] = contour(nanmean(dist'), -geosbaroc.z, ...
    squeeze(nanmean(dens(:, 1:end-1, :), 3)), cont_dens, ...
    'color', [0.5 0.5 0.5]);
clabel(cc, hh, 'labelspacing', 140, 'fontsize', 7, 'fontname', 'SansSerif');
caxis(ax5, [cont_sal(1)  cont_sal(end)]);
colormap(ax5, cmap_sal);
set(ax5, 'TickLength', [0.01 0.01], 'TickDir', 'out', ...
    'fontsize', 7);
ylabel('Depth [m]', 'fontsize', 7, 'fontname', 'SansSerif');
col5 = colorbar(ax5);
set(col5, 'fontsize', 7, 'TickLength', 0.03, 'Ticks', cont_sal(1:4:end)',...
    'TickLabels', num2str(cont_sal(1:4:end)'), 'fontname', 'SansSerif', ...
    'TickDir', 'out');
ylabel(col5, '{\itS}', 'fontsize', 7, 'fontname', 'SansSerif');
xlabel('Distance from North [km]', 'fontsize', 7, 'fontname', 'SansSerif');
text(300, -700, ['February 1997 - October 2019'], 'fontsize', 7, ...
    'fontname', 'SansSerif');

ax3.Position = ax3Pos;
col3.Position = [ax3.Position(1)+ax3.Position(3)+0.007 ...
    ax3.Position(2) 0.015 ax3.Position(4)];
cbarrow(col3, cmap_vel);
ax4.Position = ax4Pos;
col4.Position = [ax4.Position(1)+ax4.Position(3)+0.007 ...
    ax4.Position(2) 0.015 ax4.Position(4)];
cbarrow(col4, cmap_pt);
ax5.Position = ax5Pos;
col5.Position = [ax5.Position(1)+ax5.Position(3)+0.007 ...
    ax5.Position(2) 0.015 ax5.Position(4)];
cbarrow(col5, cmap_sal);
annotation(fig1, 'textbox', [0.5503    0.9367    0.0480    0.0353],...
    'String', 'SAF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', [0.6894    0.9367    0.0480    0.0353],...
    'String', 'PF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', [0.8473    0.9367    0.0480    0.0353],...
    'String', 'SACCF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', [0.5503    0.6359    0.0480    0.0353],...
    'String', 'SAF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', [0.6894    0.6359    0.0480    0.0353],...
    'String', 'PF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', [0.8473    0.6359    0.0480    0.0353],...
    'String', 'SACCF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', [0.5503    0.3350    0.0480    0.0353],...
    'String', 'SAF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', [0.6894    0.3350    0.0480    0.0353],...
    'String', 'PF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', [0.8473    0.3350    0.0480    0.0353],...
    'String', 'SACCF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', [ax1.Position(1)-0.05 sum(ax1.Position([2,4])) ...
    0.0480 0.0353], 'String', 'a', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig1, 'textbox', [ax3.Position(1)-0.06 sum(ax3.Position([2,4])) ...
    0.0480 0.0353], 'String', 'b', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig1, 'textbox', [ax4.Position(1)-0.06 sum(ax4.Position([2,4])) ...
    0.0480 0.0353], 'String', 'c', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig1, 'textbox', [ax5.Position(1)-0.06 sum(ax5.Position([2,4])) ...
    0.0480 0.0353], 'String', 'd', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
print('-dpdf', './FiguresPaper/fig1_map_mean.pdf', '-r500');

