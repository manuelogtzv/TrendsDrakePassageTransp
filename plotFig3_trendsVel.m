% Code for plotting Figure 3 - Trends in the cross-transect velocity.
%
% Oct 10/2022 Manuel O. Gutierrez-Villanueva
%
% 2023/06/27 - Rename some variables for consistency.
% 2023/06/28 - Uses the most up-to-date files. Reference velocity was
% corrected for interpolation issues with griddata.
% 2023/07/15 - Removed trends in reference

clear all;
close all;

% pycmd = '/Users/manuelgutierrez/anaconda3/bin/python';

% Define variables
datein = datenum(1996, 01, 01, 0, 0 ,0);
alpha = 0.05;
gsize = 25;

cont = [-1.8:0.2:1.8];%vel contours

% Loads colormaps
load ./colormaps/Fig3_cmap.mat

% cmap = getPyPlot_cMap('seismic', length(cont)-1, [], pycmd);

din1 = datenum(1999, 11, 1);
din2 = datenum(2005, 10, 1);
dfin1 = datenum(2019, 4, 30);

season_opt = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loads os38 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['./Datasets/lmgvelos38nb_50_300.mat']);
% Finds south and northbound transects
indsouth = find(lmgdata.northbound == 1);
indnorth = find(lmgdata.northbound == 0);

% Adds the uoffs
load('./Datasets/trandsmiss_theta.mat');
lmgdata.u(:, :, indsouth) = lmgdata.u(:, :, indsouth) - uoffs_south;
lmgdata.u(:, :, indnorth) = lmgdata.u(:, :, indnorth) - uoffs_north;

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
    
% dist = [lmgdata.dist(1):gsize:size(lmgdata.dist,1)*gsize ];
mask = lmgdata.dist./lmgdata.dist;

total_760 = lmgdata;
% total_760.dist = dist; 
clear lmgdata uoffs* indn* inds* i ii theta* dist;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% GEOSTROPHIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['./Datasets/geostrophic_noaver.mat']);

geosbci_760 = geosbaroc;
clear geosbaroc 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% REFERENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['./Datasets/reference_os38nb_760.mat' ]);
geosref_760 = geosref; clear geosref;
clear geosbaroc distref;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Removes spurious data due to bathymetry mask
total_760.u(total_760.z>694, 1, :) = NaN;
geosbci_760.u(geosbci_760.z>694, 1, :) = NaN;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Trends velocity  %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Renames variables for simplicity
% total
utotal_760 = total_760.u(:, :, total_760.time>=din2 & total_760.time<=dfin1);
timetotal_760 = total_760.time(1, total_760.time>=din2 & ...
    total_760.time<=dfin1);

% geostrophic
ugeosbci_760 = geosbci_760.u(:, :, geosbci_760.time>=din1 & ...
    geosbci_760.time<=dfin1);
timegeosbci_760 = geosbci_760.time(1, geosbci_760.time>=din1 & ...
    geosbci_760.time<=dfin1);
ugeosbci_760_2 = geosbci_760.u(:, :, geosbci_760.time>=din2 & ...
    geosbci_760.time<=dfin1);
timegeosbci_760_2 = geosbci_760.time(1, geosbci_760.time>=din2 & ...
    geosbci_760.time<=dfin1);

% reference
ugeosref_760_2 = geosref_760.u(:, :, geosref_760.time>=din2 & ...
    geosref_760.time<=dfin1);
timegeosref_760_2 = geosbci_760.time(1, geosref_760.time>=din2 & ...
    geosref_760.time<=dfin1);


% Finds coincident xbt and adcp transects
% os38nb 
% [x y] = meshgrid(timegeosbci_760, timetotal_760);
% [t1 t2] = meshgrid(1:length(timegeosbci_760), 1:length(timetotal_760));
% ind_tot_a = t2(abs(x-y)<=4);

% Loads coincident transects list
listcoin = load('./Datasets/coinc_xbt_adcp_trans.txt');
timecoin = datenum(listcoin);

rr = 0;
for i = 1:length(timecoin);
    iix = find(abs(timecoin(i) - timetotal_760)<2);

    if ~isempty(iix)
        rr = rr + 1;
        ind_tot_a(rr) = iix;
    end
end

% Estimating trends
[trend_tot confint_tot sigtr_tot trend_np_tot sigtr_np_tot] = ...
    calctrends2D(utotal_760, [timetotal_760-datein]/365.25, alpha, season_opt);

[trend_geo760 confint_geo760 sigtr_geo760 trend_np_geo760 sigtr_np_geo760] = ...
    calctrends2D(ugeosbci_760_2, [timegeosbci_760_2-datein]/365.25, alpha, season_opt);

[trend_ref760 confint_ref760 sigtr_ref760 trend_np_ref760 sigtr_np_ref760] = ...
    calctrends2D(ugeosref_760_2, [timegeosref_760_2-datein]/365.25, alpha, season_opt);


% ADCP at XBT times
[trend_tot_xbt confint_tot_xbt sigtr_tot_xbt trend_np_tot_xbt ...
    sigtr_np_tot_xbt] = calctrends2D(utotal_760(:, :, ind_tot_a), ...
    [timetotal_760(1, ind_tot_a)-datein]/365.25, alpha, season_opt);


fig1 = figure('color', 'w');
figure_width = 8.8;
figure_height = 8;
figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
set(fig1,'Visible', figuresVisible)
set(fig1, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
set(fig1, 'PaperPositionMode', 'auto');
ax = tight_subplot(3, 1, [0.05 0.05], [0.25 0.05], [0.15 0.18]);


axes(ax(1));
ax1Pos = ax(1).Position;
contourf(total_760.dist(1:34), -total_760.z, ...
    trend_np_tot(:, 1:34)*1e2, cont, 'edgecolor', 'none');
caxis(ax(1), [cont(1) cont(end)]); colormap(ax(1), cmap); hold on;
[aa bb] = contour(total_760.dist(1:34), -total_760.z, ...
    trend_np_tot(:, 1:34)*1e2, [0 0], 'color', 'k', 'linewidth', 2);
clabel(aa, bb, 'labelspacing', 500, 'fontsize', 7, ...
    'fontname', 'SansSerif');
axis(ax(1), [0 837.50 -760 0]);
hat = trend_np_tot*1e2;
hat(sigtr_np_tot==0) = 10;
[~, gh] = contourf(total_760.dist(1:34), -total_760.z, ...
    hat(:, 1:34), [10 10]);
set(gh, 'linestyle', 'none', 'Tag', 'HatchingRegion');
hp = findobj(ax(1), 'Tag', 'HatchingRegion');
hh = hatchfill2(hp, 'cross', 'HatchAngle', 45, 'HatchDensity', 50, ...
    'HatchColor', [0.5 0.5 0.5], 'HatchLineWidth', 0.5, 'Fill', 'off');
[cc hh] = contour(total_760.dist(1:34), -total_760.z, ...
    nanmean(utotal_760(:, 1:34, :), 3)*1e2, [10 15 20 25], 'linecolor', rgb('Cyan'), ...
    'linewidth', 0.5);
clabel(cc, hh, 'labelspacing', 1000, 'fontsize', 7, 'color', rgb('Cyan'), ...
    'fontname', 'SansSerif');
set(ax(1), 'XTickLabel', [], 'TickLength', [0.020 0.020], ...
    'TickDir', 'out', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'YTick', [-800:200:0]', 'YTickLabel', num2str([-800:200:0]'));
ylabel('Depth [m]', 'fontsize', 7, 'fontname', 'SansSerif');
p1 = plot([820 745], [0 0], '-', 'linewidth', 4, 'color', 'm');
p2 = plot([300 500], [0 0], '-', 'linewidth', 4, 'color', rgb('DarkGray'));
p3 = plot([35 135], [0 0], '-', 'linewidth', 4, 'color', rgb('Goldenrod'));
text(420, -640, ['Total velocity u_{tot} '], 'fontsize', 7,...
    'fontname', 'SansSerif', 'fontweight', 'bold');


axes(ax(2));
ax2Pos = ax(2).Position;
contourf(total_760.dist(1:34), -total_760.z, ...
    trend_np_tot_xbt(:, 1:34)*1e2, cont, 'edgecolor', 'none');
caxis(ax(2), [cont(1) cont(end)]); colormap(ax(2), cmap); hold on;
[aa bb] = contour(total_760.dist(1:34), -total_760.z, ...
    trend_np_tot_xbt(:, 1:34)*1e2, [0 0], 'color', 'k', 'linewidth', 2);
clabel(aa, bb, 'labelspacing', 500, 'fontsize', 7, 'fontname', 'SansSerif');
axis(ax(2), [0 837.50 -760 0]);
hat2 = trend_np_tot_xbt*1e2;
hat2(sigtr_np_tot_xbt==0) = 10;
[~, gh2] = contourf(total_760.dist(1:34), -total_760.z, ...
    hat2(:, 1:34), [10 10]);
set(gh2, 'linestyle', 'none', 'Tag', 'HatchingRegion');
hp2 = findobj(ax(2), 'Tag', 'HatchingRegion');
hh2 = hatchfill2(hp2, 'cross', 'HatchAngle', 45, 'HatchDensity', 50, ...
    'HatchColor', [0.5 0.5 0.5], 'HatchLineWidth', 0.5, 'Fill', 'off');
[cc hh] = contour(total_760.dist(1:34), -total_760.z, ...
    nanmean(utotal_760(:, 1:34, ind_tot_a), 3)*1e2, [10 15 20 25], ...
    'linecolor', rgb('Cyan'), 'linewidth', 0.5);
clabel(cc, hh, 'labelspacing', 1000, 'fontsize', 7, 'color', rgb('Cyan'),...
    'fontname', 'SansSerif');
set(ax(2), 'XTickLabel', [], 'TickLength', [0.020 0.020], ...
    'TickDir', 'out', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'YTick', [-800:200:0]', 'YTickLabel', num2str([-800:200:0]'));
ylabel('Depth [m]', 'fontsize', 7, 'fontname', 'SansSerif');
p1 = plot([820 745], [0 0], '-', 'linewidth', 4, 'color', 'm');
p2 = plot([300 500], [0 0], '-', 'linewidth', 4, 'color', rgb('DarkGray'));
p3 = plot([35 135], [0 0], '-', 'linewidth', 4, 'color', rgb('Goldenrod'));
text(420, -640, ['Total velocity^{\ast} u_{tot}'], 'fontsize', 7,...
    'fontname', 'SansSerif', 'fontweight', 'bold');


axes(ax(3));
ax3Pos = ax(3).Position;
contourf(geosbci_760.dist_u(1:76, 1), -geosbci_760.z, ...
    trend_np_geo760(:, 1:76)*1e2, cont, 'edgecolor', 'none');
caxis(ax(3), [cont(1) cont(end)]); colormap(ax(3), cmap); hold on;
[aa bb] = contour(geosbci_760.dist_u(1:76, 1), -geosbci_760.z, ...
    trend_np_geo760(:, 1:76)*1e2, [0 0], 'color', 'k', 'linewidth', 2);
clabel(aa, bb, 'labelspacing', 500, 'fontsize', 7, 'fontname', 'SansSerif');
ylabel('Depth [m]', 'fontsize', 7, 'fontname', 'SansSerif');
set(ax(3), 'TickLength', [0.020 0.020], 'TickDir', 'out', ...
    'fontsize', 7, 'fontname', 'SansSerif', ...
    'YTick', [-800:200:0]', 'YTickLabel', num2str([-800:200:0]'));
col = colorbar(ax(3), 'Location', 'South', 'Ticks', cont(2:2:end)',...
    'TickLabels', num2str(cont(2:2:end)'), 'fontsize', 7, ...
    'fontname', 'SansSerif', 'TickDirection', 'out', 'TickLength', 0.015, ...
    'AxisLocation', 'out');
col.Position(1) = ax3Pos(1);
col.Position(2) = 0.16;
col.Position(3) = ax3Pos(3);
col.Position(4) = 0.015;
cbarrow(col, cmap);
ylabel(col, 'Velocity trend [cm s^{-1} year^{-1}]', 'fontsize', 7, ...
    'fontname', 'SansSerif');
hat4 = trend_np_geo760*1e2;
hat4(sigtr_np_geo760==0) = 10;
[~, gh4] = contourf(geosbci_760.dist_u(1:76, 1), -geosbci_760.z, ...
    hat4(:, 1:76), [10 10]);
set(gh4, 'linestyle', 'none', 'Tag', 'HatchingRegion');
hp4 = findobj(ax(3), 'Tag', 'HatchingRegion');
hh4 = hatchfill2(hp4, 'cross', 'HatchAngle', 45, 'HatchDensity', 50, ...
    'HatchColor', [0.5 0.5 0.5], 'HatchLineWidth', 0.5, 'Fill', 'off');
[cc hh] = contour(geosbci_760.dist_u(1:70, 1), -geosbci_760.z, ...
    nanmean(ugeosbci_760(:, 1:70, :), 3)*1e2, [5 8 10], ...
    'linecolor', rgb('Cyan'), 'linewidth', 0.5);
clabel(cc, hh, 'fontsize', 7, 'color', rgb('Cyan'), 'fontname', 'SansSerif');
axis(ax(3), [0 837.50 -760 0]);
p1 = plot([820 745], [0 0], '-', 'linewidth', 4, 'color', 'm');
p2 = plot([300 500], [0 0], '-', 'linewidth', 4, 'color', rgb('DarkGray'));
p3 = plot([35 135], [0 0], '-', 'linewidth', 4, 'color', rgb('Goldenrod'));
text(220, -640, ['Geostrophic velocity u_{geo}'], 'fontsize', 7, ...
    'fontname', 'SansSerif', 'fontweight', 'bold');
xlabel('Distance from North [km]', 'fontsize', 7, 'fontname', 'SansSerif')


annotation(fig1, 'textbox', ...
    [ax1Pos(1)+0.01 sum(ax1Pos([2, 4]))+0.025 0.0480 0.0353],...
    'String', 'SAF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', ...
    [ax1Pos(1)+ax1Pos(3)/2-0.05 sum(ax1Pos([2, 4]))+0.025 0.0480 0.0353],...
    'String', 'PF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', [...
    sum(ax1Pos([1, 3]))-0.12 sum(ax1Pos([2, 4]))+0.025 0.0480 0.0353],...
    'String', 'SACCF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', ...
    [ax2Pos(1)+0.01 sum(ax2Pos([2, 4]))+0.025 0.0480 0.0353],...
    'String', 'SAF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', ...
    [ax2Pos(1)+ax2Pos(3)/2-0.05 sum(ax2Pos([2, 4]))+0.025 0.0480 0.0353],...
    'String', 'PF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', ...
    [sum(ax2Pos([1, 3]))-0.12 sum(ax2Pos([2, 4]))+0.025 0.0480 0.0353],...
    'String', 'SACCF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', ...
    [ax3Pos(1)+0.01 sum(ax3Pos([2, 4]))+0.025 0.0480 0.0353],...
    'String', 'SAF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', ...
    [ax3Pos(1)+ax3Pos(3)/2-0.05 sum(ax3Pos([2, 4]))+0.025 0.0480 0.0353],...
    'String', 'PF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', ...
    [sum(ax3Pos([1, 3]))-0.12 sum(ax3Pos([2, 4]))+0.025 0.0480 0.0353],...
    'String', 'SACCF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', [0 sum(ax(1).Position([2,4]))+0.025 ...
    0.0480 0.0353], 'String', 'a', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig1, 'textbox', [0 sum(ax(2).Position([2,4]))+0.025 ...
    0.0480 0.0353], 'String', 'b', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig1, 'textbox', [0 sum(ax(3).Position([2,4]))+0.025 ...
    0.0480 0.0353], 'String', 'c', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');

print('-dpdf', './FiguresPaper/fig3_DPvel_trend.pdf', '-r500');
