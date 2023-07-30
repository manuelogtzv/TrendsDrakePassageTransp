% Code for plotting Figure SI - Trends in the cross-transect velocity as
% a function of depth and distance across Passage. Trends are estimated
% using the least-squares fit after removing the seasonal cycle.
% Significance tests are done using the modified Mann-Kendall test. Data
% employed is the geostrohpic velocity u_{geo} and total u_{tot} using all 
% available transects. Trends are also estimated for u_{tot} using only
% those transect that have a coincident XBT/XCTD transect.
%
%
% Dec 10/2022 - Created by Manuel. O. Gutierrez-Villanueva

clear all;
close all;


% Define variables
datein = datenum(1996, 01, 01, 0, 0 ,0);
alpha = 0.05;
gsize = 25;

cont = [-1.8:0.2:1.8];%vel contours
cmap = getPyPlot_cMap('seismic', length(cont)-1);

din1 = datenum(1998, 1, 1);
din2 = datenum(2005, 10, 1);
dfin1 = datenum(2019, 4, 30);

season_opt = 1;
op_mostreptran = 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loads os38 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['lmgvelos38nb_50_300.mat']);

% Finds south and northbound transects
indsouth = find(lmgdata.northbound == 1);
indnorth = find(lmgdata.northbound == 0);

% Adds the uoffs
load('trandsmiss_theta.mat');
uoffs = mean(abs([uoffs_north, uoffs_south]));
lmgdata.u(:, :, indsouth) = lmgdata.u(:, :, indsouth) - uoffs;
lmgdata.u(:, :, indnorth) = lmgdata.u(:, :, indnorth) + uoffs;
nomz = 970;
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
    
dist = [lmgdata.dist(1):gsize:size(lmgdata.dist,1)*gsize ];
mask = lmgdata.dist./lmgdata.dist;

% Option for using the most repeated transects
if op_mostreptran == 1;
    fprintf('\n\n total number of transects: %s\n',...
        num2str(length(lmgdata.time)));

    box1 = [ -64.606712283408200 -56.324379209501110;
             -64.534099360692622 -56.568973205588243;
             -64.373586584163448 -56.541796094911895;
             -64.438556041330017 -56.351556320177458];

    box2 = [-63.009254066180596 -61.670918367346943;
            -62.908300616937744 -61.915816326530610;
            -62.504486819966345 -61.752551020408163;
            -62.625630959057766 -61.426020408163268];

    ii1 = inpolygon(lmgdata.lon-360, lmgdata.lat, box1(:, 1), ...
        box1(:, 2));

    [~, yj] = find(ii1 == 1);
    yj = unique(yj);

    mm1 = inpolygon(lmgdata.lon(:, yj)-360, lmgdata.lat(:, yj), ...
        box2(:, 1), box2(:, 2));

    [~, zj] = find(mm1 == 1);
    zj = unique(zj);

    % HOw many transects along the most repeated line
    fprintf('\n\n total number of transects (most repeated): %s\n',...
        num2str(length(lmgdata.time(yj(zj)))));
end

os38_970 = lmgdata;
os38_970.dist = dist; clear lmgdata uoffs* indn* inds* i ii theta* dist;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% GEOSTROPHIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['geosbarocDP_noaver.mat']);
geosbaroc.u = geosbaroc.u(:, 1:end-1, :);
distubci = 0.5*(geosbaroc.dist(1:end-1, :) + ...
            geosbaroc.dist(2:end, :));
distubci = distubci - distubci(1, :);

geosbci_780 = geosbaroc;
geosbci_780.distu = distubci + 12.5; clear geosbaroc 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Removes spurious data due to bathymetry mask
os38_970.u(os38_970.z>694, 1, :) = NaN;
geosbci_780.u(geosbci_780.z>694, 1, :) = NaN;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Trends velocity  %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Renames variables for simplicity
% os38
uos38_970 = os38_970.u(:, :, os38_970.time>=din2 & os38_970.time<=dfin1);
timeos38_970 = os38_970.time(1, os38_970.time>=din2 & ...
    os38_970.time<=dfin1);


% geostrophic
ugeosbci_760 = geosbci_780.u(:, :, geosbci_780.time>=din1 & ...
    geosbci_780.time<=dfin1);
timegeosbci_760 = geosbci_780.time(1, geosbci_780.time>=din1 & ...
    geosbci_780.time<=dfin1);

% Finds coincident xbt and adcp transects
% os38nb 
[x y] = meshgrid(timegeosbci_760, timeos38_970);
[t1 t2] = meshgrid(1:length(timegeosbci_760), 1:length(timeos38_970));
ind_os38_a = t2(abs(x-y)<=4);

% Estimating trends
[trend_os38_mrl confint_os38_mrl sigtr_os38_mrl trend_np_os38_mrl sigtr_np_os38_mrl] = ...
    calctrends2D(os38_970.u(:, :, yj(zj)), ...
    [os38_970.time(yj(zj))-datein]/365.25, alpha, season_opt);

[trend_geo760 confint_geo760 sigtr_geo760 trend_np_geo760 sigtr_np_geo760] = ...
    calctrends2D(ugeosbci_760, [timegeosbci_760-datein]/365.25, alpha, season_opt);

% ADCP at XBT times
[trend_os38_xbt confint_os38_xbt sigtr_os38_xbt trend_np_os38_xbt ...
    sigtr_np_os38_xbt] = calctrends2D(uos38_970(:, :, ind_os38_a), ...
    [timeos38_970(1, ind_os38_a)-datein]/365.25, alpha, season_opt);


fig1 = figure('color', 'w');
figure_width = 8.8;
figure_height = 12;
figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
set(fig1,'Visible', figuresVisible)
set(fig1, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
set(fig1, 'PaperPositionMode', 'auto');
ax = tight_subplot(3, 1, [0.08 0.08], [0.25 0.05], [0.15 0.18]);


axes(ax(1));
ax1Pos = ax(1).Position;
contourf(os38_970.dist(1:34), -os38_970.z, ...
    trend_np_os38_mrl(:, 1:34)*1e2, cont, 'edgecolor', 'none');
caxis(ax(1), [cont(1) cont(end)]); colormap(ax(1), cmap); hold on;
[aa bb] = contour(os38_970.dist(1:34), -os38_970.z, ...
    trend_np_os38_mrl(:, 1:34)*1e2, [0 0], 'color', 'k', 'linewidth', 2);
clabel(aa, bb, 'labelspacing', 500, 'fontsize', 7, ...
    'fontname', 'SansSerif');
axis(ax(1), [0 837.50 -760 0]);
hat = trend_np_os38_mrl*1e2;
hat(sigtr_np_os38_mrl==0) = 10;
[~, gh] = contourf(os38_970.dist(1:34), -os38_970.z, ...
    hat(:, 1:34), [10 10]);
set(gh, 'linestyle', 'none', 'Tag', 'HatchingRegion');
hp = findobj(ax(1), 'Tag', 'HatchingRegion');
hh = hatchfill2(hp, 'cross', 'HatchAngle', 45, 'HatchDensity', 50, ...
    'HatchColor', [0.5 0.5 0.5], 'HatchLineWidth', 0.5, 'Fill', 'off');
[cc hh] = contour(os38_970.dist(1:34), -os38_970.z, ...
    nanmean(uos38_970(:, 1:34, :), 3)*1e2, [10 15 20 25], 'linecolor', rgb('Cyan'), ...
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
text(420, -640, ['Total velocity^{a} u_{tot}'], 'fontsize', 7,...
    'fontname', 'SansSerif', 'fontweight', 'bold');


axes(ax(2));
ax2Pos = ax(2).Position;
contourf(os38_970.dist(1:34), -os38_970.z, ...
    trend_np_os38_xbt(:, 1:34)*1e2, cont, 'edgecolor', 'none');
caxis(ax(2), [cont(1) cont(end)]); colormap(ax(2), cmap); hold on;
[aa bb] = contour(os38_970.dist(1:34), -os38_970.z, ...
    trend_np_os38_xbt(:, 1:34)*1e2, [0 0], 'color', 'k', 'linewidth', 2);
clabel(aa, bb, 'labelspacing', 500, 'fontsize', 7, 'fontname', 'SansSerif');
axis(ax(2), [0 837.50 -760 0]);
hat2 = trend_np_os38_xbt*1e2;
hat2(sigtr_np_os38_xbt==0) = 10;
[~, gh2] = contourf(os38_970.dist(1:34), -os38_970.z, ...
    hat2(:, 1:34), [10 10]);
set(gh2, 'linestyle', 'none', 'Tag', 'HatchingRegion');
hp2 = findobj(ax(2), 'Tag', 'HatchingRegion');
hh2 = hatchfill2(hp2, 'cross', 'HatchAngle', 45, 'HatchDensity', 50, ...
    'HatchColor', [0.5 0.5 0.5], 'HatchLineWidth', 0.5, 'Fill', 'off');
[cc hh] = contour(os38_970.dist(1:34), -os38_970.z, ...
    nanmean(uos38_970(:, 1:34, ind_os38_a), 3)*1e2, [10 15 20 25], ...
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
text(420, -640, ['Total velocity^{b} u_{tot}'], 'fontsize', 7,...
    'fontname', 'SansSerif', 'fontweight', 'bold');


axes(ax(3));
ax3Pos = ax(3).Position;
contourf(geosbci_780.distu(1:76, 1), -geosbci_780.z, ...
    trend_np_geo760(:, 1:76)*1e2, cont, 'edgecolor', 'none');
caxis(ax(3), [cont(1) cont(end)]); colormap(ax(3), cmap); hold on;
[aa bb] = contour(geosbci_780.distu(1:76, 1), -geosbci_780.z, ...
    trend_np_geo760(:, 1:76)*1e2, [0 0], 'color', 'k', 'linewidth', 2);
clabel(aa, bb, 'labelspacing', 500, 'fontsize', 7, 'fontname', 'SansSerif');
ylabel('Depth [m]', 'fontsize', 7, 'fontname', 'SansSerif');
set(ax(3), 'TickLength', [0.020 0.020], 'TickDir', 'out', ...
    'fontsize', 7, 'fontname', 'SansSerif', ...
    'XTick', ax(1).XTick, 'XTickLabel', num2str(ax(1).XTick(:)), ...
    'YTick', [-800:200:0]', 'YTickLabel', num2str([-800:200:0]'));
col = colorbar(ax(3), 'Location', 'South', 'Ticks', cont(2:2:end)',...
    'TickLabels', num2str(cont(2:2:end)'), 'fontsize', 7, ...
    'fontname', 'SansSerif', 'TickDirection', 'out', 'TickLength', 0.015, ...
    'AxisLocation', 'out');
col.Position(1) = ax3Pos(1);
col.Position(2) = 0.18;
col.Position(3) = ax3Pos(3);
col.Position(4) = 0.02;
cbarrow(col, cmap);
ylabel(col, 'Velocity trend [cm s^{-1} year^{-1}]', 'fontsize', 7, ...
    'fontname', 'SansSerif');
hat4 = trend_np_geo760*1e2;
hat4(sigtr_np_geo760==0) = 10;
[~, gh4] = contourf(geosbci_780.distu(1:76, 1), -geosbci_780.z, ...
    hat4(:, 1:76), [10 10]);
set(gh4, 'linestyle', 'none', 'Tag', 'HatchingRegion');
hp4 = findobj(ax(3), 'Tag', 'HatchingRegion');
hh4 = hatchfill2(hp4, 'cross', 'HatchAngle', 45, 'HatchDensity', 50, ...
    'HatchColor', [0.5 0.5 0.5], 'HatchLineWidth', 0.5, 'Fill', 'off');
[cc hh] = contour(geosbci_780.distu(1:70, 1), -geosbci_780.z, ...
    nanmean(ugeosbci_760(:, 1:70, :), 3)*1e2, [5 8 10], ...
    'linecolor', rgb('Cyan'), 'linewidth', 0.5);
clabel(cc, hh, 'fontsize', 7, 'color', rgb('Cyan'), 'fontname', 'SansSerif');
axis(ax(3), [0 837.50 -760 0]);
p1 = plot([820 745], [0 0], '-', 'linewidth', 4, 'color', 'm');
p2 = plot([300 500], [0 0], '-', 'linewidth', 4, 'color', rgb('DarkGray'));
p3 = plot([35 95], [0 0], '-', 'linewidth', 4, 'color', rgb('Goldenrod'));
text(220, -640, ['Geostrophic velocity u_{geo}'], 'fontsize', 7, ...
    'fontname', 'SansSerif', 'fontweight', 'bold');
xlabel('Distance from North [km]', 'fontsize', 7, 'fontname', 'SansSerif')

annotation(fig1, 'textbox', ...
    [ax1Pos(1)+0.01 sum(ax1Pos([2, 4]))+0.01 0.0480 0.0353],...
    'String', 'SAF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', ...
    [ax1Pos(1)+ax1Pos(3)/2-0.05 sum(ax1Pos([2, 4]))+0.01 0.0480 0.0353],...
    'String', 'PF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', [...
    sum(ax1Pos([1, 3]))-0.12 sum(ax1Pos([2, 4]))+0.01 0.0480 0.0353],...
    'String', 'SACCF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', ...
    [ax2Pos(1)+0.01 sum(ax2Pos([2, 4]))+0.01 0.0480 0.0353],...
    'String', 'SAF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', ...
    [ax2Pos(1)+ax2Pos(3)/2-0.05 sum(ax2Pos([2, 4]))+0.01 0.0480 0.0353],...
    'String', 'PF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', ...
    [sum(ax2Pos([1, 3]))-0.12 sum(ax2Pos([2, 4]))+0.01 0.0480 0.0353],...
    'String', 'SACCF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', ...
    [ax3Pos(1)+0.01 sum(ax3Pos([2, 4]))+0.01 0.0480 0.0353],...
    'String', 'SAF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', ...
    [ax3Pos(1)+ax3Pos(3)/2-0.05 sum(ax3Pos([2, 4]))+0.01 0.0480 0.0353],...
    'String', 'PF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', ...
    [sum(ax3Pos([1, 3]))-0.12 sum(ax3Pos([2, 4]))+0.01 0.0480 0.0353],...
    'String', 'SACCF', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none');
annotation(fig1, 'textbox', [0 sum(ax(1).Position([2,4]))+0.03 ...
    0.0480 0.0353], 'String', 'a', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig1, 'textbox', [0 sum(ax(2).Position([2,4]))+0.03 ...
    0.0480 0.0353], 'String', 'b', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig1, 'textbox', [0 sum(ax(3).Position([2,4]))+0.03 ...
    0.0480 0.0353], 'String', 'c', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
% print('-dpdf', './FiguresPaper/fig3_DPvel_trend.pdf', '-r500');
