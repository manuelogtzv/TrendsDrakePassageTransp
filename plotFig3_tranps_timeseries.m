% Code for estimating time series of transport per distance bin for the
% geostrophic, total and reference components. 
%
% 2023/06/25 - Created by Manuel O. Gutierrez-Villanueva

clear all;
close all;

pycmd = '/Users/manuelgutierrez/anaconda3/bin/python';


% Define variables
datein = datenum(1996, 01, 01, 0, 0 ,0);
alpha = 0.05;
gsize = 25;
op_season = 1;

alpha = 0.05;
ti = datenum(2005, 10, 1);
tf = datenum(2019, 4, 31);
colticks1 = [-16:2:16];
colticks2 = [-3:0.5:3];
colticks3 = [-6:1:6];
cmap1 = getPyPlot_cMap('seismic', length(colticks1)-1, [], pycmd);
cmap2 = getPyPlot_cMap('seismic', length(colticks2)-1, [], pycmd);
cmap3 = getPyPlot_cMap('seismic', length(colticks3)-1, [], pycmd);
timeticks = datenum(2004:2:2020, 1, 1);

din1 = datenum(1999, 11, 1);
din2 = datenum(2005, 10, 1);
dfin1 = datenum(2019, 4, 30);
saccf = [820 745 745 820];
pf = [300 500 500 300];
saf = [35 135 135 35];
distticks = [0:200:800]';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loads os38 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads transport
load('TotalTransport_os38nbmaxz_760_2.mat'); %os38nb
tot = totaltransp;
tot.transp = tot.transp(:, totaltransp.time>=din2 & totaltransp.time<=dfin1);
tot.time = tot.time(totaltransp.time>=din2 & totaltransp.time<=dfin1);
tot.transp(tot.dist>840, :) = NaN;
tot.dist(tot.dist>840) = NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% GEOSTROPHIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['geostrophic_noaver.mat']);

geos = geosbaroc;
geos.transp = geos.transp(:, geosbaroc.time>=din2 & geosbaroc.time<=dfin1);
geos.dist_u = geos.dist_u(:, geosbaroc.time>=din2 & geosbaroc.time<=dfin1);
geos.time = geos.time(geosbaroc.time>=din2 & geos.time<=dfin1);
geos.transp(nanmean(geos.dist_u, 2)>840, :) = NaN;
geos.dist_u(nanmean(geos.dist_u, 2)>840, :) = NaN;
% geosbci.distu = distubci; 
clear geosbaroc 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% REFERENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['reference_os38nb_760.mat' ]);
ref = geosref;
ref.transp = ref.transp(:, geosref.time>=din2 & geosref.time<=dfin1);
ref.dist = ref.dist(:, geosref.time>=din2 & geosref.time<=dfin1);
ref.time = ref.time(geosref.time>=din2 & geosref.time<=dfin1);
ref.transp(nanmean(ref.dist, 2)>820, :) = NaN;
ref.dist(nanmean(ref.dist, 2)>820, :) = NaN;
clear geosbaroc geosref totaltransp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Sum of reference and geostrophic
[a b] = meshgrid(geos.time, ref.time);
[ta1 ta2] = meshgrid(1:length(geos.time), 1:length(ref.time));
ind_gbci = ta1(abs(a-b)<=1);
ind_ref = ta2(abs(a-b)<=1);

totsum = [ref.transp(:, ind_ref) + geos.transp(:, ind_gbci)]*1e-6;

% Gets only those total transp transects with a coincident xbt transect
[c d] = meshgrid(geos.time, tot.time);
[tc1 tc2] = meshgrid(1:length(geos.time), 1:length(tot.time));
ind_tot2 = tc2(abs(c-d)<=4);
ind_geo2 = tc1(abs(c-d)<=1);

totxbt = tot.transp(:, ind_tot2)*1e-6; 

tot_anom = tot.transp - nanmean(tot.transp, 2);
geos_anom = geos.transp - nanmean(geos.transp, 2);
ref_anom = ref.transp - nanmean(ref.transp, 2);

% Meshgrids of time and distance
[timetot, disttot] = meshgrid(tot.time, tot.dist);
[timegeo, ~] = meshgrid(geos.time, geos.dist_u(:, 1));
[timeref, ~] = meshgrid(ref.time, ref.dist(:, 1));

% Calculates trends
[tr_tot, ci_tot, ~, tr_np_tot sig_np_tot] = calctrends(tot.transp, ...
    [tot.time-datenum(1996, 1, 1)]/365.25, alpha, op_season);
[tr_tot2, ci_tot2, ~, tr_np_tot2 sig_np_tot2] = calctrends(tot.transp(:, ind_tot2), ...
    [tot.time(1, ind_tot2)-datenum(1996, 1, 1)]/365.25, alpha, op_season);
[tr_geo, ci_geo, ~, tr_np_geo sig_np_geo] = calctrends(geos.transp, ...
    [geos.time-datenum(1996, 1, 1)]/365.25, alpha, op_season);
[tr_ref, ci_ref, ~, tr_np_ref sig_np_ref] = calctrends(ref.transp, ...
    [ref.time-datenum(1996, 1, 1)]/365.25, alpha, op_season);

% Figure: Plots anomaly time series of transport in across-Passage distance
fig0 = figure('color', 'w');
figure_width = 18;
figure_height = 10;
figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
set(fig0,'Visible', figuresVisible)
set(fig0, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
set(fig0, 'PaperPositionMode', 'auto');


ha = tight_subplot(3, 5, [0.03 0.03], [0.15 0.05], [0.11 0.11]);
dx = sum(ha(3).Position([1, 3]));
dx = dx - ha(1).Position(1);
ha(1).Position(3) = dx;
ha(6).Position(3) = dx;
ha(11).Position(3) = dx;
ha(4).Position(3) = 0.15;
ha(9).Position(3) = 0.15;
ha(14).Position(3) = 0.15;
ha(5).Position(3) = 0.15;
ha(10).Position(3) = 0.15;
ha(15).Position(3) = 0.15;
ha(5).Position(1) = ha(5).Position(1)+0.01;
ha(10).Position(1) = ha(10).Position(1)+0.01;
ha(15).Position(1) = ha(15).Position(1)+0.01

ha(2).Visible = 'off';
ha(3).Visible = 'off';
ha(7).Visible = 'off';
ha(8).Visible = 'off';
ha(12).Visible = 'off';
ha(13).Visible = 'off';


axes(ha(1));
pos1 = ha(1).Position;
p1 = patch([ones(1, 2)*datenum(1995, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    saf, rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2); hold on;
p2 = patch([ones(1, 2)*datenum(1995, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    pf, rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2);
p3 = patch([ones(1, 2)*datenum(1995, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    saccf, rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2);
scatter(timetot(:), disttot(:), 10, tot_anom(:)*1e-6, 's', 'filled')
datetick('x');
set(ha(1), 'box', 'on', 'TickLength', [0.01 0.01], 'YTick',  distticks,  ...
    'YTickLabel', num2str(distticks, '%3.0f'), 'TickDir', 'out', ...
    'XTick', datenum(1996:2020, 1 ,1)', 'XTickLabel', [], ...
    'TickLength', [0.01 0.01], 'YDir', 'reverse');
ylabel('Distance from North [km]');
text(datenum(2005, 6, 1), 100, 'Total', 'fontweight', 'bold');
axis(ha(1), [datenum(2005, 1, 2) datenum(2020, 4, 1) -20 870]);
ha(1).Position = pos1;
colormap(ha(1), cmap1);
caxis(ha(1), [min(colticks1) max(colticks1)]);



axes(ha(4));
pos2 = ha(4).Position;
% plot(nanmean(tot.transp*1e-6), tot.psic, 'r', 'linewidth', 1); hold on
[hl1 hp1] = boundedline(nanmean(tot.transp*1e-6, 2), tot.dist, ...
    nanstd(tot.transp, 0 , 2)*1e-6./sqrt(sum(~isnan(tot.transp), 2)), '-r',...
    'alpha', 'Orientation', 'horiz');
hold on;
[hl2 hp2] = boundedline(nanmean(totxbt, 2), tot.dist, ...
    nanstd(totxbt, 0 , 2)./sqrt(sum(~isnan(totxbt), 2)), '-',...
    'alpha', 'Orientation', 'horiz');
hl2.MarkerSize = 5;
hl2.Color = rgb('Green');
hp2.FaceColor = rgb('Green');
% plot(nanmean(totxbt, 1), tot.psic, 'sr', 'markersize', 3);
axis(ha(4), [-1 6.3 -20 870]); lim = axis;
p1 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], saf, ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2); 
p2 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], pf, ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
p3 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], saccf, ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
text(1, mean(saf), 'SAF', 'HorizontalAlignment', 'center');
text(1, mean(pf), 'PF', 'HorizontalAlignment', 'center');
text(4.5, mean(saccf), 'SACCF', 'HorizontalAlignment', 'center');
text(3, 200, 'Total', 'fontweight', 'bold');
set(ha(4), 'YTick', ha(1).YTick, 'YTickLabel', [], 'TickDir', 'out', ...
    'XTick', [0:1:6]', 'XTickLabel', [], 'TickLength', [0.03 0.03], ...
    'box', 'on', 'YDir', 'reverse');
plot([0 0], [-20 870], '-k', 'linewidth', 0.5);
grid on;
leg = legend([hl1; hl2], [num2str(nansum(nanmean(tot.transp*1e-6, 2)),'%2.2f')...
    char(177) num2str(nanstd(nansum(tot.transp*1e-6, 1))/sqrt(size(tot.transp, 2)),'%2.2f') ...
    ' Sv'], [num2str(nansum(nanmean(totxbt, 2)),'%2.2f')...
    char(177) num2str(nanstd(nansum(totxbt, 1))/sqrt(size(totxbt, 2)),'%2.2f') ...
    ' Sv'], 'fontsize', 5, 'fontname', 'SansSerif');

axes(ha(5));
pos3 = ha(5).Position;
% plot(0, 10, 'r', 'linewidth', 1); hold on
[hl1 hp1] = boundedline(tr_tot*1e-6, tot.dist, ci_tot*1e-6, '-r',...
    'alpha', 'Orientation', 'horiz');
hold on;
plot(tr_tot(sig_np_tot==1)*1e-6, tot.dist(sig_np_tot==1), 'sr', ...
    'markerfacecolor', 'r', 'markersize', 5);
hl2.MarkerSize = 10;
[hl2 hp2] = boundedline(tr_tot2*1e-6, tot.dist, ci_tot*1e-6, '-r',...
    'alpha', 'Orientation', 'horiz');
pp2 = plot(tr_tot2(sig_np_tot2==1)*1e-6, tot.dist(sig_np_tot2==1), 's', ...
    'markerfacecolor', rgb('Green'), 'markeredgecolor', rgb('Green'), ...
    'markersize', 5);
hl2.Color = pp2.MarkerFaceColor;
hp2.FaceColor = pp2.MarkerFaceColor;
% plot(nanmean(totxbt, 1), tot.psic, 'sr', 'markersize', 3);
axis(ha(5), [-0.35 0.35 -20 870]); 
lim = axis;
p1 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], saf, ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2); 
p2 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], pf, ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
p3 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], saccf, ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
text(0.20, mean(saf), 'SAF', 'HorizontalAlignment', 'center');
text(0.20, mean(pf), 'PF', 'HorizontalAlignment', 'center');
text(0.20, mean(saccf), 'SACCF', 'HorizontalAlignment', 'center');
text(-0.30, 200, 'Total', 'fontweight', 'bold');
set(ha(5), 'YTick', ha(1).YTick, 'YTickLabel', [], 'TickDir', 'out', ...
    'XTick', [-0.4:0.1:0.4]', 'XTickLabel', [], 'TickLength', [0.03 0.03], ...
    'box', 'on', 'YDir', 'reverse');
plot([0 0], [-20 870], '-k', 'linewidth', 0.5);
leg1a = legend([hl1; hl2], 'All transects', 'Subsampled', ...
    'fontsize', 5, 'fontname', 'SansSerif');
grid on;


axes(ha(6));
pos4 = ha(6).Position;
p1 = patch([ones(1, 2)*datenum(1995, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    saf, rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2); hold on;
p2 = patch([ones(1, 2)*datenum(1995, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    pf, rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2);
p3 = patch([ones(1, 2)*datenum(1995, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    saccf, rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2);
scatter(timegeo(:), geos.dist_u(:), 10, geos_anom(:)*1e-6, 's', 'filled')
datetick('x');
set(ha(6), 'box', 'on', 'TickLength', [0.01 0.01], 'YTick',  distticks,  ...
    'YTickLabel', num2str(distticks, '%3.0f'), 'TickDir', 'out', ...
    'XTick', datenum(1996:2020, 1 ,1)', 'XTickLabel', [], ...
    'TickLength', [0.01 0.01], 'YDir', 'reverse');
ylabel('Distance from North [km]');
text(datenum(2005, 6, 1), 100, 'Geostrophic', 'fontweight', 'bold');
axis(ha(6), [datenum(2005, 1, 2) datenum(2020, 4, 1) -20 870]);
ha(6).Position = pos4;
colormap(ha(6), cmap2);
caxis(ha(6), [min(colticks2) max(colticks2)]);

axes(ha(9));
pos5 = ha(9).Position;
% plot(nanmean(geos.transp*1e-6), geos.psic, 'k', 'linewidth', 1); hold on
[hl3 hp3] = boundedline(nanmean(geos.transp*1e-6, 2), nanmean(geos.dist_u, 2), ...
    nanstd(geos.transp, 0 , 2)*1e-6./sqrt(sum(~isnan(geos.transp), 2)), '-k',...
    'alpha', 'Orientation', 'horiz');
hold on;
axis(ha(9), [-1 6.3 -20 870]); lim = axis;
p1 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], saf, ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2); 
p2 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], pf, ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
p3 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], ...
    saccf, ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
text(4, mean(saf), 'SAF', 'HorizontalAlignment', 'center');
text(4, mean(pf), 'PF', 'HorizontalAlignment', 'center');
text(4.5, mean(saccf), 'SACCF', 'HorizontalAlignment', 'center');
text(1, 200, 'Geostrophic', 'fontweight', 'bold');
set(ha(9), 'YTick', ha(5).YTick, 'YTickLabel', [], 'TickDir', 'out', ...
    'XTick', [0:1:18]', 'XTickLabel', [], 'TickLength', [0.03 0.03], ...
    'box', 'on', 'YDir', 'reverse');
plot([0 0], [-20 870], '-k', 'linewidth', 0.5);
grid on; 
leg1 = legend([hl3], [num2str(nansum(nanmean(geos.transp*1e-6, 2)),'%2.2f')...
    char(177) num2str(nanstd(nansum(geos.transp*1e-6, 1))/sqrt(size(geos.transp, 2)),'%2.2f') ...
    ' Sv'], 'fontsize', 5, 'fontname', 'SansSerif');

axes(ha(10));
pos6 = ha(10).Position;
dgeo = nanmean(geos.dist_u, 2);
% plot(0, 10, 'r', 'linewidth', 1); hold on
[hl1 hp1] = boundedline(tr_geo*1e-6, dgeo, ci_geo*1e-6, '-k',...
    'alpha', 'Orientation', 'horiz');
hold on;
plot(tr_geo(sig_np_geo==1)*1e-6, dgeo(sig_np_geo==1), 'sk', ...
    'markerfacecolor', 'k', 'markersize', 5);
hl2.MarkerSize = 10;
% plot(nanmean(totxbt, 1), tot.psic, 'sr', 'markersize', 3);
axis(ha(10), [-0.35 0.35 -20 870]); 
lim = axis;
p1 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], saf, ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2); 
p2 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], pf, ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
p3 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], saccf, ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
text(0.20, mean(saf), 'SAF', 'HorizontalAlignment', 'center');
text(0.20, mean(pf), 'PF', 'HorizontalAlignment', 'center');
text(0.20, mean(saccf), 'SACCF', 'HorizontalAlignment', 'center');
text(-0.34, 200, 'Geostrophic', 'fontweight', 'bold');
set(ha(10), 'YTick', ha(1).YTick, 'YTickLabel', [], 'TickDir', 'out', ...
    'XTick', [-0.4:0.1:0.4]', 'XTickLabel', [], 'TickLength', [0.03 0.03], ...
    'box', 'on', 'YDir', 'reverse');
plot([0 0], [-20 870], '-k', 'linewidth', 0.5);
grid on;


axes(ha(11));
pos7 = ha(11).Position;
p1 = patch([ones(1, 2)*datenum(1995, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    saf, rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2); hold on;
p2 = patch([ones(1, 2)*datenum(1995, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    pf, rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2);
p3 = patch([ones(1, 2)*datenum(1995, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    saccf, rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2);
scatter(timeref(:), ref.dist(:), 10, ref_anom(:)*1e-6, 's', 'filled')
datetick('x');
set(ha(11), 'box', 'on', 'TickLength', [0.01 0.01], 'YTick',  distticks,  ...
    'YTickLabel', num2str(distticks, '%3.0f'), 'TickDir', 'out', ...
    'XTick', datenum(1996:2020, 1 ,1)', ...
    'XTickLabel', datestr(datenum(1996:2020, 1 ,1)', 'yyyy'), ...
    'TickLength', [0.01 0.01], 'YDir', 'reverse');
ylabel('Distance from North [km]');
text(datenum(2005, 6, 1), 100, 'Reference', 'fontweight', 'bold');
axis(ha(11), [datenum(2005, 1, 2) datenum(2020, 4, 1) -20 870]);
ha(11).Position = pos7;
colormap(ha(11), cmap3);
caxis(ha(11), [min(colticks3) max(colticks3)]);
labels = string(ha(11).XAxis.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
ha(11).XAxis.TickLabels = labels; % set

axes(ha(14));
pos8 = ha(14).Position;
% r1 = plot(nanmean(ref.transp*1e-6), ref.psic, 'r', 'linewidth', 1); hold on
[hl4 hp4] = boundedline(nanmean(ref.transp*1e-6, 2), nanmean(ref.dist, 2), ...
    nanstd(ref.transp, 0 , 2)*1e-6./sqrt(sum(~isnan(ref.transp), 2)), '-b',...
    'alpha', 'Orientation', 'horiz');
hold on;
axis(ha(14), [-1 6.3 -20 870]); lim = axis;
p1 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], saf, ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2); 
p2 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], pf, ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
p3 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], ...
    saccf, ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
text(4, mean(saf), 'SAF', 'HorizontalAlignment', 'center');
text(4, mean(pf), 'PF', 'HorizontalAlignment', 'center');
text(4.5, mean(saccf), 'SACCF', 'HorizontalAlignment', 'center');
xlabel(ha(14), 'Mean transport [Sv]');
text(3, 200, 'Reference', 'fontweight', 'bold');
plot([0 0], [-20 870], '-k', 'linewidth', 0.5);
set(ha(14), 'YTick', ha(5).YTick, 'YTickLabel', [], 'TickDir', 'out', ...
    'XTick', [0:1:18]', 'XTickLabel', num2str([0:1:18]'), ...
    'TickLength', [0.03 0.03], 'box', 'on', 'YDir', 'reverse');
grid on;
leg3 = legend([hl4], [num2str(nansum(nanmean(ref.transp*1e-6, 2)),'%2.2f')...
    char(177) num2str(nanstd(nansum(ref.transp*1e-6, 1))/sqrt(size(ref.transp, 2)),'%2.2f') ...
    ' Sv'], 'fontsize', 5, 'fontname', 'SansSerif');
r1.Color = rgb('Blue'); 
r2.Color = r1.Color; 

axes(ha(15));
pos9 = ha(15).Position;
dref = nanmean(ref.dist, 2);
% plot(0, 10, 'r', 'linewidth', 1); hold on
[hl1 hp1] = boundedline(tr_ref*1e-6, dref, ci_ref*1e-6, '-b',...
    'alpha', 'Orientation', 'horiz');
hold on;
plot(tr_ref(sig_np_ref==1)*1e-6, dref(sig_np_ref==1), 'sb', ...
    'markerfacecolor', 'b', 'markersize', 5);
hl2.MarkerSize = 10;
% plot(nanmean(totxbt, 1), tot.psic, 'sr', 'markersize', 3);
axis(ha(15), [-0.35 0.35 -20 870]); 
lim = axis;
p1 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], saf, ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2); 
p2 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], pf, ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
p3 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], saccf, ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
text(0.20, mean(saf), 'SAF', 'HorizontalAlignment', 'center');
text(0.20, mean(pf), 'PF', 'HorizontalAlignment', 'center');
text(0.20, mean(saccf), 'SACCF', 'HorizontalAlignment', 'center');
text(-0.30, 200, 'Reference', 'fontweight', 'bold');
set(ha(15), 'YTick', ha(1).YTick, 'YTickLabel', [], 'TickDir', 'out', ...
    'XTick', [-0.4:0.1:0.4]', 'XTickLabel', num2str([-0.4:0.1:0.4]'), ...
    'TickLength', [0.03 0.03], ...
    'box', 'on', 'YDir', 'reverse');
plot([0 0], [-20 870], '-k', 'linewidth', 0.5);
xlabel(ha(15), 'Trend [Sv year^{-1}]')
grid on;

col = colorbar(ha(1)); 
col.Position(1) = sum(ha(5).Position([1, 3])) + 0.01;
col.Position(2) = ha(5).Position(2);
col.Position(4) = ha(5).Position(4);
col.Position(3) = 0.01;
cbarrow(col, cmap1);
ylabel(col, 'Transport anomaly [Sv]');
set(col, 'Ticks', colticks1', ...
    'TickLabels', num2str(colticks1'), 'TickLength', 0.03);
ha(5).Position = pos3;
labels = string(col.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
col.TickLabels = labels; % set

col = colorbar(ha(6)); 
col.Position(1) = sum(ha(10).Position([1, 3])) + 0.01;
col.Position(2) = ha(10).Position(2);
col.Position(4) = ha(10).Position(4);
col.Position(3) = 0.01;
cbarrow(col, cmap2);
ylabel(col, 'Transport anomaly [Sv]');
set(col, 'Ticks', colticks2', 'TickLabels', num2str(colticks2'), ...
    'TickLength', 0.03);
ha(10).Position = pos6;
labels = string(col.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
col.TickLabels = labels; % set

col = colorbar(ha(11)); 
col.Position(1) = sum(ha(15).Position([1, 3])) + 0.01;
col.Position(2) = ha(15).Position(2);
col.Position(4) = ha(15).Position(4);
col.Position(3) = 0.01;
cbarrow(col, cmap3);
set(col, 'Ticks', colticks3', ...
    'TickLabels', num2str(colticks3'), 'TickLength', 0.03);
ylabel(col, 'Transport anomaly [Sv]');
ha(15).Position = pos9;
labels = string(col.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
col.TickLabels = labels; % set

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 6);
set(findall(gcf, '-property', 'FontName'), 'FontName', 'SansSerif');
set(findall(gcf, '-property', 'TickDirection'), 'TickDirection', 'out');
leg.ItemTokenSize = [5 5];
leg.FontSize = 5;
leg.Position = [0.6376    0.7379    0.1029    0.0651];
leg.Box = 'off';
leg1a.ItemTokenSize = [5 5];
leg1a.FontSize = 5;
leg1a.Position = [0.8342    0.7728    0.0912    0.0618];
leg1a.Box = 'off';
leg1.ItemTokenSize = [5 5];
leg1.FontSize = 5;
leg1.Position = [0.6376    0.4638    0.1029    0.0651];
leg1.Box = 'off';
leg3.ItemTokenSize = [5 5];
leg3.FontSize = 5;
leg3.Position = [0.6376    0.1873    0.1029    0.0651];
leg3.Box = 'off';

annotation(fig0, 'textbox', [0.03 sum(ha(1).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'a', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig0, 'textbox', [ha(4).Position(1)-0.03 sum(ha(4).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'b', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig0, 'textbox', [ha(5).Position(1)-0.03 sum(ha(5).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'c', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig0, 'textbox', [0.03 sum(ha(6).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'd', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig0, 'textbox', [ha(9).Position(1)-0.03 sum(ha(9).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'e', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig0, 'textbox', [ha(10).Position(1)-0.03 sum(ha(10).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'f', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig0, 'textbox', [0.03 sum(ha(11).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'g', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig0, 'textbox', [ha(14).Position(1)-0.03 sum(ha(14).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'h', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig0, 'textbox', [ha(15).Position(1)-0.03 sum(ha(15).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'i', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
print('-dpdf', './FiguresPaper/figSI0_transpdist_mean.pdf', '-r500');
