% Code for plotting Figure 5 - Time series transport per pair of 
% Streamline. Time series of total, geostrophic and reference transport are
% estimated per pair of SSH streamlines. SSH intervals for the frontal
% regions are defined in Gutierrez-Villanueva et al. (2020). 
% Oct 13/2022
% March 3/2023 - Renames streamlines as SSH [m]
% Jun 22/2023 - Removes mean from time series and changes color of second
% line in panel b.
% Aug 28/2023 - Adds an extra column to show the trends for each transport
% component.

clear all; 
close all;

ndays = 19;

load(['./Datasets/Streamwise_binning_geos.mat']);
geos = stream;
load(['./Datasets/treamwise_binning_os38nb.mat']);
tot = stream;
load(['./Datasets/Streamwise_binning_ref.mat']);
ref = stream;
clear stream

% pycmd = '/Users/manuelgutierrez/anaconda3/bin/python';

alpha = 0.05;
ti = datenum(2005, 10, 1);
tf = datenum(2019, 4, 31);
colticks1 = [-24:2:24];
colticks2 = [-6:1:6];
colticks3 = [-24:3:24];
% cmap1 = getPyPlot_cMap('seismic', length(colticks1)-1, [], pycmd);
% cmap2 = getPyPlot_cMap('seismic', length(colticks2)-1, [], pycmd);
% cmap3 = getPyPlot_cMap('seismic', length(colticks3)-1, [], pycmd);
timeticks = datenum(2004:2:2020, 1, 1);
op_season = 1;

% Loads colormaps
load ./colormaps/Fig5_cmap.mat;

% Mean streamwise position of the ACC fronts
saf_psi = [-0.40 -0.60];%Subantarctic Front
pf_psi = [-1.00  -1.30];%Polar Front
saccf_psi = [-1.55 -1.65];%Southern ACC Front
% psiticks = [-1.7:0.2:-0.1]';
psi_saf = -0.52;
psi_pf = -1.25;
psi_saccf = -1.6;
psiticks = flip([saf_psi -0.80 pf_psi -1.45 -1.65])';


% % Calculate anomalies
[nPsi_g ntime_g] = meshgrid(geos.psic, geos.time);
% anom_g = [geos.transp-nanmean(geos.transp)]*1e-6; % geostrophic
% 
[nPsi_t ntime_t] = meshgrid(tot.psic, tot.time);
% anom_t = [tot.transp-nanmean(tot.transp)]*1e-6; % total
% 
[nPsi_r ntime_r] = meshgrid(ref.psic, ref.time);
% anom_r = [ref.transp-nanmean(ref.transp)]*1e-6; % reference

% Sum of reference and geostrophic
[a b] = meshgrid(geos.time, ref.time);
[ta1 ta2] = meshgrid(1:length(geos.time), 1:length(ref.time));
ind_gbci = ta1(abs(a-b)<=1);
ind_ref = ta2(abs(a-b)<=1);

totsum = [ref.transp(ind_ref, :) + geos.transp(ind_gbci, :)]*1e-6;

% Gets only those total transp transects with a coincident xbt transect
% [a b] = meshgrid(geos.time, tot.time);
% [ta1 ta2] = meshgrid(1:length(geos.time), 1:length(tot.time));
% ind_gbci = ta1(abs(a-b)<=1);
% ind_tot = ta2(abs(a-b)<=1);

% Loads coincident transects list
listcoin = load('./Datasets/coinc_xbt_adcp_trans.txt');
timecoin = datenum(listcoin);

rr = 0;
for i = 1:length(timecoin);
    iix = find(abs(timecoin(i) - tot.time)<2);

    if ~isempty(iix)
        rr = rr + 1;
        ind_tot(rr) = iix;
    end
end

totxbt = tot.transp(ind_tot, :)*1e-6; 

% Calculates anomalies
tot_anom = tot.transp - nanmean(tot.transp);
geos_anom = geos.transp - nanmean(geos.transp);
ref_anom = ref.transp - nanmean(ref.transp);

%%%%%% Calculate trends geostrophic %%%%%%
[tr_geo, ci_geo, ~, tr_np_geo sig_np_geo] = calctrends(geos.transp'*1e-6, ...
    [geos.time-datenum(1996, 1, 1)]/365.25, alpha, op_season);

%%%%%%% Calculate trends total %%%%%%%%
[tr_tot, ci_tot, ~, tr_np_tot sig_np_tot] = calctrends(tot.transp'*1e-6, ...
    [tot.time-datenum(1996, 1, 1)]/365.25, alpha, op_season);
[tr_tot2, ci_tot2, ~, tr_np_tot2 sig_np_tot2] = calctrends(tot.transp(ind_tot, :)'*1e-6, ...
    [tot.time(1, ind_tot)-datenum(1996, 1, 1)]/365.25, alpha, op_season);

%%%%%%%% Calculate trends reference %%%%%%%%
[tr_ref, ci_ref, ~, tr_np_ref sig_np_ref] = calctrends(ref.transp'*1e-6, ...
    [ref.time-datenum(1996, 1, 1)]/365.25, alpha, op_season);


% Figures
fig0 = figure('color', 'w');
figure_width = 18;
figure_height = 8.5;
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
ha(15).Position(1) = ha(15).Position(1)+0.01;

ha(2).Visible = 'off';
ha(3).Visible = 'off';
ha(7).Visible = 'off';
ha(8).Visible = 'off';
ha(12).Visible = 'off';
ha(13).Visible = 'off';



axes(ha(1));
pos1 = ha(1).Position;
p1 = patch([ones(1, 2)*datenum(1995, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    [saf_psi flip(saf_psi)], rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2); hold on;
p2 = patch([ones(1, 2)*datenum(1995, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    [pf_psi flip(pf_psi)], rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2);
p3 = patch([ones(1, 2)*datenum(1995, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    [saccf_psi flip(saccf_psi)], rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2);
scatter(ntime_t(:), nPsi_t(:), 40, tot_anom(:)*1e-6, 's', 'filled')
datetick('x');
set(ha(1), 'box', 'on', 'TickLength', [0.01 0.01], 'YTick',  psiticks,  ...
    'YTickLabel', num2str(psiticks, '%1.2f'), 'TickDir', 'out', ...
    'XTick', datenum(1996:2020, 1 ,1)', 'XTickLabel', [], ...
    'TickLength', [0.01 0.01]);
ylabel('SSH [m]');
text(datenum(2005, 6, 1), -0.27, 'Total', 'fontweight', 'bold');
axis(ha(1), [datenum(2005, 6, 1) datenum(2019, 09, 1) -1.75 -0.15]);
ha(1).Position = pos1;
colormap(ha(1), cmap1);
caxis(ha(1), [min(colticks1) max(colticks1)]);



axes(ha(4));
pos2 = ha(4).Position;
% plot(nanmean(tot.transp*1e-6), tot.psic, 'r', 'linewidth', 1); hold on
[hl1 hp1] = boundedline(nanmean(tot.transp*1e-6), tot.psic, ...
    nanstd(tot.transp, 0 ,1)*1e-6./sqrt(sum(~isnan(tot.transp), 1)), '-r',...
    'alpha', 'Orientation', 'horiz');
hold on;
[hl2 hp2] = boundedline(nanmean(totsum), tot.psic, ...
    nanstd(totsum, 0 ,1)./sqrt(sum(~isnan(totsum), 1)), '-',...
    'alpha', 'Orientation', 'horiz');
hl2.MarkerSize = 10;
hl2.Color = rgb('Green');
hp2.FaceColor = rgb('Green');
% plot(nanmean(totxbt, 1), tot.psic, 'sr', 'markersize', 3);
axis(ha(4), [-1 8.8 -1.75 -0.15]); lim = axis;
p1 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], [saf_psi flip(saf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2); 
p2 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], [pf_psi flip(pf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
p3 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], ...
    [saccf_psi flip(saccf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
text(8.05, mean(saf_psi), 'SAF', 'HorizontalAlignment', 'center');
text(8.2, mean(pf_psi), 'PF', 'HorizontalAlignment', 'center');
text(7.51, psi_saccf, 'SACCF', 'HorizontalAlignment', 'center');
text(6.2, -0.27, 'Total', 'fontweight', 'bold');
set(ha(4), 'YTick', ha(1).YTick, 'YTickLabel', [], 'TickDir', 'out', ...
    'XTick', [0:3:18]', 'XTickLabel', [], 'TickLength', [0.03 0.03], ...
    'box', 'on');
plot([0 0], [-1.75 -0.15], '-k', 'linewidth', 0.5);
grid on;
leg = legend([hl1; hl2], [num2str(nansum(nanmean(tot.transp*1e-6)),'%2.2f')...
    char(177) num2str(nanstd(nansum(tot.transp'*1e-6))/sqrt(size(tot.transp, 1)),'%2.2f') ...
    ' Sv'], [num2str(nanmean(nansum(totxbt')),'%2.2f')...
    char(177) num2str(nanstd(nansum(totxbt'))/sqrt(size(totsum, 1)),'%2.2f') ...
    ' Sv'], 'fontsize', 5, 'fontname', 'SansSerif');

axes(ha(5));
pos3 = ha(5).Position;
% plot(0, 10, 'r', 'linewidth', 1); hold on
[hl1 hp1] = boundedline(tr_tot, tot.psic, ci_tot*1e-6, '-r',...
    'alpha', 'Orientation', 'horiz');
hold on;
plot(tr_tot(sig_np_tot==1), tot.psic(sig_np_tot==1), 'sr', ...
    'markerfacecolor', 'r', 'markersize', 5);
hl2.MarkerSize = 10;
[hl2 hp2] = boundedline(tr_tot2, tot.psic, ci_tot2, '-b',...
    'alpha', 'Orientation', 'horiz');
pp2 = plot(tr_tot2(sig_np_tot2==1), tot.psic(sig_np_tot2==1), 's', ...
    'markerfacecolor', rgb('Green'), 'markeredgecolor', rgb('Green'), ...
    'markersize', 5);
hl2.Color = pp2.MarkerFaceColor;
hp2.FaceColor = pp2.MarkerFaceColor;
% plot(nanmean(totxbt, 1), tot.psic, 'sr', 'markersize', 3);
axis(ha(5), [-0.35 0.35 -1.75 -0.34]); 
lim = axis;
p1 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], [saf_psi flip(saf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2); 
p2 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], [pf_psi flip(pf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
p3 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], [saccf_psi flip(saccf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
text(0.20, mean(psi_saf), 'SAF', 'HorizontalAlignment', 'center');
text(0.20, mean(psi_pf), 'PF', 'HorizontalAlignment', 'center');
text(0.20, mean(psi_saccf), 'SACCF', 'HorizontalAlignment', 'center');
text(-0.30, -0.51, 'Total', 'fontweight', 'bold');
set(ha(5), 'YTick', ha(1).YTick, 'YTickLabel', [], 'TickDir', 'out', ...
    'XTick', [-0.3:0.15:0.3]', 'XTickLabel', [], 'TickLength', [0.03 0.03], ...
    'box', 'on');
plot([0 0], [-1.75 -0.15], '-k', 'linewidth', 0.5);
leg1a = legend([hl1; hl2], 'All transects', 'Subsampled', ...
    'fontsize', 5, 'fontname', 'SansSerif');
grid on;


axes(ha(6));
pos4 = ha(6).Position;
p1 = patch([ones(1, 2)*datenum(1995, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    [saf_psi flip(saf_psi)], rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2); hold on;
p2 = patch([ones(1, 2)*datenum(1995, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    [pf_psi flip(pf_psi)], rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2);
p3 = patch([ones(1, 2)*datenum(1995, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    [saccf_psi flip(saccf_psi)], rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2);
scatter(ntime_g(:), nPsi_g(:), 40, geos_anom(:)*1e-6, 's', 'filled')
datetick('x');
set(ha(6), 'box', 'on', 'TickLength', [0.01 0.01], 'YTick',  psiticks,  ...
    'YTickLabel', num2str(psiticks, '%1.2f'), 'TickDir', 'out', ...
    'XTick', datenum(1996:2020, 1 ,1)', 'XTickLabel', [], ...
    'TickLength', [0.01 0.01]);
ylabel('SSH [m]');
text(datenum(2005, 6, 1), -0.27, 'Geostrophic', 'fontweight', 'bold');
axis(ha(6), [datenum(2005, 6, 1) datenum(2019, 09, 1) -1.75 -0.15]);
ha(6).Position = pos4;
colormap(ha(6), cmap2);
caxis(ha(6), [min(colticks2) max(colticks2)]);

axes(ha(9));
pos5 = ha(9).Position;
% plot(nanmean(geos.transp*1e-6), geos.psic, 'k', 'linewidth', 1); hold on
[hl3 hp3] = boundedline(nanmean(geos.transp*1e-6), geos.psic, ...
    nanstd(geos.transp, 0 ,1)*1e-6./sqrt(sum(~isnan(geos.transp), 1)), '-k',...
    'alpha', 'Orientation', 'horiz');
hold on;
axis(ha(9), [-1 8.8 -1.75 -0.15]); lim = axis;
p1 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], [saf_psi flip(saf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2); 
p2 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], [pf_psi flip(pf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
p3 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], ...
    [saccf_psi flip(saccf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
text(8.05, mean(saf_psi), 'SAF', 'HorizontalAlignment', 'center');
text(8.2, mean(pf_psi), 'PF', 'HorizontalAlignment', 'center');
text(7.51, psi_saccf, 'SACCF', 'HorizontalAlignment', 'center');
text(4.0, -0.27, 'Geostrophic', 'fontweight', 'bold');
set(ha(9), 'YTick', ha(5).YTick, 'YTickLabel', [], 'TickDir', 'out', ...
    'XTick', [0:3:18]', 'XTickLabel', [], 'TickLength', [0.03 0.03], ...
    'box', 'on');
plot([0 0], [-1.75 -0.15], '-k', 'linewidth', 0.5);
grid on; 
leg2 = legend([hl3], [num2str(nansum(nanmean(geos.transp*1e-6)),'%2.2f')...
    char(177) num2str(nanstd(nansum(geos.transp'*1e-6))/sqrt(size(geos.transp, 1)),'%2.2f') ...
    ' Sv'], 'fontsize', 5, 'fontname', 'SansSerif');

axes(ha(10));
pos6 = ha(10).Position;
% dgeo = nanmean(geos.dist_u, 2);
% plot(0, 10, 'r', 'linewidth', 1); hold on
[hl1 hp1] = boundedline(tr_geo, geos.psic, ci_geo*1e-6, '-k',...
    'alpha', 'Orientation', 'horiz');
hold on;
plot(tr_geo(sig_np_geo==1), geos.psic(sig_np_geo==1), 'sk', ...
    'markerfacecolor', 'k', 'markersize', 5);
hl2.MarkerSize = 10;
% plot(nanmean(totxbt, 1), tot.psic, 'sr', 'markersize', 3);
axis(ha(10), [-0.35 0.35 -1.75 -0.34]); 
lim = axis;
p1 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], [saf_psi flip(saf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2); 
p2 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], [pf_psi flip(pf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
p3 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], ...
    [saccf_psi flip(saccf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
text(0.20, mean(psi_saf), 'SAF', 'HorizontalAlignment', 'center');
text(0.20, mean(psi_pf), 'PF', 'HorizontalAlignment', 'center');
text(0.20, mean(psi_saccf), 'SACCF', 'HorizontalAlignment', 'center');
text(-0.34, -0.51, 'Geostrophic', 'fontweight', 'bold');
set(ha(10), 'YTick', ha(1).YTick, 'YTickLabel', [], 'TickDir', 'out', ...
    'XTick', [-0.3:0.15:0.3]', 'XTickLabel', [], 'TickLength', [0.03 0.03], ...
    'box', 'on');
plot([0 0], [-1.75 -0.15], '-k', 'linewidth', 0.5);
grid on;


axes(ha(11));
pos7 = ha(11).Position;
p1 = patch([ones(1, 2)*datenum(1995, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    [saf_psi flip(saf_psi)], rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2); hold on;
p2 = patch([ones(1, 2)*datenum(1995, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    [pf_psi flip(pf_psi)], rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2);
p3 = patch([ones(1, 2)*datenum(1995, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    [saccf_psi flip(saccf_psi)], rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2);
scatter(ntime_r(:), nPsi_r(:), 40, ref_anom(:)*1e-6, 's', 'filled')
datetick('x');
set(ha(11), 'box', 'on', 'TickLength', [0.01 0.01], 'YTick',  psiticks,  ...
    'YTickLabel', num2str(psiticks, '%1.2f'), 'TickDir', 'out', ...
    'XTick', datenum(1996:2020, 1 ,1)', ...
    'XTickLabel', datestr(datenum(1996:2020, 1 ,1)', 'yyyy'), ...
    'TickLength', [0.01 0.01]);
ylabel('SSH [m]');
text(datenum(2005, 6, 1), -0.27, 'Reference', 'fontweight', 'bold');
axis(ha(11), [datenum(2005, 6, 1) datenum(2019, 09, 1) -1.75 -0.15]);
ha(11).Position = pos7;
colormap(ha(11), cmap3);
caxis(ha(11), [min(colticks3) max(colticks3)]);
labels = string(ha(11).XAxis.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
ha(11).XAxis.TickLabels = labels; % set

axes(ha(14));
pos8 = ha(14).Position;
% r1 = plot(nanmean(ref.transp*1e-6), ref.psic, 'r', 'linewidth', 1); hold on
[hl4 hp4] = boundedline(nanmean(ref.transp*1e-6), ref.psic, ...
    nanstd(ref.transp, 0 ,1)*1e-6./sqrt(sum(~isnan(ref.transp), 1)), '-b',...
    'alpha', 'Orientation', 'horiz');
hold on;
axis(ha(14), [-1 8.8 -1.75 -0.15]); lim = axis;
p1 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], [saf_psi flip(saf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2); 
p2 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], [pf_psi flip(pf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
p3 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], ...
    [saccf_psi flip(saccf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
text(8.05, mean(saf_psi), 'SAF', 'HorizontalAlignment', 'center');
text(8.2, mean(pf_psi), 'PF', 'HorizontalAlignment', 'center');
text(7.51, psi_saccf, 'SACCF', 'HorizontalAlignment', 'center');
xlabel(ha(14), 'Mean transport [Sv]');
text(4.5, -0.27, 'Reference', 'fontweight', 'bold');
plot([0 0], [-1.75 -0.15], '-k', 'linewidth', 0.5);
set(ha(14), 'YTick', ha(5).YTick, 'YTickLabel', [], 'TickDir', 'out', ...
    'XTick', [0:3:18]', 'XTickLabel', num2str([0:3:18]'), ...
    'TickLength', [0.03 0.03], 'box', 'on');
grid on;
leg3 = legend([hl4], [num2str(nansum(nanmean(ref.transp*1e-6)),'%2.2f')...
    char(177) num2str(nanstd(nansum(ref.transp'*1e-6))/sqrt(size(ref.transp, 1)),'%2.2f') ...
    ' Sv'], 'fontsize', 5, 'fontname', 'SansSerif');
txt1 = text(1.43, -1.27, ['Ekman = ' num2str(nansum(nanmean(ref.ek*1e-6)), '%2.2f') char(177) ...
     num2str(nanstd(nansum(ref.ek'*1e-6))/sqrt(size(ref.transp, 2)),'%2.2f') ' Sv'], ...
    'fontsize', 5, 'fontname', 'SansSerif');
r1.Color = rgb('Blue'); 
r2.Color = r1.Color; 

axes(ha(15));
pos9 = ha(15).Position;
% plot(0, 10, 'r', 'linewidth', 1); hold on
[hl1 hp1] = boundedline(tr_ref, ref.psic, ci_ref*1e-6, '-b',...
    'alpha', 'Orientation', 'horiz');
hold on;
plot(tr_ref(sig_np_ref==1), ref.psic(sig_np_ref==1), 'sb', ...
    'markerfacecolor', 'b', 'markersize', 5);
hl2.MarkerSize = 10;
% plot(nanmean(totxbt, 1), tot.psic, 'sr', 'markersize', 3);
axis(ha(15), [-0.35 0.35 -1.75 -0.34]); 
lim = axis;
p1 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], [saf_psi flip(saf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2); 
p2 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], [pf_psi flip(pf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
p3 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], ...
    [saccf_psi flip(saccf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
text(0.20, mean(psi_saf), 'SAF', 'HorizontalAlignment', 'center');
text(0.20, mean(psi_pf), 'PF', 'HorizontalAlignment', 'center');
text(0.20, mean(psi_saccf), 'SACCF', 'HorizontalAlignment', 'center');
text(-0.30, -0.51, 'Reference', 'fontweight', 'bold');
set(ha(15), 'YTick', ha(1).YTick, 'YTickLabel', [], 'TickDir', 'out', ...
    'XTick', [-0.3:0.15:0.3]', 'XTickLabel', num2str([-0.3:0.15:0.3]'), ...
    'TickLength', [0.03 0.03], 'box', 'on');
plot([0 0], [-1.75 -0.15], '-k', 'linewidth', 0.5);
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
ha(4).Position = pos2;
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
ha(9).Position = pos5;
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
leg.Position = [0.5952 0.7689 0.1029 0.0686];
leg.Box = 'off';
leg1a.ItemTokenSize = [5 5];
leg1a.FontSize = 5;
leg1a.Position = [0.7644 0.7594 0.0912 0.0686];
leg1a.Box = 'off';
leg2.ItemTokenSize = [5 5];
leg2.FontSize = 5;
leg2.Position = [0.6364 0.5316 0.1029 0.0651];
leg2.Box = 'off';
leg3.ItemTokenSize = [5 5];
leg3.FontSize = 5;
leg3.Position = [0.6481 0.1651 0.1029 0.0651];
leg3.Box = 'off';
txt1.FontSize = 5;

annotation(fig0, 'textbox', [0.05 sum(ha(1).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'a', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig0, 'textbox', [ha(4).Position(1)-0.03 sum(ha(4).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'b', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig0, 'textbox', [ha(5).Position(1)-0.02 sum(ha(5).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'c', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig0, 'textbox', [0.05 sum(ha(6).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'd', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig0, 'textbox', [ha(9).Position(1)-0.03 sum(ha(9).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'e', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig0, 'textbox', [ha(10).Position(1)-0.02 sum(ha(10).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'f', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig0, 'textbox', [0.05 sum(ha(11).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'g', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig0, 'textbox', [ha(14).Position(1)-0.03 sum(ha(14).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'h', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig0, 'textbox', [ha(15).Position(1)-0.02 sum(ha(15).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'i', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
print('-dpdf', ['./FiguresPaper/fig5_SSH_timeseries_trends.pdf'], '-r500');
