% Code for plotting Figure 6 - Trends in transport per pair of sea surface
% height streamlines. Trends estimated for the total, geostrophic and
% reference transport per pair of SSH streamlines. Trends estimated in the
% least-squares sense after removing the seasonal cycle. Significance tests
% are performed using the modified Mann-Kendall test. 
%
% 2022/10/14 - Created by Manuel O. Gutierrez-Villanueva
% 2022/03/03 - Changes xlabel to SSH.
% 2023/06/30 - Fixed an issue with the binning, transport calculation and
% dates picked up for the analysis. 


clear all;
close all;

pycmd = '/Users/manuelgutierrez/anaconda3/bin/python';


ndays = 19;%choose between 5, 20, 50, 100, 180 days
load(['Streamwise_binning_geos_' num2str(ndays) 'days.mat']);
geos = stream;
load(['Streamwise_binning_os38nb_' num2str(ndays) 'days.mat']);
tot = stream;
load(['Streamwise_binning_ref_' num2str(ndays) 'days.mat']);
ref = stream;

op_season = 1;
alpha = 0.05;
ti = datenum(2005, 10, 1);
tf = datenum(2019, 4, 31);
colticks1 = [-20:2:20];
colticks2 = [-5:5];
colticks3 = [-16:2:16];
cmap1 = getPyPlot_cMap('seismic', length(colticks1)-1, [], pycmd);
cmap2 = getPyPlot_cMap('seismic', length(colticks2)-1, [], pycmd);
cmap3 = getPyPlot_cMap('seismic', length(colticks3)-1, [], pycmd);

% Mean streamwise position of the ACC fronts
saf_psi = [-0.40 -0.60];%Subantarctic Front
pf_psi = [-1.00  -1.30];%Polar Front
saccf_psi = [-1.55 -1.65];%Southern ACC Front
% psiticks = [-1.7:0.2:-0.1]';
psi_saf = -0.52;
psi_pf = -1.25;
psi_saccf = -1.6;
psiticks = flip([saf_psi -0.80 pf_psi -1.40 saccf_psi])';

%%%%%% Calculate trends geostrophic %%%%%%
[tr_Ug, ci_Ug, ~, tr_np_Ug sig_np_Ug] = calctrends(geos.transp'*1e-6, ...
    [geos.time-datenum(1996, 1, 1)]/365.25, alpha, op_season);

%%%%%%% Calculate trends total %%%%%%%%
[tr_Ut, ci_Ut, ~, tr_np_Ut sig_np_Ut] = calctrends(tot.transp'*1e-6, ...
    [tot.time-datenum(1996, 1, 1)]/365.25, alpha, op_season);
[tr_Utmr, ci_Utmr, ~, tr_np_Utmr sig_np_Utmr] = calctrends(tot.transp(tot.indx_rep, :)'*1e-6, ...
    [tot.time(1, tot.indx_rep)-datenum(1996, 1, 1)]/365.25, alpha, op_season);

%%%%%%%% Calculate trends reference %%%%%%%%
[tr_Ur, ci_Ur, ~, tr_np_Ur sig_np_Ur] = calctrends(ref.transp'*1e-6, ...
    [ref.time-datenum(1996, 1, 1)]/365.25, alpha, op_season);

% Figures
fig1 = figure('color', 'w');
figure_width = 8;
figure_height = 5;
figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
set(fig1,'Visible', figuresVisible)
set(fig1, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
set(fig1, 'PaperPositionMode', 'auto');

ha = tight_subplot(1, 1, [0.04 0.04], [0.25 0.05], [0.20 0.03]);


axes(ha(1));
posa = ha(1).Position;
ylim(ha(1), [-0.31 0.31]); xlim(ha(1), [-1.75 -0.34]);
ax = axis;
p1 = patch([saf_psi flip(saf_psi)], [ones(1, 2)*ax(3) ones(1, 2)*ax(4)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2); hold on;
p2 = patch([pf_psi flip(pf_psi)], [ones(1, 2)*ax(3) ones(1, 2)*ax(4)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
p3 = patch([saccf_psi flip(saccf_psi)], [ones(1, 2)*ax(3) ones(1, 2)*ax(4)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
[l1 p1] = boundedline(tot.psic(2:end), tr_Ut(2:end), ci_Ut(2:end), '-r', 'alpha', ...
    'orientation', 'vert', 'nan', 'remove');
plot(tot.psic(sig_np_Ut==1), tr_Ut(sig_np_Ut==1), 'sr', ...
    'markerfacecolor', 'r', 'markersize', 5);
[l2 p2] = boundedline(geos.psic, tr_Ug, ci_Ug, '-k', 'alpha', ...
    'orientation', 'vert', 'nan', 'remove');
plot(geos.psic(sig_np_Ug==1), tr_Ug(sig_np_Ug==1), 'sk', ...
    'markerfacecolor', 'k', 'markersize', 5);
[l3 p3] = boundedline(ref.psic(2:end), tr_Ur(2:end), ci_Ur(2:end), '-b', 'alpha', ...
    'orientation', 'vert', 'nan', 'remove');
plot(ref.psic(sig_np_Ur==1), tr_Ur(sig_np_Ur==1), 'sb', ...
    'markerfacecolor', rgb('Blue'), 'markersize', 5);
plot([-1.75 -0.15], [0 0], '-k', 'linewidth', 0.5);
set(ha(1), 'Box', 'on', 'TickDir', 'out', 'TickLength', [0.02 0.02], ...
    'YTick', [-0.6:0.1:0.6]', 'YTickLabel', num2str([-0.6:0.1:0.6]'), ...
    'XDir', 'reverse', 'XTick', psiticks, 'XTickLabel', num2str(psiticks),...
    'fontsize', 7, 'fontname', 'SansSerif', 'XTickLabelRotation', -45);
% text(-0.35, 0.40, 'October 2005 - April 2019', 'fontweight', 'bold',...
%     'fontsize', 7, 'fontname', 'SanSerif');
text(mean(saf_psi), 0.25, 'SAF', 'HorizontalAlignment', 'center',...
    'fontsize', 7, 'fontname', 'SanSerif');
text(mean(pf_psi), 0.25, 'PF', 'HorizontalAlignment', 'center',...
    'fontsize', 7, 'fontname', 'SanSerif');
text(psi_saccf, 0.25, 'SACCF', 'HorizontalAlignment', 'center',...
    'fontsize', 7, 'fontname', 'SanSerif');
ylabel('Trend [Sv year^{-1}]');
xlabel('SSH [m]');
grid on; 
leg1 = legend([l1; l2; l3], 'Total', 'Geostrophic', 'Reference');
set(leg1, 'location', 'southeast', 'fontsize', 7, 'box', 'off', ...
    'NumColumns', 2);
leg1.Position = [0.2184 0.3239 0.5516 0.0353];
leg1.ItemTokenSize = [8, 8, 8];
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 6);
set(findall(gcf, '-property', 'FontName'), 'FontName', 'SansSerif');
set(findall(gcf, '-property', 'TickDirection'), 'TickDirection', 'out');

print('-dpdf', ['./FiguresPaper/fig5_trends_streamlines_' num2str(ndays) ...
    'days.pdf'], '-r500');