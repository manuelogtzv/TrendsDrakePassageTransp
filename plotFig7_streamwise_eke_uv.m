% Code to plot SI Fig - Calculates the eke and uv_prime in a synoptic 
% geostrophic streamwise coordinate system. Loads EKE and U'V' time series
% calculated by removing the time-mean velocities estimated as in
% Gutierrez-Villanueva et al. (2020) for each transects and bins and 
% averages them per pair of SSH streamlines. Trends are estimated in the
% least-squares sense after removing the seasonal cycle. The modified
% Mann-Kendall test is employed to test significance.
% 
% 2023/03/10 - Created by Manuel O. Gutierrez-Villanueva
% 2023/07/15 - Colormaps changed
% 2023/08/27 - Merged EKE and EMF figures into a single one



clear all; 
close all;

lat0 = -55;
lon0 = -65;
omega = 2*pi/86400;
gsize = 25000;% Grid-box size in meters
overlapping = 0.0;% Overlapping 0.5 = 50%
depths_psi = [100 200];% Averages mean geostrophic velocities and psi over these two depths
psi_saf = -0.52;
psi_pf = -1.25;
psi_saccf = -1.6;
op_season= 1;



emf_cont = [-0.2:0.025:0.2];
eke_cont = [0:0.01:0.2];
% eke_cmap = getPyPlot_cMap('Reds', length(eke_cont)-1, [], pycmd);
% emf_cmap = getPyPlot_cMap('bwr', length(emf_cont)-1, [], pycmd);

% Loads colormaps
load ./colormaps/Fig7_cmap.mat

gsize = 25;
ngsize = 11;
alpha = 0.05;
op_season = 1; % removes seasonal cycle
op_vid = 0;
op_mostreptran = 1;
dPsi = 0.10;

% Mean streamwise position of the ACC fronts
saf_psi = [-0.40 -0.60];%Subantarctic Front
pf_psi = [-1.00 -1.30];%Polar Front
saccf_psi = [-1.55 -1.65];%Southern ACC Front
psiticks = flip([saf_psi -0.80 pf_psi -1.45 -1.65])';


S_unique = -[0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90...
        1.0 1.15 1.25 1.35 1.45 1.50 1.55 1.60 1.65 1.70];

%%%% Loads datasets %%%%
load ./Datasets/Streamwise_binning_eke.mat

Zmax = 798;

% Loads coincident transects list
listcoin = load('./Datasets/coinc_xbt_adcp_trans.txt');
timecoin = datenum(listcoin);

rr = 0;
for i = 1:length(timecoin);
    iix = find(abs(timecoin(i) - stream.time)<2);

    if ~isempty(iix)
        rr = rr + 1;
        ind_tot2(rr) = iix;
    end
end


% Calculate trends EKE
[tr_eke, ci_eke, ~, tr_np_eke sig_np_eke] = calctrends(stream.eke'/Zmax, ...
    [stream.time-datenum(1996, 1, 1)]/365.25, alpha, op_season);
[tr_ekemr, ci_ekemr, ~, tr_np_ekemr sig_np_ekemr] = calctrends(stream.eke(ind_tot2, :)'/Zmax, ...
    [stream.time(1, ind_tot2)-datenum(1996, 1, 1)]/365.25, alpha, op_season);

% Calculate trends uv
[tr_uv, ci_uv, ~, tr_np_uv sig_np_uv] = calctrends(stream.uv_prime'/Zmax, ...
    [stream.time-datenum(1996, 1, 1)]/365.25, alpha, op_season);
[tr_uvr, ci_uvr, ~, tr_np_uvr sig_np_uvr] = calctrends(stream.uv_prime(ind_tot2, :)'/Zmax, ...
    [stream.time(1, ind_tot2)-datenum(1996, 1, 1)]/365.25, alpha, op_season);

ppsi = repmat(stream.psic, size(stream.eke, 1), 1);
ttime = repmat(stream.time', 1, size(stream.eke, 2));


fig0 = figure('color', 'w');
figure_width = 18.5;
figure_height = 8;
figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
set(fig0,'Visible', figuresVisible)
set(fig0, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
set(fig0, 'PaperPositionMode', 'auto');


ax1 = subplot(2, 6, [1, 2, 3, 4]);
pos1 = ax1.Position;
p1 = patch([ones(1, 2)*datenum(2004, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    [saf_psi flip(saf_psi)], rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2); hold on;
p2 = patch([ones(1, 2)*datenum(2004, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    [pf_psi flip(pf_psi)], rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2);
p3 = patch([ones(1, 2)*datenum(2004, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    [saccf_psi flip(saccf_psi)], rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2);
scatter(ttime(:), ppsi(:), 55, (stream.eke(:)/Zmax), 's', 'filled')
ax1.Box = 'on';
set(ax1, 'FontSize', 7, 'FontName', 'SansSerif', 'TickDir', 'out', ...
    'YTick',  psiticks,  'YTickLabel', num2str(psiticks, '%1.2f'),...
    'TickLen', [0.01 0.01]);
datetick('x');
ylabel('SSH [m]', 'FontSize', 7, 'FontName', 'SansSerif')
xlabel('');
axis(ax1, [datenum(2005, 6, 1) datenum(2019, 09, 1) min(stream.psic)-0.09 ...
    max(stream.psic)+0.09]);
ax1.Position = pos1;


ax2 = subplot(2, 6, 5);
ax2.Position(1) = sum(ax1.Position([1, 3])) + 0.02;
ax2.Position(3) = ax2.Position(3)+ 0.02;
r1 = plot((nanmean(stream.eke/Zmax)), stream.psic, 'r', 'linewidth', 1); hold on
r2 = plot((nanmean(stream.eke(ind_tot2, :)/Zmax)), stream.psic, 'linewidth', 1); 
r2.Color = rgb('Green');
% plot(nanmean(transp_lmg*1e-6), stream.psic, '--b', 'linewidth', 1);
axis(ax2, [0 0.080 min(stream.psic)-0.09 max(stream.psic)+0.09]); lim = axis;
p1 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], [saf_psi flip(saf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2); hold on;
p2 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], [pf_psi flip(pf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
p3 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], ...
    [saccf_psi flip(saccf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
text(lim(1)+0.01, mean(saf_psi), 'SAF', 'FontSize', 7, 'FontName', 'SansSerif', ...
    'HorizontalAlignment', 'center');
text(lim(1)+0.02, mean(pf_psi), 'PF', 'FontSize', 7, 'FontName', 'SansSerif', ...
    'HorizontalAlignment', 'center');
text(lim(1)+0.03, psi_saccf, 'SACCF', 'FontSize', 7, 'FontName', 'SansSerif', ...
    'HorizontalAlignment', 'center');
xlabel(ax2, ['<EKE> [m^{2} s^{-2}]'], 'FontSize', 7, ...
    'FontName', 'SansSerif');
annotation(fig0, 'line', [0.663366498851793 0.612528263557676],   ...
    [0.522328529452862 0.522328529452862]);
set(ax2, 'YTick',  psiticks,  'YTickLabel', [], 'TickDir', 'out',...
    'FontSize', 7, 'FontName', 'SansSerif', 'TickLen', [0.03 0.03], ...
    'XTick', [0:0.03:0.08]');
grid on; %grid minor;
% plot([0 0], [-1.75 -0.15], '-k', 'linewidth', 0.5);
leg = legend([r1; r2], 'All transects', 'Subsampled', 'fontsize', 5, ...
    'FontName', 'SansSerif', 'box', 'off');
leg.Position = [0.6554    0.6285    0.0991    0.0618];
leg.ItemTokenSize = [5, 5, 5];


ax3 = subplot(2, 6, 6);
xlim([-0.35 0.35]); ylim([min(stream.psic)-0.09 max(stream.psic)+0.09]);
ax = axis;
pos3 = ax3.Position;
p1 = patch([ones(1, 2)*ax(1) ones(1, 2)*ax(2)], [saf_psi flip(saf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2); hold on;
p2 = patch([ones(1, 2)*ax(1) ones(1, 2)*ax(2)], [pf_psi flip(pf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
p3 = patch([ones(1, 2)*ax(1) ones(1, 2)*ax(2)], [saccf_psi flip(saccf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
text(ax(1)+0.1, mean(saf_psi), 'SAF', 'FontSize', 7, 'FontName', 'SansSerif', ...
    'HorizontalAlignment', 'center');
text(ax(1)+0.1, mean(pf_psi), 'PF', 'FontSize', 7, 'FontName', 'SansSerif', ...
    'HorizontalAlignment', 'center');
text(ax(1)-ax(1)+0.2, psi_saccf, 'SACCF', 'FontSize', 7, 'FontName', 'SansSerif', ...
    'HorizontalAlignment', 'center');
[a4 a5] = boundedline(tr_eke*1e2, stream.psic, ci_eke, '-r', 'alpha', ...
    'orientation', 'horiz', 'nan', 'remove');
plot(tr_eke(sig_np_eke==1)*1e2, stream.psic(sig_np_eke==1), 'sr', 'markerfacecolor', 'r', ...
    'markersize', 5);
[b4 b5] = boundedline(tr_ekemr*1e2, stream.psic, ci_ekemr*1e2, '-r', 'alpha', ...
    'orientation', 'horiz', 'nan', 'remove');
p1 = plot(tr_ekemr(sig_np_ekemr==1)*1e2, stream.psic(sig_np_ekemr==1), 's', ...
    'markerfacecolor', rgb('Green'), 'markeredgecolor', rgb('Green'),...
    'markersize', 5);
b4.Color = rgb('Green');
b5.FaceColor = rgb('Green');
set(ax3, 'Box', 'on', 'TickDir', 'out', 'TickLength', [0.02 0.02], ...
    'YTick',  psiticks,  'YTickLabel', [], 'FontSize', 7, ...
    'FontName', 'SansSerif', 'XTick', [-0.2 0 0.2], 'TickLen', [0.03 0.03]);
xlabel(ax3, 'Trend \times 10^{-2} [m^{2} s^{-2} year^{-1}]', ...
    'FontSize', 7, 'FontName', 'SansSerif');
grid on; %grid minor;
plot([0 0], [-1.75 -0.15], '-k', 'linewidth', 0.5);
xlim(ax3, [-0.35 0.35]);
ax3.Position = pos3;

ax4 = subplot(2, 6, [7, 8, 9, 10]);
pos4 = ax4.Position;
p1 = patch([ones(1, 2)*datenum(2004, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    [saf_psi flip(saf_psi)], rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2); hold on;
p2 = patch([ones(1, 2)*datenum(2004, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    [pf_psi flip(pf_psi)], rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2);
p3 = patch([ones(1, 2)*datenum(2004, 6, 1) ones(1, 2)*datenum(2020, 4, 1)],...
    [saccf_psi flip(saccf_psi)], rgb('Goldenrod'), 'edgecolor', 'none', ...
    'facealpha', 0.2);
scatter(ttime(:), ppsi(:), 55, (stream.uv_prime(:)/Zmax), 's', 'filled')
ax4.Box = 'on';
set(ax4, 'FontSize', 7, 'FontName', 'SansSerif',...
    'YTick',  psiticks,  'YTickLabel', num2str(psiticks, '%1.2f'), ...
    'TickDir', 'out', 'TickLen', [0.01 0.01]);
datetick('x');
ylabel('SSH [m]', 'FontSize', 7, 'FontName', 'SansSerif')
xlabel('');
axis(ax4, [datenum(2005, 6, 1) datenum(2019, 09, 1) min(stream.psic)-0.09 ...
    max(stream.psic)+0.09]);
ax4.Position = pos4;


ax5 = subplot(2, 6, 11);
pos5 = ax5.Position
r1 = plot((nanmean(stream.uv_prime/Zmax)), stream.psic, 'r', 'linewidth', 1); hold on
r2 = plot((nanmean(stream.uv_prime(ind_tot2, :)/Zmax)), stream.psic, 'linewidth', 1);
r2.Color = rgb('Green');
% plot(nanmean(transp_lmg*1e-6), stream.psic, '--b', 'linewidth', 1);
axis(ax5, [-0.016 0.016 min(stream.psic)-0.09 max(stream.psic)+0.09]); lim = axis;
p1 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], [saf_psi flip(saf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2); hold on;
p2 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], [pf_psi flip(pf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
p3 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], ...
    [saccf_psi flip(saccf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
text(lim(1)+0.01, mean(saf_psi), 'SAF', 'FontSize', 7, 'FontName', 'SansSerif', ...
    'HorizontalAlignment', 'center');
text(lim(1)+0.01, mean(pf_psi), 'PF', 'FontSize', 7, 'FontName', 'SansSerif', ...
    'HorizontalAlignment', 'center');
text(lim(1)+0.01, psi_saccf, 'SACCF', 'FontSize', 7, 'FontName', 'SansSerif', ...
    'HorizontalAlignment', 'center');
xlabel(ax5, '<EMF> [m^{2} s^{-2}]', 'FontSize', 7, 'FontName', 'SansSerif');
annotation(fig0, 'line', [0.663366498851793 0.612528263557676],   ...
    [0.048494569967907 0.048494569967907]);
set(ax5, 'YTick',  psiticks,  'YTickLabel', [],...
    'TickDir', 'out', 'FontSize', 7, 'FontName', 'SansSerif', ...
    'XTick', [-0.01 0 0.01], 'TickLen', [0.03 0.03]);
grid on; %grid minor;
plot([0 0], [-1.75 -0.15], '-k', 'linewidth', 0.5);
% leg = legend([r1; r2], 'All transects', 'Subsampled', 'fontsize', 5, ...
%     'FontName', 'SensSerif', 'box', 'off');
% leg.Position = [0.6864 0.5749 0.0782 0.1549];
% leg.ItemTokenSize = [5, 5, 5];
ax5.Position = pos5;


ax6 = subplot(2, 6, 12);
xlim([-0.48 0.48]); ylim([min(stream.psic)-0.09 max(stream.psic)+0.09]);
ax = axis;
% ax6 = axis;
pos6 = ax6.Position;
p1 = patch([ones(1, 2)*ax(1) ones(1, 2)*ax(2)], [saf_psi flip(saf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2); hold on;
p2 = patch([ones(1, 2)*ax(1) ones(1, 2)*ax(2)], [pf_psi flip(pf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
p3 = patch([ones(1, 2)*ax(1) ones(1, 2)*ax(2)], [saccf_psi flip(saccf_psi)], ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
text(ax(1)+0.2, mean(saf_psi), 'SAF', 'FontSize', 7, 'FontName', 'SansSerif', ...
    'HorizontalAlignment', 'center');
text(ax(1)+0.2, mean(pf_psi), 'PF', 'FontSize', 7, 'FontName', 'SansSerif', ...
    'HorizontalAlignment', 'center');
text(ax(1)-ax(1)+0.25, psi_saccf, 'SACCF', 'FontSize', 7, 'FontName', 'SansSerif', ...
    'HorizontalAlignment', 'center');
[a4 a5] = boundedline(tr_uv*1e2, stream.psic, ci_uv, '-r', 'alpha', ...
    'orientation', 'horiz', 'nan', 'remove');
plot(tr_uv(sig_np_uv==1)*1e2, stream.psic(sig_np_uv==1), 'sr', 'markerfacecolor', 'r', ...
    'markersize', 5);
[b4 b5] = boundedline(tr_uvr*1e2, stream.psic, ci_uvr*1e2, '-r', 'alpha', ...
    'orientation', 'horiz', 'nan', 'remove');
pp0 = plot(tr_uvr(sig_np_uvr==1)*1e2, stream.psic(sig_np_uvr==1), 'sr', ...
    'markerfacecolor', rgb('Green'), 'markeredgecolor', rgb('Green'),...
    'markersize', 5);
b4.Color = pp0.MarkerEdgeColor;
b5.FaceColor = pp0.MarkerEdgeColor;
set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [0.02 0.02], ...
    'YTick',  psiticks,  'YTickLabel', [],...
    'FontSize', 7, 'FontName', 'SansSerif', 'TickLen', [0.03 0.03]);
xlabel(ax6, 'Trend \times 10^{-2} [m^{2} s^{-2} year^{-1}]', ...
    'FontSize', 7, 'FontName', 'SansSerif');
grid on; %grid minor;
plot([0 0], [-1.75 -0.15], '-k', 'linewidth', 0.5);
xlim(ax6, [-0.48 0.48]);
ax6.Position = pos6;



% ax1.Position(2) = ax1.Position(2)+0.11;
% ax1.Position(4) = ax1.Position(4)-0.07;
% ax1.Position(2) = ax1.Position(2); ax2.Position(4) = ax2.Position(4);
% ax1.Position(2) = ax1.Position(2); ax3.Position(4) = ax2.Position(4);




ax1.Position(1) = ax1.Position(1)-0.05;
ax2.Position(1) = sum(ax1.Position([1, 3]))+0.015;
ax2.Position(3) = ax2.Position(3)+0.020;
ax3.Position(1) = ax3.Position(1)-0.035;
ax3.Position(3) = ax3.Position(3)+0.020;
ax4.Position(1) = ax4.Position(1)-0.05;
ax5.Position(1) = ax2.Position(1);
ax5.Position(3) = ax2.Position(3);
ax6.Position(1) = ax3.Position(1);
ax6.Position(3) = ax3.Position(3);

% ax1.Position(4) = ax1.Position(4)-0.07;
% ax2.Position(2) = ax1.Position(2); ax2.Position(4) = ax1.Position(4);
% ax3.Position(2) = ax1.Position(2); ax3.Position(4) = ax1.Position(4);
col = colorbar(ax1); 
col.Position(3) = 0.01; col.TickDirection = 'out'; col.TickLength = 0.02;
col.Position(1) = sum(ax3.Position([1, 3])) + 0.01;
col.Position(2) = ax3.Position(2);
col.Position(4) = ax3.Position(4);
col.Position(3) = 0.01;
ylabel(col, '<EKE> [m^{2} s^{-2}]', ...
    'FontSize', 7, 'FontName', 'SansSerif');
col.FontSize = 7;
col.FontName = 'SansSerif';
caxis(ax1, [eke_cont(1) eke_cont(end)]); 
colormap(ax1, eke_cmap);
cbarrow(col, eke_cmap, 'up');


col2 = colorbar(ax4); 
caxis(ax4, [emf_cont(1) emf_cont(end)]); 
colormap(ax4, emf_cmap);
col2.Position(3) = 0.01; col.TickDirection = 'out'; col.TickLength = 0.02;
col2.Position(1) = sum(ax6.Position([1, 3])) + 0.01;
col2.Position(2) = ax6.Position(2);
col2.Position(4) = ax6.Position(4);
col2.Position(3) = 0.01;
ylabel(col2, '<EMF> [m^2 s^{-2}]', 'FontSize', 7, 'FontName', 'SansSerif');
col2.FontSize = 7;
col2.FontName = 'SansSerif';
cbarrow(col2, emf_cmap);

annotation(fig0, 'textbox', [ax1.Position(1)-0.05 sum(ax1.Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'a', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig0, 'textbox', [ax2.Position(1)-0.02 sum(ax2.Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'b', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig0, 'textbox', [ax3.Position(1)-0.025 sum(ax3.Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'c', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig0, 'textbox', [ax4.Position(1)-0.05 sum(ax4.Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'd', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig0, 'textbox', [ax5.Position(1)-0.02 sum(ax5.Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'e', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig0, 'textbox', [ax6.Position(1)-0.025 sum(ax6.Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'f', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');

% print('-dpdf', ['./FiguresPaper/fig_trendsEKE_streamlines_' num2str(nday) ...
%     'days.pdf'], '-r500');


% fig2 = figure('color', 'w');
% figure_width = 18;
% figure_height = 4;
% figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
% set(fig2,'Visible', figuresVisible)
% set(fig2, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
% set(fig2, 'PaperPositionMode', 'auto');

print('-dpdf', ['./FiguresPaper/fig7_trends_SSH_eke_emf.pdf'], '-r500');
