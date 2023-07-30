% Code to calculate trends in total transport per pair of streamlines 
% using only those transects that fall along the most repeated crossing.
% Trends are estimated using least-squares after removing the seasonal
% cycle. Significance tests are done using the modified Mann-Kendall test.
%
% 2023/07/11 - Created by Manuel O. Gutierrez-Villanueva


clear all;
close all;

pycmd = '/Users/manuelgutierrez/anaconda3/bin/python';


% Define variables
datein = datenum(1996, 01, 01, 0, 0 ,0);
alpha = 0.05;
gsize = 25;
op_season = 1;
op_mostreptran = 1;

alpha = 0.05;

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
tot.lon = tot.lon(:, totaltransp.time>=din2 & totaltransp.time<=dfin1);
tot.lat = tot.lat(:, totaltransp.time>=din2 & totaltransp.time<=dfin1);
tot.time = tot.time(totaltransp.time>=din2 & totaltransp.time<=dfin1);
tot.transp(tot.dist>840, :) = NaN;
tot.dist(tot.dist>840) = NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Option for using the most repeated transects
if op_mostreptran == 1;
    fprintf('\n\n total number of transects: %s\n',...
        num2str(length(tot.time)));

    box1 = [ -64.606712283408200 -56.324379209501110;
             -64.534099360692622 -56.568973205588243;
             -64.373586584163448 -56.541796094911895;
             -64.438556041330017 -56.351556320177458];

    box2 = [-63.009254066180596 -61.670918367346943;
            -62.908300616937744 -61.915816326530610;
            -62.504486819966345 -61.752551020408163;
            -62.625630959057766 -61.426020408163268];

    ii1 = inpolygon(tot.lon-360, tot.lat, box1(:, 1), ...
        box1(:, 2));

    [~, yj] = find(ii1 == 1);
    yj = unique(yj);

    mm1 = inpolygon(tot.lon(:, yj)-360, tot.lat(:, yj), ...
        box2(:, 1), box2(:, 2));

    [~, zj] = find(mm1 == 1);
    zj = unique(zj);

    % HOw many transects along the most repeated line
    fprintf('\n\n total number of transects (most repeated): %s\n',...
        num2str(length(tot.time(yj(zj)))));
end

% Calculates trends
[tr_tot, ci_tot, ~, tr_np_tot sig_np_tot] = calctrends(tot.transp, ...
    [tot.time-datenum(1996, 1, 1)]/365.25, alpha, op_season);
[tr_tot_mr, ci_tot_mr, ~, tr_np_tot_mr, sig_np_tot_mr] = calctrends(tot.transp(:, yj(zj)), ...
    [tot.time(yj(zj))-datenum(1996, 1, 1)]/365.25, alpha, op_season);

% Figure: Plots anomaly time series of transport in across-Passage distance
fig0 = figure('color', 'w');
figure_width = 6;
figure_height = 6;
figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
set(fig0,'Visible', figuresVisible)
set(fig0, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
set(fig0, 'PaperPositionMode', 'auto');


ha = tight_subplot(1, 1, [0.03 0.03], [0.15 0.03], [0.18 0.10]);


axes(ha(1));
% plot(0, 10, 'r', 'linewidth', 1); hold on
[hl1 hp1] = boundedline(tr_tot*1e-6, tot.dist, ci_tot*1e-6, '-r',...
    'alpha', 'Orientation', 'horiz');
hold on;
plot(tr_tot(sig_np_tot==1)*1e-6, tot.dist(sig_np_tot==1), 'sr', ...
    'markerfacecolor', 'r', 'markersize', 8);
[hl2 hp2] = boundedline(tr_tot_mr*1e-6, tot.dist, ci_tot_mr*1e-6, '-c',...
    'alpha', 'Orientation', 'horiz');
hold on;
p1 = plot(tr_tot_mr(sig_np_tot_mr==1)*1e-6, tot.dist(sig_np_tot_mr==1), 'oc', ...
    'markerfacecolor', 'c', 'markersize', 5);
hl2.Color = rgb('DarkCyan');
hp2.FaceColor = hl2.Color;
p1.MarkerFaceColor = hl2.Color;
p1.MarkerEdgeColor = hl2.Color;
% plot(nanmean(totxbt, 1), tot.psic, 'sr', 'markersize', 3);
axis(ha(1), [-0.21 0.21 -20 870]); 
lim = axis;
p1 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], saf, ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2); 
p2 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], pf, ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
p3 = patch([ones(1, 2)*lim(1) ones(1, 2)*lim(2)], saccf, ...
    rgb('Goldenrod'), 'edgecolor', 'none', 'facealpha', 0.2);
text(0.15, mean(saf), 'SAF', 'HorizontalAlignment', 'center');
text(0.15, mean(pf), 'PF', 'HorizontalAlignment', 'center');
text(-0.15, mean(saccf), 'SACCF', 'HorizontalAlignment', 'center');
text(-0.15, 200, 'Total', 'fontweight', 'bold');
set(ha(1), 'YTick', [0:200:1000]', 'YTickLabel', num2str([0:200:1000]'), 'TickDir', 'out', ...
    'XTick', [-0.2:0.1:0.2]', 'XTickLabel', num2str([-0.2:0.1:0.2]'), ...
    'TickLength', [0.03 0.03], 'box', 'on', 'YDir', 'reverse');
plot([0 0], [-20 870], '-k', 'linewidth', 0.5);
xlabel(ha(1), 'Trend [Sv year^{-1}]')
ylabel(ha(1), 'Distance from North [km]')
leg = legend([hl1; hl2], 'All transects', 'MR line');
leg.ItemTokenSize = [5 5];
leg.FontSize = 5;
leg.Position = [0.5899    0.2330    0.3088    0.1147];
leg.Color = 'None';
leg.Box = 'off';
grid on;
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 6);
set(findall(gcf, '-property', 'FontName'), 'FontName', 'SansSerif');
set(findall(gcf, '-property', 'TickDirection'), 'TickDirection', 'out');
print('-dpdf', ['./FiguresPaper/figR3_trendsTotalMRLine.pdf'], '-r500');
