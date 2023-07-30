% Code for calculating the trends using long-term averages by distance.
% Calculates trends for total, geostrophic and reference transport 
% averaged per month and transport averaged per year. Trends are evaluated 
% using two different methods:
% least-squares fit (ordinary and weighted for monthly and annual averages)
% and the Thein-Seil estimator. Confidence limits of trends are calculated
% using the t-student distribution at 95% confidence (least-squares fit)
% and the modified Mann-Kendall test. 
%
%
% 2022/10/05 - Manuel O. Gutierrez-Villanueva
% 2023/06/28 - Uses the corrected .mat files. Reference transport is now
% corrected and has no gaps due to interpolation error.

clear all;
close all;

% Define variables
datein = datenum(1996, 01, 01, 0, 0 ,0);
alpha = 0.05;
gsize = 25;

dateticks = datenum(1996:2:2020, 1, 1);
season_opt = 1;

din1 = datenum(1999, 11, 1);
din2 = datenum(2005, 10, 1);
dfin1 = datenum(2019, 4, 30);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loads os38 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads transport
load('TotalTransport_os38nbmaxz_760_2.mat'); %os38nb
transptotal_760 = totaltransp;
% transptotal_760.transpdist = diff(transptotal_760.cumtransp);

% load('TransportDP_adcp_os38nbmaxz210_2.mat');
% transpos38_210 = transpDP;
% transpos38_210.transpdist = diff(transpos38_210.cumtransp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%% GEOSTROPHIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removes incomplete transects
load('geostrophic_noaver.mat'); 
transpgeosbci_760 = geosbaroc; 
transpgeosbci_760.transp(:, isnan(transpgeosbci_760.cumtransp(1, :))) = NaN;
transpgeosbci_760.dist(:, find(isnan(transpgeosbci_760.cumtransp(1, :))==1)-1) = NaN;
transpgeosbci_760.meandist = nanmean(transpgeosbci_760.dist, 2);

clear barocDP totaltransp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% REFERENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('reference_os38nb_760.mat');
transpref_760 = geosref;
transpref_760.time = transpref_760.time(1, ...
    transpref_760.time>datenum(2005, 1, 1));
transpref_760.transp = transpref_760.transp(1, ...
    transpref_760.time>datenum(2005, 1, 1));
clear transpref
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%% Annual averages %%%%%%%%%%%%%%%%%%%%%%%%%%%
antime = [0:23]';
xx1 = [transpgeosbci_760.time - datein]/365.25;
xx3 = [transpref_760.time - datein]/365.25;
xx5 = [transptotal_760.time - datein]/365.25; 

% Removes seasonal cycle
[~, yy1, ~] = fitseasoncycle(transpgeosbci_760.nettransp*1e-6, xx1);
[~, yy3, ~] = fitseasoncycle(transpref_760.nettransp, xx3);
[~, yy5, ~] = fitseasoncycle(transptotal_760.nettransp*1e-6, xx5);


bci_anaver = nan(size(antime));
bci_anastd = bci_anaver;
ref1_anaver = bci_anaver;
ref1_anastd = bci_anaver;
ref2_anaver = bci_anaver;
ref2_anastd = bci_anaver;
ref3_anaver = bci_anaver;
ref3_anastd = bci_anaver;
nb150_anaver = bci_anaver;
nb150_anastd = bci_anaver;
os38_anaver = bci_anaver;
os38_anastd = bci_anaver;
os382_anaver = bci_anaver;
os382_anastd = bci_anaver;

for j = 1:length(antime);
    if sum(floor(xx1) == antime(j))>3
        [bci_anaver(j, 1) bci_anastd(j, 1) np] = nanmstd(...
            yy1(floor(xx1) == antime(j)));
    end
    if sum(floor(xx3) == antime(j))>3
        [ref2_anaver(j, 1) ref2_anastd(j, 1) np] = nanmstd(...
            yy3(floor(xx3) == antime(j)));
    end
    if sum(floor(xx5) == antime(j))>3
        [os38_anaver(j, 1) os38_anastd(j, 1) np] = nanmstd(...
            yy5(floor(xx5) == antime(j)));
    end
end

% Finds coincident XBT/ADCP transects
% os38
[a b] = meshgrid(transpgeosbci_760.time, transptotal_760.time);
[t1 t2] = meshgrid(1:length(transpgeosbci_760.time), ...
    1:length(transptotal_760.time));
ind_tot = t2(abs(a-b)<=4);



%%%%%%%%%%%%%%%%%%%%%%%%%% Trends %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geostrophic
[tr_bci conf_bci sigtls_bci trnp_bci signp_bci] = ...
    calctrends(transpgeosbci_760.nettransp*1e-6, xx1, alpha, season_opt);

[tr_bci2 conf_bci2 sigtls_bci2 trnp_bci2 signp_bci2] = ...
    calctrends(transpgeosbci_760.nettransp(xx1>=[din2-datein]/365.25)*1e-6, ...
    xx1(xx1>=[din2-datein]/365.25), alpha, season_opt);

% Reference os38 (760 m)
[tr_ref2 conf_ref2 sigtls_ref2 trnp_ref2 signp_ref2] = ...
    calctrends(transpref_760.nettransp, xx3, alpha, season_opt);

% Total (760 m)
[tr_tot conf_tot sigtls_tot trnp_tot signp_tot] = ...
    calctrends(transptotal_760.nettransp*1e-6, xx5, alpha, season_opt);

% Linear least-squares
yy = transptotal_760.nettransp*1e-6;
xx = [transptotal_760.time - datein]/365.25;
Xfull = [ones(length(xx), 1) xx(:)];
bregfull1 = regress(yy(:), Xfull);
fitline1 = Xfull*bregfull1;

yy = transpgeosbci_760.nettransp*1e-6;
xx = [transpgeosbci_760.time - datein]/365.25;
Xfull = [ones(length(xx), 1) xx(:)];
bregfull2 = regress(yy(:), Xfull);
fitline2 = Xfull*bregfull2;

yy = transpref_760.nettransp;
xx = [transpref_760.time - datein]/365.25;
Xfull = [ones(length(xx), 1) xx(:)];
bregfull3 = regress(yy(:), Xfull);
fitline3 = Xfull*bregfull3;

yy = transpgeosbci_760.nettransp(transpgeosbci_760.time>=din2)*1e-6;
xx = [transpgeosbci_760.time(transpgeosbci_760.time>=din2) - datein]/365.25;
Xfull = [ones(length(xx), 1) xx(:)];
bregfull4 = regress(yy(:), Xfull);
fitline4 = Xfull*bregfull4;


fig1 = figure('color','w');
figure_width = 7;
figure_height = 8;
figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
set(fig1,'Visible', figuresVisible);
set(fig1, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height]);
set(fig1, 'PaperPositionMode', 'auto');
ha = tight_subplot(3, 1, [0.08 0.01], [0.15 0.05],[0.20 0.05]);


axes(ha(1));
p1 = plot(transptotal_760.time, ...
    transptotal_760.nettransp*1e-6, '-r', 'linewidth', 1);
hold on;
[p4 p5] = boundedline([antime(~isnan(os38_anastd))+0.5]*365.25 + datein, ...
    os38_anaver(~isnan(os38_anastd)), ...
    os38_anastd(~isnan(os38_anastd)), 'sr', ...
    'alpha');
p4.LineWidth = 1.5;
p4.MarkerFaceColor = 'none';
p4.MarkerEdgeColor = p4.MarkerFaceColor;
p5.FaceAlpha = 0.2;
p2 = plot(transptotal_760.time, fitline1, '-r', 'linewidth', 2.5);
p2.Color = rgb('DarkRed');
ylabel('Transport [Sv]', 'FontSize', 7, 'FontName', 'SansSerif');
text(datenum(1998, 1, 1), 112, ['Trend = ' ...
    num2str(tr_tot, '%3.2f') char(177) num2str(abs(conf_tot), '%1.2f')...
    ' Sv year^{-1} ({\itp}>0.05)'], 'color', 'r', 'FontSize', 7, 'FontName', 'SansSerif');
% text(datenum(2015, 12, 18), 114, '', 'color', 'r', ...
%     'FontSize', 7, 'FontName', 'SansSerif');
axis(ha(1), [datenum(1996, 6 ,1) datenum(2020, 5, 1) 45 119]);
text(datenum(1997, 1, 1), 53, 'Total {\it{U_{tot}}} (0-760 m)',...
    'FontSize', 7, 'FontName', 'SansSerif');
datetick('x', 'yyyy', 'keeplimits');
set(ha(1), 'TickDir', 'out', 'TickLength', [0.015 0.015], ...
    'XTick', dateticks', 'XTickLabels', [], ...
    'FontSize', 7, 'FontName', 'SansSerif');
grid on;

axes(ha(2));
p1 = plot(transpgeosbci_760.time, ...
    transpgeosbci_760.nettransp*1e-6, '-k', 'linewidth', 1);
hold on;
[p4 p5] = boundedline([antime(~isnan(bci_anastd))+0.5]*365.25 + datein, ...
    bci_anaver(~isnan(bci_anastd)), ...
    bci_anastd(~isnan(bci_anastd)), 'sk', ...
    'alpha');
p4.LineWidth = 1.5; p4.Color(4) = 0.3; p4.MarkerFaceColor = 'none';
p5.FaceAlpha = 0.2;
p4.MarkerEdgeColor = p4.MarkerFaceColor;
p2 = plot(transpgeosbci_760.time, fitline2, '-k', 'linewidth', 2.5);
p3 = plot(transpgeosbci_760.time(transpgeosbci_760.time>=din2),...
    fitline4, '-k', 'linewidth', 1.5);
p3.Color = [0.6 0.6 0.6];
axis(ha(2), [datenum(1996, 6 ,1) datenum(2020, 5, 1) 6 23]);
datetick('x', 'yyyy', 'keeplimits');
ylabel('Transport [Sv]', 'FontSize', 7, 'FontName', 'SansSerif');
text(datenum(1998, 1, 1), 21, ['Trend = ' ...
    num2str(tr_bci, '%3.2f') char(177) num2str(abs(conf_bci), '%1.2f')...
    ' Sv year^{-1} ({\itp}>0.05)'], 'color', 'k', 'FontSize', 7, 'FontName', 'SansSerif');
% text(datenum(2015, 9, 1), 21.5, '', 'color', 'k', ...
%     'FontSize', 7, 'FontName', 'SansSerif');
text(datenum(1997, 1, 1), 8.2, 'Geostrophic {\it{U_{geo}}} (0-760 m)',...
    'FontSize', 7, 'FontName', 'SansSerif');
set(ha(2), 'TickDir', 'out', 'TickLength', [0.015 0.015], ...
    'XTick', dateticks', 'XTickLabel', [], 'FontSize', 7, ...
    'FontName', 'SansSerif');
grid on;

axes(ha(3));
a1 = plot(transpref_760.time, ...
    transpref_760.nettransp*1e-6, '-k', 'linewidth', 1);hold on;
[a4 a5] = boundedline([antime(~isnan(ref2_anastd))+0.5]*365.25 + datein, ...
    ref2_anaver(~isnan(ref2_anastd))*1e-6, ...
    ref2_anastd(~isnan(ref2_anastd))*1e-6, 'sk', ...
    'alpha');
a1.Color = rgb('Blue');
a4.LineWidth = 1.5; a4.MarkerFaceColor = 'none'; 
a4.MarkerEdgeColor = a4.MarkerFaceColor;
a5.FaceAlpha = 0.2; 
a5.FaceColor = rgb('Blue');
a2 = plot(transpref_760.time, fitline3*1e-6, '-k', 'linewidth', 2.5);
a2.Color = rgb('Blue');
ax = gca;
text(datenum(1998, 1, 1), 125, ['Trend = ' ...
    num2str(tr_ref2*1e-6, '%3.2f') char(177) num2str(abs(conf_ref2)*1e-6, '%1.2f')...
    ' Sv year^{-1} ({\itp}>0.05)'], 'color', rgb('Blue'), ...
    'FontSize', 7, 'FontName', 'SansSerif');
% text(datenum(2015, 11, 1), 99, '', 'color', rgb('Blue'), ...
%     'FontSize', 7, 'FontName', 'SansSerif');
axis(ha(3), [datenum(1996, 6 ,1) datenum(2020, 5, 1) -15 138]);
datetick('x', 'yyyy', 'keeplimits');
set(gca, 'TickDir', 'out', 'TickLength', [0.015 0.015], ...
    'XTick', dateticks', 'XTickLabel', datestr(dateticks', 'yyyy'),...
    'FontSize', 7, 'FontName', 'SansSerif');
ylabel(gca, 'Transport [Sv]', 'FontSize', 7, 'FontName', 'SansSerif');
text(datenum(1997, 1, 1), 1.5, 'Reference {\it{U_{ref}}} (0-760 m)',...
    'FontSize', 7, 'FontName', 'SansSerif');
set(ha(3), 'TickDir', 'out', 'TickLength', [0.015 0.015], ...
    'XTick', dateticks', 'XTickLabel', datestr(dateticks', 'yyyy'),...
    'FontSize', 7, 'FontName', 'SansSerif');
grid on;


% set(findall(fig1, '-property', 'FontSize'), 'FontSize', 7)
% set(findall(fig1, '-property', 'FontName'), 'FontName', 'SainsSerif')
annotation(fig1, 'textbox', [0 sum(ha(1).Position([2,4]))+0.03 ...
    0.0480 0.0353], 'String', 'a', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig1, 'textbox', [0 sum(ha(2).Position([2,4]))+0.03 ...
    0.0480 0.0353], 'String', 'b', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig1, 'textbox', [0 sum(ha(3).Position([2,4]))+0.03 ...
    0.0480 0.0353], 'String', 'c', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
print('-dpdf', './FiguresPaper/fig2_DPtransp_trend.pdf', '-r500');
