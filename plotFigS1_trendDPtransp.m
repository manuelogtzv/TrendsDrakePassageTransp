% Code for calculating and plotting trends in net Drake Passage transport
% using different number of transects.
%
%
% 2021/09/17 - Manuel O. Gutierrez-Villanueva
%
% 2022/07/13 - Fair comparisons in time for total transport using all
% different time series.
% 2022/07/18 - Option to remove seasonal cycle before calculating trends.
% 2023/08/28 - Dropped some coincident XBT/ADCP transects due to quality
% flags.

clear all;
close all;

% pycmd = '/Users/manuelgutierrez/anaconda3/bin/python';

% Define variables
datein = datenum(1996, 01, 01, 0, 0 ,0);
alpha = 0.05;
gsize = 25;
dateticks = datenum(1996:2:2020, 1, 1);
season_opt = 1;
op_mostreptran = 1;

din2 = datenum(2005, 10, 1);
dfin1 = datenum(2019, 4, 30);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loads os38 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads transport
load('./Datasets/TransportDP_adcp_os38nbmaxz780_2.mat'); %os38nb
tot = transpDP;
tot.tottransp = tot.tottransp';
tot.transpdist = diff(tot.cumtransp);

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
    indr = yj(zj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%% GEOSTROPHIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removes incomplete transects
load('./Datasets/TransportDP_baroc_noaverg.mat'); %baroclinic
geo = barocDP; 
geo.transpdist(:, isnan(geo.cumtransp(1, :))) = NaN;
geo.dist(:, find(isnan(geo.cumtransp(1, :))==1)-1) = NaN;
geo.meandist = nanmean(geo.dist, 2);

clear barocDP transpDP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% REFERENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('./Datasets/ReferenceTransp_DP.mat');
ref = transpref;
ref.time_os38 = ref.time_os38(1, ...
    ref.time_os38>datenum(2005, 1, 1));
ref.transp_os38 = ref.transp_os38(1, ...
    ref.time_os38>datenum(2005, 1, 1));
ref.transp_os38_2 = ref.transp_os38_2(1, ...
    ref.time_os38>datenum(2005, 1, 1));
clear transpref
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Renames variables for simplicity
% os38
tot760 = tot.tottransp(tot.time>=din2 & ...
    tot.time<=dfin1);
timetot760 = tot.time(tot.time>=din2 & ...
    tot.time<=dfin1);

% geostrophic
geo760 = geo.transp(1, geo.time>=din2 & ...
    geo.time<=dfin1);
timegeo760 = geo.time(1, geo.time>=din2 & ...
    geo.time<=dfin1);

% reference
ref760 = ref.transp_os38(1, ref.time_os38>=din2 & ...
    ref.time_os38<=dfin1);
timeref760 = ref.time_os38(1, ref.time_os38>=din2 & ...
    ref.time_os38<=dfin1);


% Finds coincident xbt and adcp transects
% [x y] = meshgrid(timegeo760, timetot760);
% [t1 t2] = meshgrid(1:length(timegeo760), 1:length(timetot760));
% ind_tot = t2(abs(x-y)<=4);

% Loads coincident transects list
listcoin = load('./Datasets/coinc_xbt_adcp_trans.txt');
timecoin = datenum(listcoin);

rr = 0;
for i = 1:length(timecoin);
    iix = find(abs(timecoin(i) - timetot760)<2);

    if ~isempty(iix)
        rr = rr + 1;
        ind_tot(rr) = iix;

%         iig = find(abs(timecoin(i) - geos.time)<7);
%         ind_geo2(rr) = iig;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%% Trends %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geostrophic
[tr_geo conf_geo sigtls_geo trnp_geo signp_geo] = ...
    calctrends(geo.transp*1e-6, [geo.time-datein]/365.25, alpha, season_opt);
[tr_geo_1 conf_geo_1 sigtls_geo_1 trnp_geo_1 signp_geo_1] = ...
    calctrends(geo760*1e-6, [timegeo760-datein]/365.25, alpha, season_opt);


% Reference 
[tr_ref conf_ref sigtls_ref trnp_ref signp_ref] = ...
    calctrends(ref.transp_os38, [ref.time_os38-datein]/365.25, alpha, season_opt);


% Total
[tr_tot conf_tot sigtls_tot trnp_tot signp_tot] = ...
    calctrends(tot.tottransp*1e-6, [tot.time-datein]/365.25, alpha, season_opt);
[tr_tot_1 conf_tot_1 sigtls_tot_1 trnp_tot_1 signp_tot_1] = ...
    calctrends(tot760*1e-6, [timetot760-datein]/365.25, alpha, season_opt);
[tr_tot_2 conf_tot_2 sigtls_tot_2 trnp_tot_2 signp_tot_2] = ...
    calctrends(tot760(ind_tot(:))*1e-6, ...
    [timetot760(ind_tot(:))-datein]/365.25, alpha, season_opt);
[tr_tot_3 conf_tot_3 sigtls_tot_3 trnp_tot_3 signp_tot_3] = ...
    calctrends(tot.tottransp(indr)*1e-6, ...
    [tot.time(indr)-datein]/365.25, alpha, season_opt);




% Figure trends all-together
fig3 = figure('color','w');
figure_width = 8;
figure_height = 7;
figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
set(fig3,'Visible', figuresVisible);
set(fig3, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height]);
set(fig3, 'PaperPositionMode', 'auto');
ha = tight_subplot(1, 1, [0.18 0.16], [0.25 0.05],[0.20 0.03]);

axes(ha(1))
errorbar(1, tr_tot, conf_tot, 'r', 'linewidth', 1); 
hold on;
a1 = plot(1, tr_tot, '^r', 'markersize', 6, 'markerfacecolor', 'none', 'DisplayName', 'All transects');
errorbar(1, tr_tot_1, conf_tot_1, 'r', 'linewidth', 1); 
a2 = plot(1, tr_tot_1, 'sr', 'markersize', 6, 'markerfacecolor', 'none', 'DisplayName', 'Oct 2005 - Apr 2019');
errorbar(1, tr_tot_2, conf_tot_2, 'r', 'linewidth', 1); 
a3 = plot(1, tr_tot_2, 'or', 'markersize', 6, 'markerfacecolor', 'none', 'DisplayName', 'Subsampled');
errorbar(1, tr_tot_3, conf_tot_3, 'r', 'linewidth', 1); 
a4 = plot(1, tr_tot_3, 'dr', 'markersize', 6, 'markerfacecolor', 'none', 'DisplayName', 'Most repeated line');

errorbar(2, tr_geo, conf_geo, 'k', 'linewidth', 1); 
hold on;
b1 = plot(2, tr_geo, '^k', 'markersize', 6, 'markerfacecolor', 'none');
errorbar(2, tr_geo_1, conf_geo_1, 'k', 'linewidth', 1); 
b2 = plot(2, tr_geo_1, 'sk', 'markersize', 6, 'markerfacecolor', 'none');

errorbar(3, tr_ref, conf_ref*1e-1, 'b', 'linewidth', 1); 
hold on;
c1 = plot(3, tr_ref*1e-1, 'sb', 'markersize', 6, 'markerfacecolor', 'none');

axis(ha(1), [0.5 4.5 -0.38 0.38]);
ax = axis;
plot(ax(1:2), [0 0], '-k', 'linewidth', 1.5)
grid on;
set(ha(1), 'XTick', [1:3]', ...
    'XTickLabel', {'Total'; 'Geostrophic'; 'Reference \times 10^-^1'}, ...
    'TickDir', 'out', 'TickLength', [0.02 0.02], 'fontsize', 7, ...
    'fontname', 'SansSerif');
ylabel('Trend Drake Passage transport [Sv year^{-1}]', 'fontsize', 7, ...
    'fontname', 'SansSerif');
text(2, -0.15, '{\itp}>0.05', 'fontsize', 7, 'fontname', 'SansSerif')

ax2 = axes;
ax2.Position = [0.64 0.30 0.38 0.38];
ax2.FontName = 'SansSerif';
m_proj('azimuthal equal-area', 'radius', 55, 'lat', -70, ...
    'long', -65, 'rot', 0);
m_proj('lambert', 'lon', [-68 -56], 'lat', [-65.5 -54]); 
hold on;
m_gshhs_i('patch', [.4 .4 .4], 'edgecolor', 'none');
m_grid('linestyle', 'none', 'tickdir', 'out', 'linewidth', 0.5,...
    'xtick', [], 'fontsize', 7, 'ytick', []);
lmg = m_plot(tot.lon(:)-360, tot.lat(:), '-', 'linewidth', 0.5, ...
    'Color', rgb('Black'));
mrl = m_plot([-64.9422 -62.4422], [-55.0384 -62.4384], '-r', ...
    'linewidth', 3);
mrl.Color(4) = 0.85;

%%%% Tricks to make desired legend. Multiple symbols for a single entry %%%
hL(1) = a1;
hL(2) = b1;
hB(1) = a2;
hB(2) = b2;
hB(3) = c1;
hA(1) = a3;
hC(1) = a4;

% Add legend for the first/main plot handle
hLegend = legend([hL(1); hB(1); hA(1); hC(1)]);
hLegend.Position = [0.4707    0.7404    0.4824    0.1894];
hLegend.FontSize = 6;
hLegend.FontName = 'SansSerif';
drawnow(); 

% have to render the internal nodes before accessing them
% Extract legend nodes/primitives
hLegendEntry = hLegend.EntryContainer.NodeChildren(4); % first/bottom row of legend
iconSet = hLegendEntry.Icon.Transform.Children.Children; % array of first/bottom row's icons (marker+line)
% Create a new icon marker to add to the icon set
newLegendIcon = copy(iconSet(1)); % copy the object (or look into making a matlab.graphics.primitive.world.Marker)
newLegendIcon.get % list all properties
newLegendIcon.Parent = iconSet(1).Parent; % set the parent, adding to the legend's icon draw set
% Mess with the new icon's properties to show how you want
newLegendIcon.EdgeColorData = uint8([0;0;0;0]); % rgba uint8
newLegendIcon.VertexData(1) = 0.23; % [0-1] within the icon's boundaries (not the +0.02)
newLegendIcon.VertexData(2) = 0.50; % [0-1] within the icon's boundaries (not the +0.02)

% Extract legend nodes/primitives
hLegendEntry = hLegend.EntryContainer.NodeChildren(3); % first/bottom row of legend
iconSet = hLegendEntry.Icon.Transform.Children.Children; % array of first/bottom row's icons (marker+line)
% Create a new icon marker to add to the icon set
newLegendIcon = copy(iconSet(1)); % copy the object (or look into making a matlab.graphics.primitive.world.Marker)
newLegendIcon.get % list all properties
newLegendIcon.Parent = iconSet(1).Parent; % set the parent, adding to the legend's icon draw set
% Mess with the new icon's properties to show how you want
newLegendIcon.EdgeColorData = uint8([0;0;0;0]); % rgba uint8
newLegendIcon.VertexData(1) = 0.23; % [0-1] within the icon's boundaries (not the +0.02)
newLegendIcon.VertexData(2) = 0.50; % [0-1] within the icon's boundaries (not the +0.02)

% Extract legend nodes/primitives
hLegendEntry = hLegend.EntryContainer.NodeChildren(3); % first/bottom row of legend
iconSet = hLegendEntry.Icon.Transform.Children.Children; % array of first/bottom row's icons (marker+line)
% Create a new icon marker to add to the icon set
newLegendIcon = copy(iconSet(1)); % copy the object (or look into making a matlab.graphics.primitive.world.Marker)
newLegendIcon.get % list all properties
newLegendIcon.Parent = iconSet(1).Parent; % set the parent, adding to the legend's icon draw set
% Mess with the new icon's properties to show how you want
newLegendIcon.EdgeColorData = uint8([0;0;255;0]); % rgba uint8
newLegendIcon.VertexData(1) = 0.73; % [0-1] within the icon's boundaries (not the +0.02)
newLegendIcon.VertexData(2) = 0.50; % [0-1] within the icon's boundaries (not the +0.02)

hLegend.Color = 'none';
set(0, 'DefaultFigureRendererMode', 'manual');

print(gcf,'./FiguresPaper/figS1_trendtottransp.png','-dpng','-r1000');
