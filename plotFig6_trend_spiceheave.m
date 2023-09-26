% Code for estimating trends in the temperature and salinity as a function
% of gamma (neutral density). Trends are estimated using the ordinary 
% least-squares fit by fitting y* = a1 + a2*t + a3*cos(omega_an*t) 
% + a4*sin(omega_an*t) + a5*cos(omega_sem*t) + a6*sin(omega_sem*t), 
% where omega_sem and omega_an are the semiannual and annual 
% frequencies and t is the time vector. a2 provides the trend. 
% Code uses also the Thein-Seil regression estimate assuming 
% that data does not conform to any knon distribution function,
% and its significance is assessed with the modified Mann-Kendall test.
% 
% Manuel O. Gutierrez Villanueva
% 
% 2022/06/15
% 2022/11/4 - Plots trends only
clear all;
close all;



% Contours and colormaps
cont_pttrend = [-0.16:0.02:0.16];
cont_saltrend = [-16:2:16];
cont_gamgrtrend = [-3:0.5:3]*1e-7;
% cmap_pttrend = flipud(getPyPlot_cMap('RdYlBu', length(cont_pttrend) - 1, [], pycmd));
% cmap_saltrend = flipud(getPyPlot_cMap('PiYG', length(cont_saltrend) - 1, [], pycmd));
% cmap_gamgrtrend = getPyPlot_cMap('Spectral_r', length(cont_gamgrtrend)-1, [], pycmd);
colugeo = rgb('Aqua');
cont_ugeo = [0.025 0.05 0.10]*1e2;

% Loads colormaps
load ./colormaps/Fig6_cmap.mat

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

% Loads baroclinic velocities
load(['./Datasets/geostrophic_noaver.mat']);


% Renames T, S, gamma and other variables
sal = geosbaroc.sal(:, 1:end-1, :);
temp = geosbaroc.temp(:, 1:end-1, :);
dist = geosbaroc.dist(1:end-1, :);
gamma_surfaces = [26.6:0.05:28.1];
gpot = geosbaroc.gpot(:, 1:end-1, :);
ubci = geosbaroc.u(:, 1:end-1, :);
dist_u = geosbaroc.dist_u(1:end-1, :);
z = geosbaroc.z;
lon = geosbaroc.lon(1:end-1, :);
lat = geosbaroc.lat(1:end-1, :);
% gamma_surfaces = gamma_surfaces(1:28);

% Pre-allocates memory
pt0 = nan(size(sal));
dens0 = pt0;
gamma_n = pt0;
sal_gamma = nan(length(gamma_surfaces), size(sal, 2), size(sal, 3));
temp_gamma = sal_gamma;
z_gamma = sal_gamma;

for i = 1:length(geosbaroc.time)
    % Estimates potential temperature
    pt0(:, :, i) = gsw_pt0_from_t(squeeze(sal(:, :, i)), ...
        squeeze(temp(:, :, i)), z(:));
end


% Calculates potential temperature, density and neutral density
for i = 1:length(geosbaroc.time)
    if gamma ~= 1; %calculates potential density
        dens0(:, :, i) = sw_pden(squeeze(sal(:, :, i)), ...
            squeeze(temp(:, :, i)), z(:), 0) - 1000;
    else % calculates gamma
        gamma_n(:, :, i) = eos80_legacy_gamma_n(squeeze(sal(:, :, i)), ...
            squeeze(temp(:, :, i)), repmat(z, 1, size(sal, 2)),...
            repmat(lon(:, i)', size(z, 1), 1),...
            repmat(lat(:, i)', size(z, 1), 1));
        
        % Caculates temperature, salinity and depth for each given
        % isopycnal
        for m = 1:size(sal, 2)
            if sum(~isnan(gamma_n(:, m, i)))~=0;
                ind = ~isnan(gamma_n(:, m, i));
                sal_gamma(: , m, i) = interp1(squeeze(gamma_n(ind, m, i)), ...
                    squeeze(sal(ind, m, i)), gamma_surfaces(:)); %sal
                temp_gamma(: , m, i) = interp1(squeeze(gamma_n(ind, m, i)), ...
                    squeeze(temp(ind, m, i)), gamma_surfaces(:)); %temp
                z_gamma(: , m, i) = interp1(squeeze(gamma_n(ind, m, i)), ...
                    z(:), gamma_surfaces(:)); %gamma
                
            end
        end
    end
end

% Estimates means
meanPt = nanmean(temp_gamma, 3);
meanPt(sum(~isnan(temp_gamma), 3)<50) = nan;

meanSal = nanmean(sal_gamma, 3);
meanSal(sum(~isnan(sal_gamma), 3)<50) = nan;

meanZ = nanmean(z_gamma, 3);
meanZ(sum(~isnan(z_gamma), 3)<50) = nan;

% Mask for data that has less than 50 transects 
mask(1, :, :) = meanPt./meanPt;
mask = repmat(mask, [size(temp_gamma, 3) 1 1]);
mask = permute(mask, [2 3 1]);

% Estimates stratication 
dgamma = mean(diff(gamma_surfaces)); % delta gamma
dPt = diff(meanPt)/dgamma;
dSal = diff(meanSal)/dgamma;
dZ = diff(meanZ)/dgamma;

dPtdgam = nan(length(gamma_surfaces), size(meanPt, 2));
dSaldgam = dPtdgam;
dZdgam = dPtdgam;

% Interpolates stratification to neutral surfaces
for i=1:size(dPt, 2);
    if sum(~isnan(dPt(:, i)))>0
        dPtdgam(:, i) = interp1(0.5*[gamma_surfaces(1:end-1) + ...
            gamma_surfaces(2:end)]', dPt(:, i), gamma_surfaces(:), ...
            'spline', 'extrap');
        dSaldgam(:, i) = interp1(0.5*[gamma_surfaces(1:end-1) + ...
            gamma_surfaces(2:end)]', dSal(:, i), gamma_surfaces(:), ...
            'spline', 'extrap');
        dZdgam(:, i) = interp1(0.5*[gamma_surfaces(1:end-1) + ...
            gamma_surfaces(2:end)]', dZ(:, i), gamma_surfaces(:), ...
            'spline', 'extrap');
    end
end

dPtdgam(sum(~isnan(temp_gamma), 3)<50) = NaN;
dSaldgam(sum(~isnan(sal_gamma), 3)<50) = NaN;
dZdgam(sum(~isnan(z_gamma), 3)<50) = NaN;

% Interpolates instantaneous gamma to mean depth of each gamma surface
% Calculates potential temperature, density and neutral density
gamman = nan(size(meanZ, 1), size(meanZ, 2), size(gamma_n, 3));

for i = 1:length(geosbaroc.time)
    for m = 1:size(sal, 2)
        if sum(~isnan(meanZ(:, m)))>1 & sum(~isnan(gamma_n(:, m, i)))>1;
            ind = ~isnan(meanZ(:, m));
            gamman(ind, m, i) = interp1(z, squeeze(gamma_n(:, m, i)), ...
                meanZ(ind, m), 'spline');
        end
    end
end

%%%%%%%%%%%%%%%% T R E N D S %%%%%%%%%%%%%%%%%%
din2 = datenum(2005, 10, 1);
dfin1 = datenum(2019, 4, 30);

indx = find(geosbaroc.time>= din2 & geosbaroc.time <=dfin1);

%%% Spice or along isopycnals
% potential temperature
[trendgam_pt confintgam_pt sigtrgam_ls_pt trendgam_thse_pt ...
    sigtrgam_mk_pt] = calctrends2D(temp_gamma(:, :, indx).*mask(:, :, indx), [geosbaroc.time(indx)-datein]/365.25,...
    alpha, season_opt);

% salinity
[trendgam_sal confintgam_sal sigtrgam_ls_sal trendgam_thse_sal ...
    sigtrgam_mk_sal] = calctrends2D(sal_gamma(:, :, indx).*mask(:, :, indx), [geosbaroc.time(indx)-datein]/365.25,...
    alpha, season_opt);

% depth
[trendgam_z confintgam_z sigtrgam_ls_z trendfull_thse_z ...
    sigtrend_mk_z] = calctrends2D(z_gamma(:, :, indx).*mask(:, :, indx), [geosbaroc.time(indx)-datein]/365.25,...
    alpha);

%%% Heave
% gamma
[trendz_gam confintz_gam sigtrz_ls_gam trendz_thse_gam ...
    sigtrz_mk_gam] = calctrends2D(gamman(:, :, indx), [geosbaroc.time(indx)-datein]/365.25,...
    alpha, season_opt);

%%% Total
% potential temperature
[trendtot_pt confinttot_pt sigtrtot_ls_pt trendtot_thse_pt ...
    sigtrtot_mk_pt] = calctrends2D(pt0(:, :, indx), [geosbaroc.time(indx)-datein]/365.25,...
    alpha, season_opt);

% salinity
[trendtot_sal confinttot_sal sigtrtot_ls_sal trendtot_thse_sal ...
    sigtrtot_mk_sal] = calctrends2D(sal(:, :, indx), [geosbaroc.time(indx)-datein]/365.25,...
    alpha, season_opt);

%%% Spice + Heave
% Potential temperature
[trendtot2_pt confinttot2_pt sigtrtot2_ls_pt trendtot2_thse_pt ...
    sigtrtot2_mk_pt] = calctrends2D(temp_gamma.*mask + ...
    dPtdgam+gamman.*repmat(dPtdgam, [1 1 size(temp_gamma, 3)]), ...
    [geosbaroc.time-datein]/365.25, alpha, season_opt);

% salinity
[trendtot2_sal confinttot2_sal sigtrtot2_ls_sal trendtot2_thse_sal ...
    sigtrtot2_mk_sal] = calctrends2D(sal_gamma.*mask + ...
    dSaldgam+gamman.*repmat(dSaldgam, [1 1 size(sal_gamma, 3)]), ...
    [geosbaroc.time-datein]/365.25, alpha, season_opt);


d1 = repmat(nanmax(dist'), [size(meanZ, 1) 1]);


% Plots of spice, heave and total
fig1 = figure('color', 'w');
figure_width = 18;
figure_height = 9;
figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
set(gcf,'Visible', figuresVisible)
set(gcf, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
set(gcf, 'PaperPositionMode', 'auto');

h4 = tight_subplot(2, 3, [0.06 0.05], [0.15 0.05], [0.10 0.12]);

axPos = [h4(1).Position; h4(2).Position; h4(3).Position; h4(4).Position; ...
    h4(5).Position; h4(6).Position];
% d1 = nanmax(geosbaroc.dist);


axes(h4(1))
contourf(d1, -meanZ, trendgam_pt, cont_pttrend, ...
    'edgecolor', 'none')
caxis(h4(1), [cont_pttrend(1) cont_pttrend(end)]);
hold on
colormap(h4(1), cmap_pttrend)
[k1 k2] = contour(d1, -meanZ, trendgam_pt, ...
    [0 0], '-', 'LineColor', [0 0 0], 'LineWidth', 1.5);
clabel(k1, k2, 'labelspacing', 400, 'fontsize', 7, 'fontname', 'SansSerif');
axis(h4(1), [0 890 -760 0]);
% xlabel('Distance from North [km]');
ylabel('Depth [m]');
set(h4(1), 'TickLength', [0.025 0.025], 'TickDir', 'out', ... 
    'XTickLabel', [], 'TickDir', 'out');

% Hatched regions
hat0 = trendgam_pt;
hat0(sigtrgam_mk_pt==0) = 10;
[~, gh0] = contourf(d1, -meanZ, hat0, [10 10]);
set(gh0, 'linestyle', 'none', 'Tag', 'HatchingRegion');
hp = findobj(h4(1), 'Tag', 'HatchingRegion');
hh = hatchfill2(hp, 'cross', 'HatchAngle', 45, 'HatchDensity', 60, ...
    'HatchColor', [0.5 0.5 0.5], 'HatchLineWidth', 0.5, 'Fill', 'off');
text(h4(1), 50, -700, 'Spice \theta', 'BackgroundColor', 'w');
[zz1 zz2] = contour(nanmax(dist_u, [], 2)', -z, ...
    nanmean(ubci(:, :, indx), 3)*1e2, cont_ugeo, ...
    'linecolor', colugeo, 'linestyle', '-');
clabel(zz1, zz2, 'labelspacing', 400, 'color', colugeo, 'fontsize', 6);
h4(1).Position = axPos(1, :);
p1 = plot([820 745], [0 0], '-', 'linewidth', 4, 'color', 'm');
p2 = plot([300 500], [0 0], '-', 'linewidth', 4, 'color', rgb('DarkGray'));
p3 = plot([35 135], [0 0], '-', 'linewidth', 4, 'color', rgb('Goldenrod'));



axes(h4(2))
contourf(d1, -meanZ, trendz_gam.*dPtdgam, ...
    cont_pttrend, 'edgecolor', 'none')
caxis(h4(2), [cont_pttrend(1) cont_pttrend(end)]);
hold on
colormap(h4(2), cmap_pttrend)
[k1 k2] = contour(d1, -meanZ, trendz_gam.*dPtdgam, ...
    [0 0], '-', 'LineColor', [0 0 0], 'LineWidth', 1.5);
clabel(k1, k2, 'labelspacing', 400, 'fontsize', 7, 'fontname', 'SansSerif');
axis(h4(2), [0 890 -760 0]);
% xlabel('Distance from North [km]');
% ylabel('Depth [m]');
set(h4(2), 'TickLength', [0.025 0.025], 'TickDir', 'out', ... 
    'XTickLabel', [], 'YTickLabel', []);

% Hatched regions
hat1 = trendz_gam;
hat1(sigtrz_mk_gam==0) = 10;
[~, gh1] = contourf(d1, -meanZ, hat1, [10 10]);
set(gh1, 'linestyle', 'none', 'Tag', 'HatchingRegion');
hp = findobj(h4(2), 'Tag', 'HatchingRegion');
hh = hatchfill2(hp, 'cross', 'HatchAngle', 45, 'HatchDensity', 60, ...
    'HatchColor', [0.5 0.5 0.5], 'HatchLineWidth', 0.5, 'Fill', 'off');
text(h4(2), 50, -700, 'Heave \theta', 'BackgroundColor', 'w');
[zz1 zz2] = contour(nanmax(dist_u, [], 2)', -z, ...
    nanmean(ubci(:, :, indx), 3)*1e2, cont_ugeo, ...
    'linecolor', colugeo, 'linestyle', '-');
clabel(zz1, zz2, 'labelspacing', 400, 'color', colugeo, 'fontsize', 6);
h4(2).Position = axPos(2, :);
p1 = plot([820 745], [0 0], '-', 'linewidth', 4, 'color', 'm');
p2 = plot([300 500], [0 0], '-', 'linewidth', 4, 'color', rgb('DarkGray'));
p3 = plot([35 135], [0 0], '-', 'linewidth', 4, 'color', rgb('Goldenrod'));


axes(h4(3))
contourf(nanmax(dist, [], 2)', -z, ...
    trendtot_pt, cont_pttrend, 'edgecolor', 'none');
caxis(h4(3), [cont_pttrend(1) cont_pttrend(end)]);
hold on
colormap(h4(3), cmap_pttrend)
[k1 k2] = contour(nanmax(dist'), -z, ...
    trendtot_pt, [0 0], '-', ...
    'LineColor', [0 0 0], 'LineWidth', 1.5);
clabel(k1, k2, 'labelspacing', 400, 'fontsize', 7, 'fontname', 'SansSerif');
axis(h4(3), [0 890 -760 0]);
set(h4(3), 'TickLength', [0.025 0.025], 'TickDir', 'out', ... 
    'XTickLabel', [], 'YTickLabel', []);
col3 = colorbar(h4(3), 'Position', [sum(axPos(3, [1,3]))+0.01 ...
    axPos(3, 2) 0.01 axPos(3, 4)], 'Ticks', cont_pttrend(1:2:end)',...
    'TickLabels', num2str(cont_pttrend(1:2:end)'), 'TickLength', 0.025, ...
    'TickDirection', 'out');
ylabel(col3, 'Trend \theta [^oC year^-^1]');
cbarrow(col3, cmap_pttrend);



% Hatched regions
hat2= trendtot_pt;
hat2(sigtrtot_mk_pt==0) = 10;
[~, gh2] = contourf(nanmax(dist, [], 2)',...
    -z, hat2, [10 10]);
set(gh2, 'linestyle', 'none', 'Tag', 'HatchingRegion');
hp = findobj(h4(3), 'Tag', 'HatchingRegion');
hh = hatchfill2(hp, 'cross', 'HatchAngle', 45, 'HatchDensity', 60, ...
    'HatchColor', [0.5 0.5 0.5], 'HatchLineWidth', 0.5, 'Fill', 'off');
text(h4(3), 50, -700, 'Depth \theta', 'BackgroundColor', 'w');
[zz1 zz2] = contour(nanmax(dist_u, [], 2)', -z, ...
    nanmean(ubci(:, :, indx), 3)*1e2, cont_ugeo, ...
    'linecolor', colugeo, 'linestyle', '-');
clabel(zz1, zz2, 'labelspacing', 400, 'color', colugeo, 'fontsize', 6);
h4(3).Position = axPos(3, :);
p1 = plot([820 745], [0 0], '-', 'linewidth', 4, 'color', 'm');
p2 = plot([300 500], [0 0], '-', 'linewidth', 4, 'color', rgb('DarkGray'));
p3 = plot([35 135], [0 0], '-', 'linewidth', 4, 'color', rgb('Goldenrod'));



axes(h4(4))
contourf(d1, -meanZ, trendgam_sal*1e3, cont_saltrend, ...
    'edgecolor', 'none')
caxis(h4(4), [cont_saltrend(1) cont_saltrend(end)]);
hold on
colormap(h4(4), cmap_saltrend)
[k1 k2] = contour(d1, -meanZ, trendgam_sal*1e3, ...
    [0 0], '-', 'LineColor', [0 0 0], 'LineWidth', 1.5);
clabel(k1, k2, 'labelspacing', 400, 'fontsize', 7, 'fontname', 'SansSerif');
axis(h4(4), [0 890 -760 0]);
xlabel('Distance from North [km]');
ylabel('Depth [m]');
set(h4(4), 'TickLength', [0.025 0.025], 'TickDir', 'out');


% Hatched regions
hat3 = trendgam_sal;
hat3(sigtrgam_mk_sal==0) = 10;
[~, gh3] = contourf(d1, -meanZ, hat3, [10 10]);
set(gh3, 'linestyle', 'none', 'Tag', 'HatchingRegion');
hp = findobj(h4(4), 'Tag', 'HatchingRegion');
hh = hatchfill2(hp, 'cross', 'HatchAngle', 45, 'HatchDensity', 60, ...
    'HatchColor', [0.5 0.5 0.5], 'HatchLineWidth', 0.5, 'Fill', 'off');
text(h4(4), 50, -700, 'Spice {\itS}', 'BackgroundColor', 'w');
[zz1 zz2] = contour(nanmax(dist_u, [], 2)', -z, ...
    nanmean(ubci(:, :, indx), 3)*1e2, cont_ugeo, ...
    'linecolor', colugeo, 'linestyle', '-');
clabel(zz1, zz2, 'labelspacing', 400, 'color', colugeo, 'fontsize', 6);
h4(4).Position = axPos(4, :);
p1 = plot([820 745], [0 0], '-', 'linewidth', 4, 'color', 'm');
p2 = plot([300 500], [0 0], '-', 'linewidth', 4, 'color', rgb('DarkGray'));
p3 = plot([35 135], [0 0], '-', 'linewidth', 4, 'color', rgb('Goldenrod'));


axes(h4(5))
contourf(d1, -meanZ, trendz_gam.*dSaldgam*1e3, ...
    cont_saltrend, 'edgecolor', 'none')
caxis(h4(5), [cont_saltrend(1) cont_saltrend(end)]);
hold on
colormap(h4(5), cmap_saltrend)
[k1 k2] = contour(d1, -meanZ, trendz_gam.*dSaldgam*1e3, ...
    [0 0], '-', 'LineColor', [0 0 0], 'LineWidth', 1.5);
clabel(k1, k2, 'labelspacing', 400, 'fontsize', 7, 'fontname', 'SansSerif');
axis(h4(5), [0 890 -760 0]);
xlabel('Distance from North [km]');
% ylabel('Depth [m]');
set(h4(5), 'TickLength', [0.025 0.025], 'TickDir', 'out', ... 
    'YTickLabel', []);

% Hatched regions
hat4 = trendz_gam;
hat4(sigtrz_mk_gam==0) = 10;
[~, gh4] = contourf(d1, -meanZ, hat4, [10 10]);
set(gh4, 'linestyle', 'none', 'Tag', 'HatchingRegion');
hp = findobj(h4(5), 'Tag', 'HatchingRegion');
hh = hatchfill2(hp, 'cross', 'HatchAngle', 45, 'HatchDensity', 60, ...
    'HatchColor', [0.5 0.5 0.5], 'HatchLineWidth', 0.5, 'Fill', 'off');
text(h4(5), 50, -700, 'Heave {\itS}', 'BackgroundColor', 'w');
[zz1 zz2] = contour(nanmax(dist_u, [], 2)', -z, ...
    nanmean(ubci(:, :, indx), 3)*1e2, cont_ugeo, ...
    'linecolor', colugeo, 'linestyle', '-');
clabel(zz1, zz2, 'labelspacing', 400, 'color', colugeo, 'fontsize', 6);
h4(5).Position = axPos(5, :);
p1 = plot([820 745], [0 0], '-', 'linewidth', 4, 'color', 'm');
p2 = plot([300 500], [0 0], '-', 'linewidth', 4, 'color', rgb('DarkGray'));
p3 = plot([35 135], [0 0], '-', 'linewidth', 4, 'color', rgb('Goldenrod'));


axes(h4(6))
contourf(nanmax(dist'), -z, ...
    trendtot_sal*1e3, cont_saltrend, 'edgecolor', 'none');
caxis(h4(6), [cont_saltrend(1) cont_saltrend(end)]);
hold on
colormap(h4(6), cmap_saltrend)
[k1 k2] = contour(nanmax(dist'), -z, ...
    trendtot_sal*1e3, [0 0], '-', ...
    'LineColor', [0 0 0], 'LineWidth', 1.5);
clabel(k1, k2, 'labelspacing', 400, 'fontsize', 7, 'fontname', 'SansSerif');
axis(h4(6), [0 890 -760 0]);
set(h4(6), 'TickLength', [0.025 0.025], 'TickDir', 'out', ... 
    'YTickLabel', []);
xlabel('Distance from North [km]');
col4 = colorbar(h4(6), 'Position', [sum(axPos(6, [1,3]))+0.01 ...
    axPos(6, 2) 0.01 axPos(6, 4)], 'Ticks', cont_saltrend(1:2:end)',...
    'TickLabels', num2str(cont_saltrend(1:2:end)'), 'TickLength', 0.025, ...
    'TickDirection', 'out');
ylabel(col4, 'Trend {\itS} \times 10^-^3 [year^-^1]');
cbarrow(col4, cmap_saltrend);


% Hatched regions
hat5= trendtot_sal;
hat5(sigtrtot_mk_sal==0) = 10;
[~, gh5] = contourf(nanmax(dist'), -z, hat5, [10 10]);
set(gh5, 'linestyle', 'none', 'Tag', 'HatchingRegion');
hp = findobj(h4(6), 'Tag', 'HatchingRegion');
hh = hatchfill2(hp, 'cross', 'HatchAngle', 45, 'HatchDensity', 60, ...
    'HatchColor', [0.5 0.5 0.5], 'HatchLineWidth', 0.5, 'Fill', 'off');
text(h4(6), 50, -700, 'Depth {\itS}', 'BackgroundColor', 'w');
[zz1 zz2] = contour(nanmax(dist_u, [], 2)', -z, ...
    nanmean(ubci(:, :, indx), 3)*1e2, cont_ugeo, ...
    'linecolor', colugeo, 'linestyle', '-');
clabel(zz1, zz2, 'labelspacing', 400, 'color', colugeo, 'fontsize', 6);
h4(6).Position = axPos(6, :);
p1 = plot([820 745], [0 0], '-', 'linewidth', 4, 'color', 'm');
p2 = plot([300 500], [0 0], '-', 'linewidth', 4, 'color', rgb('DarkGray'));
p3 = plot([35 135], [0 0], '-', 'linewidth', 4, 'color', rgb('Goldenrod'));


set(findall(gcf, '-property', 'FontSize'), 'FontSize', 7);
set(findall(gcf, '-property', 'FontName'), 'FontName', 'SansSerif');
set(findall(gcf, '-property', 'TickDirection'), 'TickDirection', 'out');
for i = 1:length(h4);
    annotation(fig1, 'textbox', ...
        [h4(i).Position(1) sum(h4(i).Position([2,4]))+0.02 0.0480 0.0353],...
        'String', 'SAF', 'fontsize', 7, 'fontname', 'SansSerif', ...
        'EdgeColor', 'none');
    annotation(fig1, 'textbox', ...
        [h4(i).Position(1)+0.085 sum(h4(i).Position([2,4]))+0.02 0.0480 0.0353],...
        'String', 'PF', 'fontsize', 7, 'fontname', 'SansSerif', ...
        'EdgeColor', 'none');
    annotation(fig1, 'textbox', ...
        [h4(i).Position(1)+0.17 sum(h4(i).Position([2,4]))+0.02 0.0480 0.0353],...
        'String', 'SACCF', 'fontsize', 7, 'fontname', 'SansSerif', ...
        'EdgeColor', 'none');
end
annotation(fig1, 'textbox', [h4(1).Position(1)-0.04 sum(h4(1).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'a', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig1, 'textbox', [h4(2).Position(1)-0.04 sum(h4(2).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'b', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig1, 'textbox', [h4(3).Position(1)-0.04 sum(h4(3).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'c', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig1, 'textbox', [h4(4).Position(1)-0.04 sum(h4(4).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'd', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig1, 'textbox', [h4(5).Position(1)-0.04 sum(h4(5).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'e', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
annotation(fig1, 'textbox', [h4(6).Position(1)-0.04 sum(h4(6).Position([2,4]))+0.01 ...
    0.0480 0.0353], 'String', 'f', 'fontsize', 7, 'fontname', 'SansSerif', ...
    'EdgeColor', 'none', 'fontweight', 'bold');
print('-dpdf', './FiguresPaper/fig6_trend_spiceheave.pdf', '-r500');

