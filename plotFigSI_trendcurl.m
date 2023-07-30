% Code to estimate and plot the mean wind stress curl and trends in the
% Southern Ocean (<40ºS). Wind stress curl is estimated using the daily 
% outputs European Center for Medium-Range Weather Forecast (ECMWF) ERA5 
% Reanalysis. Trends are estimated using least-squares after removing the
% seasonal cycle. Trends significance test is performed using the modified
% Mann-Kendall test.
%
% Manuel O. Gutierrez Villanueva 2022/07/22

clear all;
close all;

pycmd = '/Users/manuelgutierrez/anaconda3/bin/python';


% Loads data
load TrendsERA5_1996_2019_monthly.mat


curlcont = [-11:11];
trcurlcont = [-0.1:0.01:0.1];
cmapcurl = getPyPlot_cMap('BrBG', length(curlcont) - 1, [], pycmd);
cmapcurltrend = getPyPlot_cMap('BrBG', length(trcurlcont) - 1, [], pycmd);
alpha = 0.5;
datein = datenum(1996, 1, 1);

[llat_era llon_era] = meshgrid(latera, lonera);
% [llat_merra llon_merra] = meshgrid(latmerra, lonmerra);

% StrCurl_era = reshape(StrCurl_era, ...
%     [size(llon_era, 1), size(llon_era, 2), size(timeera, 2)]);

% StrCurl_merra = reshape(StrCurl_merra, ...
%     [size(llon_merra, 1), size(llon_merra, 2), size(timeera, 2)]);

% time series or zonally averaged windstress curl in Drake Passage
latDP = [-62.5 -55];
lonDP = [-70 -56];

drake_lat = latera(latera>= latDP(1) & latera <= latDP(end));
drake_lon = lonera(lonera>= lonDP(1) & lonera<= lonDP(end));

% Period of study
ind = timeera>=datenum(2005, 10, 1) & timeera<=datenum(2019, 12, 31);
timeera = timeera(ind);
StrCurl_era = StrCurl_era(:, :, ind);

[drake_time drake_llat] = meshgrid(timeera, drake_lat);

drake_curl = StrCurl_era(lonera>= lonDP(1) & lonera<= lonDP(end), ...
    latera>= latDP(1) & latera <= latDP(end), :);
weights = sqrt(cosd(repmat(drake_lat', ...
    [size(drake_curl, 1) 1 size(drake_curl, 3)])));
avrdrake_curl = squeeze(nanmean(drake_curl, 1));

for i = 1:length(drake_lat);
    [~, avrDPcurl_noseas(i, :)] = fitseasoncycle(avrdrake_curl(i, :)*1e8, ...
        [timeera - datein]/365.25);
end

[tr_curl, int_curl, ~, trnp_curl signp_curl] = ...
    calctrends(avrDPcurl_noseas, [timeera - datein]/365.25, alpha);

figure('color', 'w')
figure_width = 17;
figure_height = 8;
% --- setup plot windows
figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
set(gcf,'Visible', figuresVisible)
set(gcf, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Renderer','opengl');

ha = tight_subplot(1, 2, [0.1 0.1], [0.15 0.15], [0.15 0.15]);

axes(ha(1));
m_proj('stereographic', 'lat', -90, 'long', -60, 'radius',40);
m_contourf(lonera, latera, squeeze(nanmean(StrCurl_era*1e8, 3))', ...
    [-225:225], 'edgecolor', 'none');
hold on;
caxis(ha(1), [curlcont(1) curlcont(end)]);
colormap(ha(1), cmapcurl)
m_coast('patch', [0.5 0.5 0.5], 'edgecolor', 'none');
m_grid('xtick', 12, 'tickdir', 'out', 'ytick', [-80 -70 -60 -50 -40],...
    'xtick', [-180:40:180], 'linest', ':', 'xaxislocation', 'top', ...
    'yticklabel', []);
tit = title('a) Mean {$\mathbf{k}$ \boldmath$\cdot$}{\boldmath${\nabla}$} $\times$ {\boldmath${\tau}$} (2005 - 2019)', ...
    'interpreter', 'latex');
tit.Position(2) = 0.9354;
col2 = colorbar; 
col2.Location = 'SouthOutside';
col2.Position(4) = 0.01; 
col2.Position(2) = col2.Position(2)-0.13;
col2.Ticks = curlcont(2:2:end)'; col2.TickLabels = num2str(curlcont(2:2:end)');
cbarrow(col2, cmapcurl);
ylabel(col2, ' $10^{-8}$ [N m$^{-3}$]',...
    'interpreter', 'latex', 'fontsize', 10);

axes(ha(2));
m_proj('stereographic', 'lat', -90, 'long', -60, 'radius',40);
m_contourf(lonera, latera, tr_curl_era'*1e2, [-4:0.01:4], ...
    'edgecolor', 'none');
hold on;
caxis(ha(2), [trcurlcont(1) trcurlcont(end)]);
colormap(ha(2), cmapcurltrend);
hat2 = tr_curl_era*1e4;
hat2(sig_np_curl_era==0) = 0;
[~, gh2] = m_contour(lonera, latera, hat2', [0 0], 'k', 'linewidth', 0.5);
% m_2 = m_plot([lonDP(1) lonDP(2) lonDP(2) lonDP(1) lonDP(1)], ...
%     [latDP(1), latDP(1), latDP(2), latDP(2), latDP(1)], '-b', ...
%     'linewidth', 2);
% m_3 = m_contour(lon_niimax, lat_niimax', mdt_niimax', [psi_saf psi_saf], ...
%     'color', 'c', 'linewidth', 1);
% m_4 = m_contour(lon_niimax, lat_niimax', mdt_niimax', [psi_pf psi_pf], ...
%     'color', 'c', 'linewidth', 1);
% m_5 = m_contour(lon_niimax, lat_niimax', mdt_niimax', [psi_saccf psi_saccf]-.025, ...
%     'color', [0.3 0.3 0.3], 'linewidth', 1);
m_coast('patch', [0.5 0.5 0.5], 'edgecolor', 'none');
m_grid('xtick', 12, 'tickdir', 'out', 'ytick', [-80 -70 -60 -50 -40],...
    'xtick', [-180:40:180], 'linest', ':', 'xaxislocation', 'top',...
    'yticklabel', []);
tit = title('b) Trend {$\mathbf{k}$ \boldmath$\cdot$}{\boldmath${\nabla}$} $\times$ {\boldmath${\tau}$} (2005 - 2019)', ...
    'interpreter', 'latex');
tit.Position(2) = 0.9354;
col1 = colorbar; 
col1.Location = 'SouthOutside';
col1.Position(4) = 0.01; 
col1.Position(2) = col1.Position(2)-0.13;
col1.Ticks = trcurlcont(1:2:end)'; 
col1.TickLabels = num2str(trcurlcont(1:2:end)');
cbarrow(col1, cmapcurltrend);
ylabel(col1, '$10^{-10}$ [N m$^{-3}$ year$^{-1}$]',...
    'interpreter', 'latex', 'fontsize', 10);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 6);
set(findall(gcf, '-property', 'FontName'), 'FontName', 'SansSerif');
set(findall(gcf, '-property', 'TickDirection'), 'TickDirection', 'out');
print('-dpdf', ['./FiguresPaper/figSI_trendscurl.pdf'], '-r500');