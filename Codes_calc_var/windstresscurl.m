% Code to estimate wind stress, wind stress curl and trends in the wind
% stress curl. Uses either ERA5 (ECMWF) or MERRA2 (NASA GDAC). For
% estimating the trends, two approaches are used: a) least-squares fit, and
% b) nonparemetric test (Theil-Sen estimator). Before computing the trends,
% seasonal cycle from the monthly means is removed. Confidence intervals
% for the least-squares are calculated following Fyne and Mckinley.
% Modified Mann-Kendall test is employed for the nonparametric trend.
%
%
% Created by Manuel O. Gutierrez-Villannueva
% 2022/07/22

clear all; close all;


datein = datenum(1996, 01, 01, 0, 0 ,0);
alpha = 0.05;
season_opt = 1;


%%%%% ERA5 %%%%%%
fprintf('\n Loading ERA 5 data\n');
lonera = ncread('../BigDatasets/winddata_ERA5_wind10m_monthly.nc', 'longitude');
latera = ncread('../BigDatasets/winddata_ERA5_wind10m_monthly.nc', 'latitude');
timeera = ncread('../BigDatasets/winddata_ERA5_wind10m_monthly.nc', 'time');
% initimeera = datenum(1900, 01, 01, 0 , 0 ,0);
% timeera = initimeera + hours(timeera);
uera = ncread('../BigDatasets/winddata_ERA5_wind10m_monthly.nc', 'u10');
vera = ncread('../BigDatasets/winddata_ERA5_wind10m_monthly.nc', 'v10');
%uera = squeeze(uera(:, :, 1, :));
%vera = squeeze(vera(:, :, 1, :));
timeera = datenum(1900, 1, 1, double(timeera), 0, 0);

uera = uera(:, :, timeera>=datenum(1996, 1, 1) & timeera<datenum(2020, 1, 1));
vera = vera(:, :, timeera>=datenum(1996, 1, 1) & timeera<datenum(2020, 1, 1));
timeera = timeera(timeera>=datenum(1996, 1, 1) & timeera<datenum(2020, 1, 1));

% Calculates wind stress
Tx_era = nan(size(uera));
Ty_era = nan(size(vera));
StrCurl_era = nan(size(uera));


[llat_era llon_era] = meshgrid(latera, lonera);
% [llat_merra llon_merra] = meshgrid(latmerra, lonmerra); 

for i = 1:size(Tx_era, 3);
    % Wind stress
    [Tx_era(:, :, i), Ty_era(:, :, i)] = ra_windstr(squeeze(uera(:, :, i)), ...
        squeeze(vera(:, :, i)));
%     [Tx_merra(:, :, i), Ty_merra(:, :, i)] = ra_windstr(squeeze(umerra(:, :, i)), ...
%         squeeze(vmerra(:, :, i)));

    % Wind stress curl
    StrCurl_era(:, :, i) = cdtcurl(llat_era, llon_era, ...
        squeeze(Tx_era(:, :, i)), squeeze(Ty_era(:, :, i)));
%     StrCurl_merra(:, :, i) = cdtcurl(llat_merra, llon_merra, ...
%         squeeze(Tx_merra(:, :, i)), squeeze(Ty_merra(:, :, i)));
end



% Calculate trends
% ERA5 %
tr_curl_era = nan(size(squeeze(Tx_era(:, :, 1))));
ci_curl_era = tr_curl_era;
tr_np_curl_era = tr_curl_era;
sig_np_curl_era = tr_curl_era;

tr_tx_era = tr_curl_era;
ci_tx_era = tr_curl_era;
tr_np_tx_era = tr_curl_era;
sig_np_tx_era = tr_curl_era;

StrCurl_era = reshape(StrCurl_era, ...
    [size(StrCurl_era, 1)*size(StrCurl_era, 2) size(StrCurl_era, 3)]);

% for i = 1:size(StrCurl_era, 1)
%     for j = 1:size(StrCurl_era, 2)
        [tr_curl_era, ci_curl_era, ~, tr_np_curl_era sig_np_curl_era] = ...
            calctrends(StrCurl_era*1e6, ...
            [timeera-datein]'/365.25, alpha, season_opt);

%         [tr_tx_era, ci_tx_era, ~, tr_np_tx_era sig_np_tx_era] = ...
%             calctrends2D(Tx_era(i, j, :)', ...
%             [timeera-datein]/365.25, alpha, season_opt);
%     end
% end


StrCurl_era = reshape(StrCurl_era, ...
    [size(llon_era, 1), size(llon_era, 2), size(timeera, 1)]);
timeera = timeera';

tr_curl_era = reshape(tr_curl_era, size(llat_era, 1), ...
    size(llat_era, 2));
ci_curl_era = reshape(ci_curl_era, size(llat_era, 1), ...
    size(llat_era, 2));
tr_np_curl_era = reshape(tr_np_curl_era, size(llat_era, 1), ...
    size(llat_era, 2));
sig_np_curl_era = reshape(sig_np_curl_era, size(llat_era, 1), ...
    size(llat_era, 2));
doc = {'tr_curl_era and ci_curl_era are the trends and confidence intervals';...
       'using the least-squares fit approach (seasonal cycle removed).';...
       'tr_np_curl_era and sig_np_curl_era are the trends and significance';...
       'test (modified Mann-Kendall test). Significant if test is = 1.';...
       'Trends are estimated using the nonparametric Theil-Sen estimator.';...
       'lonera and latera are the longitude and latitude vectors';...
       'StrCurl_era is the wind stress curl time series';...
       'Tx and Ty are the zonal and meridional windstress time series.'};

save('TrendsERA5_1996_2019_monthly.mat', 'lonera', 'latera', 'StrCurl_era', ...
    'tr_curl_era', 'ci_curl_era', 'tr_np_curl_era', 'sig_np_curl_era', ...
    'timeera', 'Tx_era', 'Ty_era', 'doc');




