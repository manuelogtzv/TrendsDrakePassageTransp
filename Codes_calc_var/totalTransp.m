% Program to calculate transport across Drake Passage as a function of
% distance. Corrects transport by adding or subtracting the cross-transect
% velocity offset to the measured velocity. 
%
% 2022/4/27 - Manuel O. Gutierrez Villanueva
% 2023/6/27 - Renamed some variables for consistency, removed data in the
% northern shelf that is deeper than the bottom topography

clear all; 
close all;

gsize = 25; %km
cont = [-1.0:0.1:1.0];%vel contours
contmaperr = [0.01 0.05 0.10 0.20 0.5 0.80 1.0];%mapping error contours
decorr_x = 50000; %meters distance
decorr_y = 300; %meters depth
error = 0.2;

% loads colormaps
load('velcmap.mat')

% datein = datenum(2005, 01, 01, 0, 0 ,0);
adcp = 'os38nb';
datein = datenum(1999, 01, 01, 0 ,0 ,0);
dec_nb150 = 0;

% Velocities
nomz = 760;% os38nb max z: 780 or 210 or 970
load(['../Datasets/lmgvel' adcp '_50_300.mat']);

% Finds south and northbound transects
indsouth = find(lmgdata.northbound == 1);
indnorth = find(lmgdata.northbound == 0);

% Adds the uoffs
load('../Datasets/trandsmiss_theta.mat');
lmgdata.u(:, :, indsouth) = lmgdata.u(:, :, indsouth) - uoffs_south;
lmgdata.u(:, :, indnorth) = lmgdata.u(:, :, indnorth) - uoffs_north;
arch = ['TotalTransport_' adcp 'maxz_' ...
        num2str(nomz, '%3.0f') '_2.mat'];
    


boxlim = [0    1200   -nomz   0];

% Loads bathymetry to determine which cruises are good for transport
b = load('../Datasets/drakembf_avgd'); 
b.lat = [b.lats; b.lat; b.latn]; 
b.lon = [b.lonw b.lon]; 
b.bathyf = [b.bathyfw [b.bathyfs; b.bathyf; b.bathyfn]];


% % Loads transect data
% load(['lmgvel' adcp '_50_300missalign.mat']);

% Remove data in the SA shelf that is deeper than bottom topography
lmgdata.u(lmgdata.z>694, 1, :) = NaN;

lmgdata.u = lmgdata.u(lmgdata.z<=nomz, :, :);
lmgdata.u_om = lmgdata.u_om(lmgdata.z<=nomz, :, :);
% lmgdata.dist = lmgdata.dist(:, 1:length(lmgdata.northbound));
% lmgdata.lon = lmgdata.lon(:, 1:length(lmgdata.northbound));
% lmgdata.lat = lmgdata.lat(:, 1:length(lmgdata.northbound));
lmgdata.z = lmgdata.z(lmgdata.z<=nomz);


% Finds south and northbound transects
    indsouth = find(lmgdata.northbound == 1);
    indnorth = find(lmgdata.northbound == 0);

for i = 1:length(indsouth);
    [lmgdata.dist(:, indsouth(i)) ii] = ...
        sort(abs(nanmax(lmgdata.dist(:, indsouth(i))) - ...
        lmgdata.dist(:, indsouth(i))) + 12.5);
    lmgdata.u(:, :, indsouth(i)) = lmgdata.u(:, ii, indsouth(i));
    lmgdata.u_om(:, :, indsouth(i)) = lmgdata.u_om(:, ii, indsouth(i));
    lmgdata.lon(:, indsouth(i)) = lmgdata.lon(ii, indsouth(i));
    lmgdata.lat(:, indsouth(i)) = lmgdata.lat(ii, indsouth(i));
end
    
dist = [lmgdata.dist(1):gsize:size(lmgdata.dist,1)*gsize ];
mask = lmgdata.dist./lmgdata.dist;

%%%%%% Calculates transport %%%%%%%%%%%
transp = nan(size(lmgdata.dist));
dz = abs(diff(lmgdata.z(2:3)));
dz2 = lmgdata.z(2)-dz/2;

for i = 1:size(lmgdata.u, 2)
    transp(i, :) = squeeze(nansum(lmgdata.u(2:end, i, :)*dz, 1)) ...
        + dz2*squeeze(lmgdata.u(1, i, :)); % m^2 /s
end
% transp = transp + dz2*squeeze(lmgdata.u(1,:,:));

cumtransp = cumsum(gsize*1e3*transp, 1, 'omitnan');%Sv
transpdist = transp*1e3*gsize; %Sv per unit area
cumtransp = cumtransp.*mask;

% % Checks that takes only when the max is located
% for ii = 1:length(lmgdata.time);
%     indmax = find(cumtransp(:, ii) == nanmax(cumtransp(:, ii)));
%     cumtransp(indmax(1)+1:end, ii) = NaN;
%     transpdist(indmax(1)+1:end, ii) = NaN;
% end

% Total transport
tottransp = nansum(transpdist)';
[allMean allstd alldof] = nanmstd(tottransp);

% Removes outliers (transport is larger than 3*std dev)
an = find(abs(tottransp(indnorth)-nanmean(tottransp(indnorth)))>...
    3*nanstd(tottransp(indnorth)));
as = find(abs(tottransp(indsouth)-nanmean(tottransp(indsouth)))>...
    3*nanstd(tottransp(indsouth)));

tottransp(indnorth(an)) = NaN;
tottransp(indsouth(as)) = NaN;

% Transport and times for north and southbound transects
northtransp = tottransp(indnorth); timenorth = lmgdata.time(indnorth);
southtransp = tottransp(indsouth); timesouth = lmgdata.time(indsouth);

%%%%%%% Least-squares %%%%%%%
% All transects
phianualall = 2*pi*[lmgdata.time(~isnan(tottransp))-datein]/365.25;
phisemiall = 2*phianualall;

% Northbound
phianualnorth = 2*pi*[timenorth(~isnan(northtransp))-datein]/365.25;
phiseminorth = 2*phianualnorth;

% Southbound
phianualsouth = 2*pi*[timesouth(~isnan(southtransp))-datein]/365.25;
phisemisouth = 2*phianualsouth;


% Model y_m = a_0 + a_1*t + a_2*sin(omega_anual*t) + a_3*cos(omega_anual*t)
% + a_4*sin(omega_anual*t) + a_5*cos(omega_anual*t)
northmodel = [ones(length(phianualnorth), 1) ...
    [timenorth(~isnan(northtransp))-datein]'/365.25 ...
    cos(phianualnorth(:)) sin(phianualnorth(:)) cos(phiseminorth(:)) ...
    sin(phiseminorth(:))]; %northbound
southmodel = [ones(length(phianualsouth), 1) ...
    [timesouth(~isnan(southtransp))-datein]'/365.25 ...
    cos(phianualsouth(:)) sin(phianualsouth(:)) cos(phisemisouth(:)) ...
    sin(phisemisouth(:))]; %southbound
allmodel = [ones(length(phianualall), 1) ...
    [lmgdata.time(~isnan(tottransp))-datein]'/365.25 ...
    cos(phianualall(:)) sin(phianualall(:)) cos(phisemiall(:)) ...
    sin(phisemiall(:))]; %all transects

% Least-squares fit (coefficients and errors)
[northreg bnorth] = ...
    regress(northtransp(~isnan(northtransp))*1e-6, northmodel);
[southreg bsouth] = ...
    regress(southtransp(~isnan(southtransp))*1e-6, southmodel);
[allreg ball] = ...
    regress(tottransp(~isnan(tottransp))*1e-6, allmodel);

allModeledTransp = allreg(1) + allreg(2)*[0:0.05:21] + ...
    allreg(3)*cos(2*pi*[0:0.05:21]) + ...
    allreg(4)*sin(4*pi*[0:0.05:21]) + ...
    allreg(5)*cos(2*pi*[0:0.05:21]) + ...
    allreg(6)*sin(4*pi*[0:0.05:21]);

% Arithmetic means
% Northbound transports
[northMean northstd northdof] = nanmstd(tottransp(indnorth));
northMedian = nanmedian(tottransp(indnorth));


[southMean southstd southdof] = nanmstd(tottransp(indsouth));
southMedian = nanmedian(tottransp(indsouth));

% Percentiles (uncorrected)
Ynorth = prctile(tottransp(indnorth)*1e-6, [2.5 97.5]);
Ysouth = prctile(tottransp(indsouth)*1e-6, [2.5 97.5]);
Yall = prctile(tottransp*1e-6, [2.5 97.5]);

totaltransp.adcp = adcp;
totaltransp.lon = lmgdata.lon;
totaltransp.lat = lmgdata.lat;
totaltransp.dist = dist;
totaltransp.nettransp = tottransp';
totaltransp.cumtransp = cumtransp;
totaltransp.transp = transpdist;
totaltransp.time = lmgdata.time;
totaltransp.north = lmgdata.northbound';
totaltransp.maxdepth = nomz;
totaltransp.doc = {'Total transport corrected for misalignment angle';...
               'adcp is the adcp used';...
               'transp is the  transport [m^3 s^{-1}] per distance bin';...
               'cumtransp is the cumulative transport [m^3 s^{-1}]';...
               'nettransp is the net Drake Passage transport [m^3 s^{-1}]';...
               'dist is the alongtrack distance from the north';...
               'time is the time';...
               'north = 1 is a northbound transect';...
               'lon and lat are the longitude and latitude';...
               'maxdepth is the deepest depth used [m]'};

% save(arch, 'totaltransp');

