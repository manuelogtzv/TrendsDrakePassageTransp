% Code for estimating the difference between total transport and baroclinic
% transport referenced at 760 m. Barotropic component is then calculated as
% the depth-averaged residual velocity (total - baroclinic) at each
% distance position and transect, multiplied by water depth and distance.
% 
% Manuel O. Gutierrez Villanueva
% 
% 2021/11/29
%
% 2021/12/14 - Includes nb150 option
% 2023/06/26 - No longer includes the deep estimate, changed the name of
% the output .mat file, saves velocity and transport in a single .mat file.
% 2023/06/27 - Changed the interpolation method. Now uses interp1 instead
% of griddata.

clear all;
close all;

adcp = 'os38nb';
cont_bar = [-10:2:10]; 
cont_bci = [-4:0.5:4];
cont_tot = [-30:5:30];
datein = datenum(1996, 1, 1);
ptch = [0 0 25 25 0; 1000 -10 -10 1000 1000];
bxlim = [0 25];

% Loads colormaps
load('./geos_ref_cmap.mat');



% Loads baroclinic velocities
load(['geostrophic_noaver.mat']);

% Velocities
nomz = 760;% os38nb max z: 780 or 210
load(['lmgvel' adcp '_50_300.mat']);
    
% loads total trasnport data from adcp
load(['TotalTransport_' adcp 'maxz_' num2str(nomz) '_2.mat']);
lmgtransp = totaltransp; clear totaltransp;

% Keeps nominal depth
lmgdata.u = lmgdata.u(lmgdata.z<=760, :, :);
lmgdata.u_om = lmgdata.u_om(lmgdata.z<=760, :, :);
lmgdata.z = lmgdata.z(lmgdata.z<=760);


if strcmp(adcp, 'os38nb');

    % Load velocity offset
    load('trandsmiss_theta.mat')
    
%     % Velocity offset for os38nb
    indsouth = find(lmgdata.northbound == 1);
    indnorth = find(lmgdata.northbound == 0);

    % Adds the uoffs
    load('trandsmiss_theta.mat');
    lmgdata.u(:, :, indsouth) = lmgdata.u(:, :, indsouth) - uoffs_south;
    lmgdata.u(:, :, indnorth) = lmgdata.u(:, :, indnorth) - uoffs_north;

    nomz = 760;
    lmgdata.u = lmgdata.u(lmgdata.z<=nomz, :, :);
    lmgdata.u_om = lmgdata.u_om(lmgdata.z<=nomz, :, :);
    lmgdata.z = lmgdata.z(lmgdata.z<=nomz);


    % Makes the transect that start from south to north to now be north to
    % south
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


end



r = 0;

ubc = nan(length(lmgdata.z), ...
    size(geosbaroc.u, 2), size(geosbaroc.u, 3));

% Interpolates baroclinic velocities profiles to adcp depths
for i = 1:size(geosbaroc.u, 2);
    for j = 1:size(geosbaroc.u, 3);
        
        if ~isnan(nanmean(squeeze(geosbaroc.u(:, i, j))));
            ubc(:, i, j) = interp1(geosbaroc.z(:), ...
                squeeze(geosbaroc.u(:, i, j)), lmgdata.z, 'pchip', ...
                'extrap');
        end
    end
end


% Finds coincident xbt and adcp transects
[x y] = meshgrid(geosbaroc.time, lmgdata.time);
[t1 t2] = meshgrid(1:length(geosbaroc.time), 1:length(lmgdata.time));

ind_lmg = t2(abs(x-y)<=5);
ind_geos = t1(abs(x-y)<=5);

% % Only coincident transects
distubci = geosbaroc.dist(:, ind_geos);% Cumulative along-track distance
diffdist = diff(geosbaroc.dist(:, ind_geos));
distubci = geosbaroc.dist_u(:, ind_geos);


clear x y t1 t2

figure('color', 'w');
h3 = tight_subplot(5, 3, [0.01 0.01], [0.10 0.07], [0.15 0.03]);
% ubar1 = nan(length(lmgdata.z), size(lmgdata.lon, 1), ...
%     length(ind_geos));

r = 0;

fig2 = figure;

% Substract baroclinic velocity from total velocity
for i = 1:length(lmgdata.z);
    for j = 1:length(ind_lmg);

        utot = squeeze(lmgdata.u(i, :, ind_lmg(j)))';
        ubci = squeeze(ubc(i, :, ind_geos(j)));

        lontot = lmgdata.lon(:, ind_lmg(j)) - 360;
        lattot = lmgdata.lat(:, ind_lmg(j));
        

        % Applies mask to data
        lonbci = geosbaroc.lon_u(:, ind_geos(j));
        latbci = geosbaroc.lat_u(:, ind_geos(j));
        mask = geosbaroc.lon_u(:, ind_geos(j))./...
            geosbaroc.lon_u(:, ind_geos(j));


        vtot = nan(size(lonbci));

        % Interpolates total velocity to baroclinic velocity grid
%         vtot(~isnan(latbci)) = griddata(lontot(~isnan(lontot)), ...
%             lattot(~isnan(lontot)), utot(~isnan(lontot)), ...
%             lonbci(~isnan(latbci)), latbci(~isnan(latbci)), 'linear');
        vtot(~isnan(latbci)) = interp1(lattot(~isnan(lontot)), utot(~isnan(lontot)), ...
            latbci(~isnan(latbci)), 'linear');
        vtot = vtot;


        ubci1(i, :, j) = ubci;
        utot1(i, :, j) = vtot;
        ubar1(i, :, j) = vtot(:) - ubci(:);

    end

    vbar1 = squeeze(ubar1(i, :, :));
    vbci1 = squeeze(ubci1(i, :, :));
    vtot1 = squeeze(utot1(i, :, :));
    disttot = lmgdata.dist(:, ind_lmg);
    
    
    if mod(i, 7) == 1 & r<5;
        
        r = r + 1;
        
        figure(gcf);
        axes(h3((r*3)-2));
        scatter(reshape(repmat([lmgdata.time(ind_lmg) - ...
            datenum(1996, 1, 1)]/365.25, ...
            size(distubci, 1), 1), size(vbar1(:))), ...
            reshape(distubci, size(vbar1(:))), 8,...
            vbar1(:), 'o', 'filled');
        colormap(h3((r*3)-2), cmap_tot);
        caxis(h3((r*3)-2), [-1.0 1.0]);
        text(h3((r*3)-2), 12, 900, [num2str(lmgdata.z(i)) ' m']);
        box on;
        
        axes(h3((r*3)-1));
        scatter(reshape(repmat([geosbaroc.time(ind_geos) - ...
            datenum(1996, 1, 1)]/365.25, ...
            size(distubci, 1), 1), size(distubci(:))), distubci(:), 8,...
            vbci1(:), 'o', 'filled');
        colormap(h3((r*3)-1), cmap_tot);
        caxis(h3((r*3)-1), [-1.0 1.0]);
        text(h3((r*3)-1), 12, 900, [num2str(lmgdata.z(i)) ' m']);
        box on;

        axes(h3((r*3)));
        scatter(reshape(repmat([lmgdata.time(ind_lmg) - ...
            datenum(1996, 1, 1)]/365.25, ...
            size(distubci, 1), 1), size(distubci(:))), ...
            reshape(distubci, size(vbar1(:))), 8,...
            vtot1(:), 'o', 'filled');
        colormap(h3((r*3)), cmap_tot);
        caxis(h3((r*3)), [-1.0 1.0]);
        text(h3((r*3)), 12, 900, [num2str(lmgdata.z(i)) ' m']);
        box on;
    end
    
    
end



title(h3(1), 'u_r_e_f [m/s]');
title(h3(2), 'u_g_e_o_s [m/s]');
title(h3(3), 'u_t_o_t [m/s]');
ylabel(h3(10), 'Distance from north [km]');
xlabel([h3(end-2)],...
    ['Years since ' datestr(datenum(1996, 1, 1))]);
xlabel([h3(end-1)],...
    ['Years since ' datestr(datenum(1996, 1, 1))]);
xlabel([h3(end)],...
    ['Years since ' datestr(datenum(1996, 1, 1))]);



xx = squeeze(lmgdata.u(:, :, ind_lmg));
x1 = squeeze(lmgdata.u(:, :, :));

% dist_bar = lmgdata.dist(:, ind_lmg);

% Calculates reference velocity
uref1 = nanmean(ubar1(lmgdata.z>=100 & lmgdata.z<=nomz, :, :), 1);
dprng = lmgdata.z(lmgdata.z>=100 & lmgdata.z<=nomz);
dz = abs(diff(lmgdata.z(2:3)));
dldp = (dprng(end) + dz/2) - (dprng(1) - dz/2);

mask =ubar1(1, :, :)./ubar1(1, :, :);

% Calculates transport due to the barotropic component
transp_bar = squeeze(uref1).*diffdist*dldp*1e3;
cumtransp_bar = cumsum(transp_bar, 1, 'omitnan').*squeeze(mask);
tottransp_bar = nansum(transp_bar, 1);


dzs = lmgdata.z(2)-diff(abs(lmgdata.z(2:3)))/2;
dz = abs(diff(lmgdata.z(2:3)));

% Calculates Ekman transport
uEks = squeeze(ubar1(1, :, :)).*diffdist*1e3*dzs;
uEk = squeeze(nansum(ubar1(lmgdata.z>0 & lmgdata.z<100, :, :), 1)).*diffdist*1e3*dz;
svEk = uEks + uEk;
Ektransp = nansum(svEk, 1);


geosref.adcp = adcp;
geosref.nettransp = tottransp_bar;
geosref.transp = transp_bar;
geosref.cumtransp = cumtransp_bar;
geosref.u = uref1;
geosref.ektransp = svEk;
geosref.eknettransp = Ektransp;
geosref.lon = geosbaroc.lon_u(:, ind_geos);
geosref.lat = geosbaroc.lat_u(:, ind_geos);
geosref.dist = geosbaroc.dist_u(:, ind_geos);
geosref.time = geosbaroc.time(ind_geos);
geosref.doc = {'Reference velocity and transport';...
               'adcp is the adcp used';...
               'transp is the transport [m^3 s^{-1}] per distance bin';...
               'cumtransp is the cumulative transport [m^3 s^{-1}]';...
               'nettransp is the net Drake Passage transport [m^3 s^{-1}]';...
               'ektransp is the Ekman transport [m^3 s^{-1}] per distance bin';...
               'eknettransp is the net Drake Passage Ekman transport [m^3 s^{-1}]';...
               'u is the barotropic velocity component [m s^{-1}]';...
               'dist is the alongtrack distance from the north [km]';...
               'time is the time';...
               'lon and lat are the longitude and latitude'};
           
save(['reference_' adcp '_760.mat'], 'geosref');


