function [sm3_om] = obj_map_3D_sADCP(simple_mean_3D,decorr_scale,var,error,plotm);

%%% Function that calculates the objectively mapped mean geostrophic
%%% streamfunction (Gille 2003a; Firing et al. 2014) assuming a gaussian
%%% covariance matrix for the streamfunction. Error, variance and
%%% decorrelation scales are inputs. The 3D simple-mean velocities are
%%% also an input (as a structure variable). To calculate the anomalies, we
%%% use the Maximenko et al. (2012) mean dynamic topography (MDT). If
%%% plotm = 1, then the mean streamfunction is calculated.
%%% 
%%% Inputs:
%%% + simple_mean_3D: structure variable containing the simple means of
%%% velocity.
%%% + decorr_scale: decorrelation scales [Lx Ly] in meters.
%%% + var: varianza of the signal. If var = 1, then it is assumed as a
%%% correlation coefficient.
%%% + error: noise-to-signal ratio.
%%% + plotm: if specified, plots the geostrophic streamfunction for some
%%% depth bins.
%%%
%%% 08/19/2019: Saves also lat0,lon0,theta0, which come from the
%%% simple_mean_3D structure file.

path = '/Users/manuelgutierrez/Desktop/Thesis';
path1 = [path '/drake_mean/data/all_drake/'];%all data gridded and transects

load([path1 'DPshort.mat'],'theta');
lat0 = simple_mean_3D.lat0;
lon0 = simple_mean_3D.lon0;
theta0 = simple_mean_3D.theta0;
omega = 2*pi/86400;

[cc aa bb] = size(simple_mean_3D.u);
[X Y] = meshgrid(simple_mean_3D.x,simple_mean_3D.y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% NIILER and MAXIMENKO DYNAMIC TOPOGRAPHY %%%%%%%
nii_mdot = ...ncread('cnes_cls13_global.nc','mdt');
    ncread('../niiler_maximenko_1992_2012_mdt.nc','MDOT')*1e-2;
nii_lon = ...ncread('cnes_cls13_global.nc','lon') - 360;%
    ncread('../niiler_maximenko_1992_2012_mdt.nc','LONN319_N207');
nii_lat = ...ncread('cnes_cls13_global.nc','lat');%
    ncread('../niiler_maximenko_1992_2012_mdt.nc','LAT37_109');

[r_nii c_nii] = size(nii_mdot);

nii_lon = repmat(nii_lon,1,c_nii);
nii_lat = repmat(nii_lat',r_nii,1);

% Mask of 1000 m depth
% load mask_DP_1000m
% in = inpolygon(nii_lon,nii_lat,lon_mask1000,lat_mask1000);
% nii_mdot(~in) = nan;


% compute grid positions
[xr_nii, yr_nii] = ll2xy2(nii_lat(:),nii_lon(:),lat0,lon0,'km'); 

xx_nii = reshape(xr_nii,r_nii,c_nii);
yy_nii = reshape(yr_nii,r_nii,c_nii);

% Rotates CCW
W_nii = complex(xx_nii,yy_nii).*exp(-sqrt(-1)*theta*pi/180);
xr_nii = real(W_nii); yr_nii = imag(W_nii);

% Coriolis parameter and gravity
grav = gsw_grav(nii_lat,0);% m/s
f = 2*omega*sin(nii_lat*pi/180);

% Lens and YLF code
gsuv_nii = ssh2gv(nii_lon(:,1),nii_lat(1,:)',nii_mdot'*100); 
gsuv_nii = permute(gsuv_nii,[2 1]);
Wr_nii = gsuv_nii.*exp(-sqrt(-1)*theta*pi/180);
ur_nii = real(Wr_nii); vr_nii = imag(Wr_nii);

%%%%% Now we start by doing simple interpolation of velocities from
%%%%% Maximenko and Niiler MDT %%%%%
U_float_nii = interp2(nii_lon',...
    nii_lat',real(gsuv_nii)',simple_mean_3D.lon,simple_mean_3D.lat);
V_float_nii = interp2(nii_lon',...
    nii_lat',imag(gsuv_nii)',simple_mean_3D.lon,simple_mean_3D.lat);
mdot_argo_nii = interp2(nii_lon',...
    nii_lat',nii_mdot',simple_mean_3D.lon,simple_mean_3D.lat);

% Rotates
W_float_nii = complex(U_float_nii,V_float_nii);
Wr_float_nii = W_float_nii.*exp(-sqrt(-1)*theta*pi/180);%rotates
Ur_float_nii = real(Wr_float_nii); Vr_float_nii = imag(Wr_float_nii);


% Pre-alocate memory
om_mean_psi = nan(cc,aa,bb);
om_mean_u = om_mean_psi;
om_mean_v = om_mean_psi;

om_prime_psi = om_mean_psi;
om_prime_u = om_mean_psi;
om_prime_v = om_mean_psi;

om_Error = om_prime_u;
XXint = om_prime_u;
YYint = om_prime_v;

res_Ur_om = om_prime_u;
res_Vr_om = om_prime_u;
angle_res = om_prime_v;
norm_vec = om_prime_v + sqrt(-1).*om_prime_v;


tic

% Axes and colorbar position
pos1 = [0.0310    0.5014    0.4337    0.3266];
pos2 = [0.4561    0.5014    0.4466    0.3261];
pos3 = [0.0310    0.1000    0.4337    0.3266];
pos4 = [0.4561    0.1026    0.4466    0.3261];
colb1 = [0.8382    0.5000    0.0100    0.3272];
colb2 = [0.4068    0.0991    0.0100    0.3272];
colb3 = [0.8382    0.1014    0.0100    0.3272];


% Legend position
leg1 = [0.1194    0.6917    0.0787    0.0262];
leg2 = [0.5451    0.6618    0.0936    0.0572];
leg3 = [0.1151    0.2648    0.0936    0.0507];
leg4 = [0.5452    0.2715    0.1000    0.0463];


    for m = 1:cc;

        nanfac = squeeze(real(simple_mean_3D.u(m,:,:))./...
            real(simple_mean_3D.u(m,:,:)));

        X1 = X.*nanfac; Y1 = Y.*nanfac; %makes nans where no data available
    
        uui = squeeze(simple_mean_3D.u(m,:,:));
        vvi = squeeze(simple_mean_3D.v(m,:,:));

        % Calculates anomalies
        Wr_anom_upd = complex(uui,vvi) - Wr_float_nii;
        Ur_anom_upd = real(Wr_anom_upd); 
        Vr_anom_upd = imag(Wr_anom_upd);


        %%%%%%%%%%%% Objective mapping %%%%%%%%%%
        f_sm3 = 2*omega*sind(simple_mean_3D.lat);
        grav_sm3 = gsw_grav(simple_mean_3D.lat,simple_mean_3D.z(m));

        fu = Ur_anom_upd(:).*f_sm3(:); 
        fv = Vr_anom_upd(:).*f_sm3(:); 

        xu = X(:);
        yv = Y(:);

        llon = simple_mean_3D.lon(:); llon = llon(~isnan(fu));
        llat = simple_mean_3D.lat(:); llat = llat(~isnan(fu));

        % Removes nans
        xu = xu(~isnan(fu));
        yv = yv(~isnan(fv));
 
        fu = fu(~isnan(fu));
        fv = fv(~isnan(fv));


        % New parameters for objective mapping
%         boxsize = 25000;
        L = decorr_scale;
        C = var;
        err = error;


        % Mapping Psi, f*u_vec 
        [psi umap vmap E] = obj_mapping_geos(xu(:),yv(:),fu,fv,...
            X(:),Y(:),err,L,C);

        % Convert to matrix values and Errors
        psi = reshape(psi,aa,bb);
        umap = reshape(umap,aa,bb);
        vmap = reshape(vmap,aa,bb);
        % z1 = reshape(z(:,1)',aa,bb);


        Epsi = reshape(E(:,1),aa,bb);
        Eumap = reshape(E(:,2),aa,bb);
        Evmap = reshape(E(:,2),aa,bb);

        % psi(mean_noverlap.num < 2) = 0;
        % umap(mean_noverlap.num < 2) = 0;
        % vmap(mean_noverlap.num < 2) = 0;


        % Convert to psi (m), uprime and vprime (m/s)
        MDT_prime_om = psi./grav_sm3;
        U_prime_om = umap./f_sm3;
        V_prime_om = vmap./f_sm3;
        % z1 = z1./grav_sm3;

        Xint = X; Yint = Y;
        mdot_nii_mean = mdot_argo_nii;
        ur_nii_mean = Ur_float_nii;
        vr_nii_mean = Vr_float_nii;

        % make nans outside of the area enclosed by 3 most 
        % repeated transects
        for i = 1:aa;
            iind = find(isnan(nanfac(i,:)) == 0);
    
            if ~isempty(iind)
                ind1 = [1:iind(1)-1 iind(end)+1:bb];
    
                % interpolated fields
                Xint(i,ind1) = nan;
                Yint(i,ind1) = nan;
                MDT_prime_om(i,ind1) = nan;
                U_prime_om(i,ind1) = nan;
                V_prime_om(i,ind1) = nan;
                Epsi(i,ind1) = nan;
%                 z1(i,ind1) = nan;
    
                % mean fields
                mdot_nii_mean(i,ind1) = nan;
                ur_nii_mean(i,ind1) = nan;
                vr_nii_mean(i,ind1) = nan;
            else
                Xint(i,:) = nan;
                Yint(i,:) = nan;
                MDT_prime_om(i,:) = nan;
                U_prime_om(i,:) = nan;
                V_prime_om(i,:) = nan;
                Epsi(i,:) = nan;
%                 z1(i,ind1) = nan;
    
                % mean fields
                mdot_nii_mean(i,:) = nan;
                ur_nii_mean(i,:) = nan;
                vr_nii_mean(i,:) = nan;    
            end
        end

        % Add back the mean
        mean_psi_om = MDT_prime_om + mdot_nii_mean;
        mean_ur_om = U_prime_om + ur_nii_mean;
        mean_vr_om = V_prime_om + vr_nii_mean;

        % Residuals (erros)
        res_om_ur = mean_ur_om - uui;
        res_om_vr = mean_vr_om - vvi;

        angle = atan2(res_om_vr,res_om_ur);

        angle_mean(m,1) = nanmean(angle(:));
        angle_var(m,1) = nanvar(angle(:));

        % 3D fields of objectively mapped \Psi, U, V
        om_mean_psi(m,:,:) = mean_psi_om;
        om_mean_u(m,:,:) = mean_ur_om;
        om_mean_v(m,:,:) = mean_vr_om;

        % 3D fields of objectively mapped anomalies
        om_prime_psi(m,:,:) = MDT_prime_om;
        om_prime_u(m,:,:) = U_prime_om;
        om_prime_v(m,:,:) = V_prime_om;

        % 3D fields of mapping Error
        om_Error(m,:,:) = Epsi;
        XXint(m,:,:) = Xint;
        YYint(m,:,:) = Yint;

        % 3D fields of residual
        res_Ur_om(m,:,:) = res_om_ur;
        res_Vr_om(m,:,:) = res_om_vr;
        angle_res(m,:,:) = angle;

        std_dv_realdata = sqrt(nanvar(uui(:),0,1) + nanvar(vvi(:),0,1));

        % RMS normalized by std deviation of original data
        RMS_om(m,1) = sqrt(nansum((res_om_ur(:).^2 + res_om_vr(:).^2) ...
            /length(find(isnan(res_om_ur(:))==0))))*(std_dv_realdata)^-1;

        % Computes normal vector
        psi_y = -mean_ur_om.*f_sm3./grav_sm3;
        psi_x = mean_vr_om.*f_sm3./grav_sm3;
        norm_vec(m,:,:) = complex(psi_x,psi_y)./abs(complex(psi_x,psi_y));

        %%%% Plots
        if plotm == 1;
            
            % Loads Bathymetry
            load('Bath_DrakePassage_ETOPO2','x_et2','y_et2','lon_et2',...
                'lat_et2','z_et2');

            
            if mod(m,8) == 1;
                fac = 150;    
                scale = 0.15;
                arch = [num2str(L(1)/1e3,'%2.0f') 'km_' ...
                    num2str(simple_mean_3D.z(m),'%3.0f') ...
                    'm_E' num2str(err*100,'%3.0f') '_MN12'];
                folder = '/Users/gvillanuma/Desktop/Thesis/Figures/sADCP/Objective mapping/Maximenko and Niiler/';

                figure('color','w');
                set(gcf,'PaperUnits','centimeters')
                xSize = 27; ySize = 27;
                xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
                set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
                set(gcf,'Position',[200 800 xSize*40 ySize*40]);


                dd = suptitle([' $ Lx = Ly = $ ' ...
                    num2str(L(1)/1000) ' km, $ E = $ ' ...
                    num2str(err,'%1.2f') ', MDT = MN12, z = ' ...
                    num2str(simple_mean_3D.z(m),'%4.0f') ' m '],18);
                set(dd,'interpreter','latex')

                subplot(221)
                [c1 h1] = contour(x_et2/1000,y_et2/1000,z_et2,...
                    [0 0],'k'); hold on;
%                 plot(area(:,1),area(:,2),'r','linewidth',2)
                text(-363,-140,[num2str(scale,'%1.2f') ' m/s '],...
                    'fontsize',14,'horizontalAlignment','center',...
                    'verticalAlignment','middle','color','k',...
                    'fontweight','bold')
                quiver(-400,-100,scale*fac,0*fac,0,'k',...
                    'linewidth',1,'maxheadsize',10)
                quiver(-400,-100,0*fac,scale*fac,0,'k',...
                    'linewidth',1,'maxheadsize',10)
                p1 = quiver(Xint/1000,Yint/1000,Ur_anom_upd.*fac,...
                    Vr_anom_upd.*fac,0,'k','linewidth',1);
                col0 = colorbar;
                axis equal; 
                set(gca,'fontsize',14,'tickdir','out','ticklen',...
                    [0.02 0.02]);
                ylabel(' km ','fontsize',14); %xlabel(' km ','fontsize',14); 
                axis([-500 300 -900 100]);
                set(gca,'position',pos1); colorbar(col0,'delete');
                ll1 = legend([p1],' $\vec{\textbf{u}}_{upd}^\prime$ ',...
                    'location','West');
                set(ll1,'fontsize',16,'interpreter','latex',...
                    'fontweight','bold');

                subplot(222)
                [c1 h1] = contour(x_et2/1000,y_et2/1000,z_et2,[0 0],'k'); 
                hold on;
                [p2 pp2] = contour(X/1000,Y/1000,MDT_prime_om,10,...
                    'linewidth',2);
                colormap(flipud(othercolor('RdBu10')));
%                 plot(area(:,1),area(:,2),'r','linewidth',2)
                text(-363,-140,[num2str(scale,'%1.2f') ' m/s '],...
                    'fontsize',14,'horizontalAlignment','center',...
                    'verticalAlignment','middle','color','k',...
                    'fontweight','bold')
                quiver(-400,-100,scale*fac,0*fac,0,'k','linewidth',1,...
                    'maxheadsize',10)
                quiver(-400,-100,0*fac,scale*fac,0,'k','linewidth',1,...
                    'maxheadsize',10)
                p1 = quiver(Xint/1000,Yint/1000,U_prime_om.*fac,...
                    V_prime_om.*fac,0,'k','linewidth',1);
                col1 = colorbar; set(col1,'fontsize',14,...
                    'position',colb1); 
                ylabel(col1,' $\mathbf{\psi}^\prime_{om}$   m ',...
                    'fontsize',14,'interpreter','latex');
                axis equal; caxis([-0.07 0.07]);
                set(gca,'fontsize',14,'tickdir','out','ticklen',...
                    [0.02 0.02]);
%                 xlabel(' km ','fontsize',14); ylabel(' km ','fontsize',14);
                axis([-500 300 -900 100]);
                set(gca,'position',pos2);
                ll1 = legend([p1],' $\vec{\textbf{u}}^\prime_{om}$ ',...
                    'location','West');
                set(ll1,'fontsize',16,'interpreter','latex','fontweight','bold');
%                title([' $ Lx = Ly = $ ' num2str(L(1)/1000) ' km, $ E = $ ' ...
%                     num2str(err,'%1.2f') ', MDT = CNES13, z = ' num2str(mean_noverlap.z(m),'%4.0f') ...
%                     ' m '],'fontsize',14,'interpreter','latex');

                subplot(223)
                [c1 h1] = contour(x_et2/1000,y_et2/1000,z_et2,...
                    [0 0],'k'); hold on;
                [p2 pp2] = contour(X/1000,Y/1000,mean_psi_om - 0.6,...
                    [-1.70:0.05:-0.35],'linewidth',2);
                colormap(flipud(othercolor('RdBu10')));
%                plot(area(:,1),area(:,2),'r','linewidth',2)
                text(-363,-140,[num2str(scale,'%1.2f') ' m/s '],...
                    'fontsize',14,'horizontalAlignment','center',...
                    'verticalAlignment','middle','color','k',...
                    'fontweight','bold')
                quiver(-400,-100,scale*fac,0*fac,0,'k',...
                    'linewidth',1,'maxheadsize',10)
                quiver(-400,-100,0*fac,scale*fac,0,'k',...
                    'linewidth',1,'maxheadsize',10)
                p1 = quiver(Xint/1000,Yint/1000,mean_ur_om.*fac,...
                    mean_vr_om.*fac,0,'k','linewidth',1);
                col1 = colorbar; set(col1,'fontsize',14,'position',colb2);
                ylabel(col1,' $\mathbf{\psi}_{om}$  m ',...
                    'fontsize',14,'interpreter','latex');
                axis equal; caxis([-1.7 -0.5])
                set(gca,'fontsize',14,'tickdir','out',...
                    'ticklen',[0.02 0.02]);
                xlabel(' km ','fontsize',14); ylabel(' km ','fontsize',14);
                axis([-500 300 -900 100]);
                set(gca,'position',pos3);
                ll1 = legend([p1],' $\vec{\textbf{u}}_{om}$ ',...
                    'location','West');
                set(ll1,'fontsize',16,'interpreter','latex',...
                    'fontweight','bold');

%                pos = get(ax1,'Position');
%                title([' $ Lx = Ly = $ ' num2str(L(1)/1000) ' km, $ E = $ ' ...
%                    num2str(err,'%1.2f') ', MDT = CNES13, z = ' ...
%                    num2str(mean_noverlap.z(m),'%4.0f') ...
%                    ' m '],'fontsize',14,'interpreter','latex');


                subplot(224)
                [c1 h1] = contour(x_et2/1000,y_et2/1000,z_et2,[0 0],...
                    'k'); hold on;
                [p1 pp2] = contour(X/1000,Y/1000,Epsi,...
                    [0:0.1:1],'linewidth',2);
                colormap(flipud(othercolor('RdBu10')));
%                plot(area(:,1),area(:,2),'r','linewidth',2)
                text(-363,-140,[num2str(scale,'%1.2f') ' m/s '],...
                    'fontsize',14,'horizontalAlignment','center',...
                'verticalAlignment','middle','color','k',...
                'fontweight','bold');
                quiver(-400,-100,scale*fac,0*fac,0,'k',...
                    'linewidth',1,'maxheadsize',10)
                quiver(-400,-100,0*fac,scale*fac,0,'k',...
                    'linewidth',1,'maxheadsize',10)
                p2 = quiver(Xint/1000,Yint/1000,res_om_ur.*fac,...
                    res_om_vr.*fac,0,'k','linewidth',1);
                col2 = colorbar; set(col2,'fontsize',14,'position',colb3); 
                axis equal; caxis([0 1]);
                set(gca,'fontsize',14,'tickdir',...
                    'out','ticklen',[0.02 0.02]);
                xlabel(' km ','fontsize',14); %ylabel(' km ','fontsize',14);
                axis([-500 300 -900 100]);
                set(gca,'position',pos4);
                ll1 = legend([p2],' $\vec{\textbf{u}}_{res}$ ',...
                    'location','West');
                set(ll1,'fontsize',16,'interpreter','latex',...
                    'fontweight','bold');
                text(-470,-10,['RMS = ' num2str(squeeze(RMS_om(m)),...
                    '%2.4f') ],'interpreter','latex','fontsize',13);
%                 set(ax2,'Position',...
%                     [pos(1)+pos(3)+0.15 pos(2) pos(3) pos(4)],...
%                     'YTickLabel',[]);
%                 title([' $ Lx = Ly = $ ' num2str(L(1)/1000) ...
%                     ' km, $ E = $ ' num2str(err,'%1.2f') ...
%                     ', MDT = CNES13, z = ' num2str(mean_noverlap.z(m), ...
%                     '%4.0f') ' m '],'fontsize',14,'interpreter','latex');
% 
%                 print('-depsc2','-painters',[folder arch ...
%                     num2str(overlapping*100)]);
            end
        end
    end

toc


sm3_om.uvmean = complex(om_mean_u,om_mean_v);
sm3_om.psi = om_mean_psi - 0.6;
sm3_om.err = om_Error;
% sm3_om.u_prime = ueddy;
% sm3_om.v_prime = veddy;
% sm3_om.eke = EKE;
% sm3_om.mke = MKE;
sm3_om.normal = norm_vec;
sm3_om.x = X;
sm3_om.y = Y;
sm3_om.lon = simple_mean_3D.lon;
sm3_om.lat = simple_mean_3D.lat;
sm3_om.dof = simple_mean_3D.num;
sm3_om.z = simple_mean_3D.z;
sm3_om.lat0 = lat0; 
sm3_om.lon0 = lon0;
sm3_om.theta0 = theta0;


