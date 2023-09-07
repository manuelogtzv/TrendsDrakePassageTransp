function [gadcp all_lmg] = grid_adcp_month(gsize,data,itype,varargin);
% GADCP = GRID_ADCP(GSIZE,DATA,itype,theta) grids the ADCP velocities in grid boxes
% of dimension GSIZE (km), the ADCP velocities are detided, rotated into 
% passage coordinates (given in DPshort) and averaged by cruise. 
%
% Choice of DATA determines whether the sea level anomalies are removed. 
%   DATA = 'raw' (upd not removed, bt tide is removed)
%        = 'improved' (upd  removed),   NOTE: same as YLF's old 'no_upd'
%
% itype is 'nb150' or 'os38nb'
%
% theta is optional angle to rotate by (degrees), if unspecified it will 
% use YDL's from DPshort
%
% created by YDL, june 2006
% modified by YLF, april 2008
% se changed code to match YLF's, 06-05-12, to use YLF's nanmstd.m for statistics
%
% Modified by MOGV (08/19/2019): Allows to specify the lat0, lon0, theta
% and the climatological months.


if length(varargin)==0
   load ([cd '/drake_mean/data/all_drake/', 'DPshort'] , 'theta');
   theta0 = theta; clear theta;
   overlapping = 0;
elseif length(varargin) == 1;
   theta0 = varargin{1};
   lat0 = -55;
   lon0 = -65;
elseif length(varargin) == 2;
   theta0 = varargin{1};
   overlapping = varargin{2};
   lat0 = -55;
   lon0 = -65;
elseif length(varargin) == 4;
   theta0 = varargin{1};
   overlapping = varargin{2};
   lat0 = varargin{3};
   lon0 = varargin{4};
elseif length(varargin) == 5;
   theta0 = varargin{1};
   overlapping = varargin{2};
   lat0 = varargin{3};
   lon0 = varargin{4};
   months = varargin{5};
   
end


load([cd '/drake_mean/data/all_drake/lmg_' itype]);

if strcmp(itype,'nb150') == 1;
    all_lmg.uv = [repmat(all_lmg.uv(1,:),3,1); all_lmg.uv];
    all_lmg.z = [0;10;18;all_lmg.z];
elseif strcmp(itype,'os38nb') == 1;
    all_lmg.uv = [repmat(all_lmg.uv(1,:),2,1); all_lmg.uv];
    all_lmg.z = [0;22;all_lmg.z];
end

lon_d = all_lmg.lon;
lat_d = all_lmg.lat;

% Extracts data that is within specific months
if length(varargin) == 5;% Finds data between specific dates
    
    % Get indices for each transect
    cross_lmg = interval(all_lmg.t,3,'n');
    
    % XBT transects
    for i = 1:length(cross_lmg);
        mean_time_lmg(i,1) = mean(all_lmg.t(cross_lmg{i}));
        date_vec(i,:) = datevec(mean_time_lmg(i,1));
    end

    no_months = length(months);% Number of months specified
    
    ind_months = [];
    
    for p = 1:no_months;% Loops for each month
       % Indexes for mean transect time
       ind_date = find(date_vec(:,2) ==  months(p));
       
       ind_month = [];

       % Concatenates indexes for each month
       for y = 1:length(ind_date);
          ind_month = [ind_months cross_lmg{ind_date}];   
       end
       
       % Concatenates all indexes for all specified months
       ind_months = [ind_months ind_month];
    end
    
    all_lmg.uv = all_lmg.uv(:,ind_months);
    all_lmg.uvupd = all_lmg.uvupd(ind_months);
    all_lmg.t = all_lmg.t(ind_months);
    all_lmg.lon = all_lmg.lon(ind_months);
    all_lmg.lat = all_lmg.lat(ind_months);   
    
    tide.u = tide.u(ind_months);
    tide.v = tide.v(ind_months);
    
    [all_lmg.t I] = sort(all_lmg.t,'ascend');
    all_lmg.uvupd = all_lmg.uvupd(I);
    all_lmg.uv = all_lmg.uv(:,I);
    all_lmg.lon = all_lmg.lon(I);
    all_lmg.lat = all_lmg.lat(I);
    
    clear I
    
%     ind_dates = find(temp_time >= date_i & temp_time <= date_f);
%     
%     temp_new2 = temp_new2(:,ind_dates);
%     temp_new = temp_new(:,ind_dates);
%     temp_time = temp_time(ind_dates);
%     temp_lon = temp_lon(ind_dates);
%     temp_lat = temp_lat(ind_dates);
end

nz = length(all_lmg.z);
dd = all_lmg.t; 

disp('before removing barotropic tide');
all_lmg.uv = all_lmg.uv - repmat(complex(tide.u', tide.v'),nz,1);
if instring(data, 'improved')
   disp('subtracting geostrophic velocity anomalies');
   ur = real(all_lmg.uv - repmat(shiftdim(all_lmg.uvupd,1),nz,1));%changed by MGV
   vr = imag(all_lmg.uv - repmat(shiftdim(all_lmg.uvupd,1),nz,1));
   
%    all_lmg.uv = all_lmg.uv - repmat(shiftdim(all_lmg.uvupd,1),nz,1);  % DO NOT use transpose.  MUST use shiftdim
end



% Improved ADCP velocities (no tides no geostrophic velocity anomalies)
ur = real(complex(ur,vr)*exp(-sqrt(-1)*theta0/180*pi));
vr = imag(complex(ur,vr)*exp(-sqrt(-1)*theta0/180*pi));

% ADCP data with barotropic tide removed.
all_lmg.uv = all_lmg.uv*exp(-sqrt(-1)*theta0/180*pi);

% ur = real(all_lmg.uv*exp(-sqrt(-1)*theta/180*pi));
% vr = imag(all_lmg.uv*exp(-sqrt(-1)*theta/180*pi));


[gpos, xr, yr] = makegrid2...
    (all_lmg.lat, all_lmg.lon-360,gsize,0,theta0,lat0,lon0); % compute grid positions

% Make the along-across DP grid relative to the grid obtained from all the
% data. This ensures that the grid stays the same regardless of the time
% period used. Using a smaller subset of transects affects the size of the
% along/across DP grid.
[gpos3, xr3, yr3] = makegrid2...
    (lat_d,lon_d-360,gsize,0,theta0,lat0,lon0); % compute grid positions

% saves (x,y) of data
all_lmg.x = xr;
all_lmg.y = yr;
all_lmg.lat0 = lat0;
all_lmg.lon0 = lon0;
all_lmg.theta0 = theta0;

% initialize matrices
U = cell(length(gpos3.y), length(gpos3.x)); V = U; DD = U;
llon = U;
llat = U;
gadcp = struct('x', gpos3.x, 'y', gpos3.y, 'lon', gpos3.lon,...
    'lat', gpos3.lat, 'dchoice', data, 'cruises', nan(length(gpos3.y),...
    length(gpos3.x)));

ur_all = [];
vr_all = [];
xr_all = [];
yr_all = [];
time_all = [];

for nx = 1:length(gpos3.x)
   for ny = 1:length(gpos3.y);
      I = find(abs(xr-gpos3.x(nx)) < gsize*(0.5 + overlapping) ...
          & abs(yr-gpos3.y(ny)) < gsize*(0.5 + overlapping));
      if ~isempty(I)
         % counting number of crossings per grid box 
         cr = interval(dd(I),3,'n');
         ubin  = ones(size(ur,1),length(cr))*NaN;
         vbin = ubin;
         lonbin  = ones(length(cr),1)*NaN;
         latbin  = lonbin;
         
         jd = ubin(1,:);
         % averaging observations in a single crossing
         for cn =  1:length(cr);
            if length(cr{cn}) > 1
                
               %ubin(:,cn) = nanmean(ur(:,I(cr{cn}))')';
               %vbin(:,cn) = nanmean(vr(:,I(cr{cn}))')';
               [lonbin(cn,1),s,nn] = nanmstd(xr(I(cr{cn})));
               [latbin(cn,1),s,nn] = nanmstd(yr(I(cr{cn})));
               [ubin(:,cn),s,nn] = nanmstd(ur(:,I(cr{cn})),2,'squeeze');
               [vbin(:,cn),s,nn] = nanmstd(vr(:,I(cr{cn})),2,'squeeze');
               
%                [ddbin(:,cn),s,nn] = nanmstd(dd(I(cr{cn})));
%             else
%                ubin(:,cn) = ur(:,I(cr{cn}));
%                vbin(:,cn) = vr(:,I(cr{cn}));
%                lonbin(:,cn) = xr(I(cr{cn}));
%                latbin(:,cn) = yr(I(cr{cn}));
            end
            %jd(cn) = nanmean(dd(I(cr{cn})));
            %[jd(cn),s,nn] = mstdgap(dd(I(cr{cn})));
            [jd(cn),s,nn] = nanmstd(dd(I(cr{cn})));
         end 
            
         gadcp.cruises(ny,nx) = length(cr);
         DD{ny,nx} = jd;
         U{ny,nx} = ubin;
         V{ny,nx} = vbin;
         llon{ny,nx} = lonbin;
         llat{ny,nx} = latbin;
         
%          % Concatenates
%          ur_all = [ur_all U{ny,nx}];
%          vr_all = [vr_all V{ny,nx}];
%          xr_all = [xr_all; llon{ny,nx}(:)];
%          yr_all = [yr_all; llat{ny,nx}(:)];
%          time_all = [time_all; jd(:)];
         
      end
   end
end

gadcp.u = U;
gadcp.v = V;
gadcp.timeg = DD;
gadcp.nc = length(interval(all_lmg.t,2));
gadcp.z = all_lmg.z;
gadcp.lat0 = lat0;
gadcp.lon0 = lon0;
gadcp.theta0 = theta0;
% gadcp.xrg = llon;
% gadcp.yrg = llat;
% gadcp.ur_all = ur_all;
% gadcp.vr_all = vr_all;
% gadcp.xr_all = xr_all;
% gadcp.yr_all = yr_all;
% gadcp.time_all = time_all;

gadcp.doc = char('.x and .y are the grid box locations in Drake Passage coordinates (km)', ...
            '.lon and .lat are the grid box locations in geographic coordinates', ...
            '.dchoice is the choice of data used, raw = just ADCP, no upd = ADCP - SLA', ...
            '.cruises is the number of time each grid box has been crossed', ...
            '.u and .v are the velocities(m/s) projected into Drake Passage coordinates', ...
            '.jd are the decimal days from 1 jan 1999', ...
            '.nc are total number of cruises included');
