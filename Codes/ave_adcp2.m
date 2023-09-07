function sm = ave_adcp2(nd, datatype, itype, varargin);
% sm = ave_adcp(nd, datatype, itype, varargin)
% computes averages of the gridded data. 
%
% nd determines the number of dimensions for the data, valid options are 
%         2 (depth-averaged) or 3 (not). 

% dataype may be either 'raw' or 'no upd' for the detided ADCP data 
%         or the detided ADCP data with the AVISO sea level anomalies removed, respectively. 
%
% itype is 'nb150' or 'os38nb'
%
% varargin is gadcp, maxz, where gadcp is the structure of gridded data
%         to use (see grid_adcp.m) and maxz is the maximum depth to average
%         over. 3rd is overlapping 0 < overlapping < 1
%
% maxz defaults to 300 m for nb150, 1030 for os38nb
%
% or varargin can have gadcp, otherwise will load

% mean and std.err. will be computed where there are data from 3 or more 
% cruises, principal axis information where there are 15 or more cruises
%
% modified YLF April 2008:
%        check that at least 2 values (EACH DEPTH) used to compute standard error when using 3D
%
% se changed code to match YLF's, 02-21-12, so that matlab's statistic's toolbox (nanstd) 
%   is not needed when at sea.  Using mstdgap (UHDAS distribution) instead)
%

% se added YLF code so number_data_pts (ncd) is now also being saved 02-21-12
%
% as per YLF,  changed '.use and .vse are the standard deviation (m/s) in the currents'
%     to read "standard error" which is correct. 02-21-12
%
% added line ' uv = uv+i*uv; % 02-21-12 as per YLF added this line' as per YLF
%
% changed code so nanstd,nanmean (statistics tool box) are not used.  replaced with mstdgap (UHDAS) 02-21-12

% changed lines as per YLF 02-21-12
%      %uv(iil,I(n)) = NaN;      u_ste(iil,I(n)) = NaN; v_ste(iil,I(n)) = NaN;
%      uv(iil,I(n)) = NaN+i*NaN; u_ste(iil,I(n)) = NaN; v_ste(iil,I(n)) = NaN;
%              
% changed 'for d = 1:length(gadcp.z)' to ' for d = 1:length(iz)'  as per YLF 02-21-12
%    so only depths shallower than maxz are used.  This is kinda a bug in our previous
%    code,  because all good depths were used,   not just good depths shallower than
%    maxz.
%
% se changed code to match YLF's, 06-05-12, to use YLF's nanmstd.m for statistics
% see "09-2012" below for 4 changes made, so as not to transpose.
%
% modified by MOGV (08/19/2019): saves lat0,lon0 and theta0 from the gadcp.

%================================================================================

if length(varargin) > 0
   gadcp = varargin{1};
   maxz = varargin{2};
   overlapping = varargin{3};
elseif instring(datatype,'raw')
   load([cd '/drake_mean/data/all_drake/lmg_' itype], 'gadcp')
   if instring(itype, 'nb150'); maxz = 300; else; maxz = 1030; end
elseif instring(datatype, 'no upd')
   load([cd '/drake_mean/data/all_drake/lmg_' itype], 'gadcp_improved')
   gadcp = gadcp_improved; clear gadcp_improved;
   if instring(itype, 'nb150'); maxz = 300; else; maxz = 1030; end
else
   error('invalid data type specified and no gadcp provided');
end
clear gadcp_*;

sm = struct('x', gadcp.x, 'y', gadcp.y, 'lon', gadcp.lon, 'lat', gadcp.lat, 'dchoice', gadcp.dchoice, 'num', gadcp.cruises);

sm.doc = char(  '.x and .y are the grid box locations in Drake Passage coordinates',...
                '.lon and .lat are the grid box locations in geographic coordinates',...
                '.dchoice is the choice of data used, raw = just ADCP, no upd = ADCP - SLA',...
                '.num is the number of time each grid box has been crossed',...
                '.u and .v are the mean velocities (m/s) projected into Drake Passage coordinates',...
                '.use and .vse are the standard errors (m/s) in the currents',...
                '.theta, .major and .minor are the angle of rotation(rad), major axis(m/s) and minor axis(m/s)',...
                ' of the  standard deviation ellipses',...
                '.z are the depths(m) of the velocities uesd');


iz = find(gadcp.z<=maxz);
I = find(gadcp.cruises>6);
uv = ones(length(iz),size(gadcp.u,1),size(gadcp.u,2))*NaN;
u_ste = uv; v_ste = uv;
ncd = uv ;  % added as per YLF to now also save number_data_pts (ncd)
if nd == 2; 
   theta = ones(size(gadcp.u,1),size(gadcp.u,2))*NaN;
   major = theta; minor = theta;
   M = theta; N = theta;
elseif nd == 3
   theta = uv;
   major = theta; minor = theta;
   M = theta; N = theta;
end

uv = uv+i*uv; % 02-21-12 as per YLF added this line

for n = 1:length(I)

   if gadcp.cruises(I(n))>=3
      le = sum(~isnan(gadcp.u{I(n)}(iz,:))')';

      % se changed code below as per YLF 02-21-12, to not use nanmean and nanstd
      %uv(:,I(n)) = complex(nanmean(gadcp.u{I(n)}(iz,:)'), nanmean(gadcp.v{I(n)}(iz,:)'));
      %u_ste(:,I(n)) = nanstd(gadcp.u{I(n)}(iz,:)')'./sqrt(le);
      %v_ste(:,I(n)) = nanstd(gadcp.v{I(n)}(iz,:)')'./sqrt(le);
      %[m,s,junk] = mstdgap(gadcp.u{I(n)}(iz,:)');
      [m,s,junk] = nanmstd(gadcp.u{I(n)}(iz,:),2);
      uv(:,I(n)) = m;                                % changed from m' to m 09-2012
      u_ste(:,I(n)) = s./sqrt(le);                   % changed from s' to s 09-2012
      %[m,s,junk] = mstdgap(gadcp.v{I(n)}(iz,:)');
      [m,s,junk] = nanmstd(gadcp.v{I(n)}(iz,:),2);
      uv(:,I(n)) = uv(:,I(n)) + i*m;                 % changed from m' to m 09-2012
      v_ste(:,I(n)) = s./sqrt(le);                  % changed from s' to s 09-2012


      % next two lines NOT in YDL code.  YLF includes to ensure 2+ values to make a ste.
      %    these lines make the final .mat files slightly different.
      iil = find(le<3);

      % as per YLF 02-21-12 changed lines
      %uv(iil,I(n)) = NaN;      u_ste(iil,I(n)) = NaN; v_ste(iil,I(n)) = NaN;
      uv(iil,I(n)) = NaN+i*NaN; u_ste(iil,I(n)) = NaN; v_ste(iil,I(n)) = NaN;

      ncd(:,I(n)) = le;  % added as per YLF 02-21-12

      if gadcp.cruises(I(n))>14
         if nd == 2
            %[theta(I(n)), major(I(n)), minor(I(n))] = princax(nanmean(complex(gadcp.u{I(n)}, gadcp.v{I(n)}))); % se added from YDL
            % YLF code, below, uses a max depth,  so all depths are not used as in YDL
             %[theta(I(n)), major(I(n)), minor(I(n))] = princax(nanmean(complex(gadcp.u{I(n)}(iz,:), gadcp.v{I(n)}(iz,:)),1)); 
            [uvd, s, junk] = nanmstd(complex(gadcp.u{I(n)}(iz,:), gadcp.v{I(n)}(iz,:)));
            [theta(I(n)), major(I(n)), minor(I(n)), w, M(I(n)), N(I(n))] ...
                = princax(uvd);
         elseif nd ==3
            % se changed as per YLF,  to use only depths shallower than max depth 02-21-12
            %for d = 1:length(gadcp.z)  % if 3D,  do it for each depth,
            for d = 1:length(iz)  % if 3D,  do it for each depth,
               uvd = complex(gadcp.u{I(n)}(d,:), gadcp.v{I(n)}(d,:));
               uvd(isnan(uvd)) = [];
               if length(uvd)>1
%                   [theta(d,I(n)), major(d,I(n)), minor(d,I(n)) w ...
%                       M(d,I(n)), N(d,I(n))] = ...
%                       princax(uvd);
                  [theta(d,I(n)), major(d,I(n)), minor(d,I(n)) w] = ...
                      princax(uvd);
               end
            end
         end
      end
	
   end

end

if nd == 3 ;
   sm.u = real(uv(iz,:,:));
   sm.v = imag(uv(iz,:,:));
   sm.use= u_ste(iz,:,:);
   sm.vse = v_ste(iz,:,:);
   sm.ncd = ncd;  % 02-21-12  added as per YLF
elseif nd == 2
   %sm.u = squeeze(nanmean(real(uv(iz,:,:))));
   %sm.v = squeeze(nanmean(imag(uv(iz,:,:))));
   %sm.use = squeeze(sqrt(nanmean(u_ste(iz,:,:).^2)));
   %sm.vse = squeeze(sqrt(nanmean(v_ste(iz,:,:).^2)));
   [sm.u, s, junk] = nanmstd(real(uv(iz,:,:)),1,'squeeze');
   [sm.v, s, junk] = nanmstd(imag(uv(iz,:,:)),1,'squeeze');
   [m, s, junk] = nanmstd(u_ste(iz,:,:).^2,1,'squeeze');
   sm.use = sqrt(m);
   [m, s, junk] = nanmstd(v_ste(iz,:,:).^2,1,'squeeze');
   sm.vse = sqrt(m);
end

sm.theta = theta;
sm.major = major;
sm.minor = minor;
% sm.nc = gadcp.nc;
sm.z = gadcp.z(iz);
% sm.M = real(w);
% sm.N = imag(w);
sm.lat0 = gadcp.lat0;
sm.lon0 = gadcp.lon0;
sm.theta0 = gadcp.theta0;
