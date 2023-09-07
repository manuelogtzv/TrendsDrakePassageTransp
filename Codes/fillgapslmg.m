function [data] = fillgapslmg(data, indgap);

% [data] = fillgapslmg(data);
% 
% Function to fill data gaps that span the entire sampled water column in
% the LMG transects. Fills with linearly interpolated longitudes and
% latitudes when gaps are larger than 25 km.
%
% Input:
%
% data - must be a structure variable containing time, depth, longitude,
% latitude and uv
%
% Manuel O. Gutierrez Villanueva
%
% 2021/08/31
%
% 2022/03/28 - Includes the option for distance or latitude as flags.

% if strcmp(flag, 'distance')
%     dist = [sw_dist(data.lat', data.lon', 'km')];
%     indgap = find(abs(dist)>gap);
% elseif strcmp(flag, 'latitude')
%     dist = diff(data.lat);
%     indgap = find(abs(dist)>gap);
% else
%     error('Flag not recognized. Only distance or latitude are valid options');
% end


if ~isempty(indgap);
    
    xnew = [];
    ynew = [];
    tnew = [];
    uvnew = [];
            
    for nn = 1:length(indgap);
        
        if data.lat(1)>-60;% if starts from the north
             ygap = [-data.lat(indgap(nn)):0.01:...
                 -data.lat(indgap(nn)+1)]*-1;
             xgap = interp1(data.lat, data.lon, ygap, 'linear');
             tgap = interp1(data.lat, data.time, ygap, 'linear');
        else
             ygap = data.lat(indgap(nn)):0.01:...
                 data.lat(indgap(nn)+1);
             xgap = interp1(data.lat, data.lon, ygap, 'linear');
             tgap = interp1(data.lat, data.time, ygap, 'linear');
        end
                
        if nn == 1;
             xnew = [xnew; data.lon(1:indgap(nn)); xgap(2:end-1)'];
             tnew = [tnew; data.time(1:indgap(nn)); tgap(2:end-1)'];
             ynew = [ynew; data.lat(1:indgap(nn)); ygap(2:end-1)'];
             uvnew = [uvnew data.uv(:, 1:indgap(nn)) ...
                  nan(length(data.z), length(xgap(2:end-1))) + sqrt(-1)*...
                  nan(length(data.z), length(xgap(2:end-1)))];
        else
             xnew = [xnew; data.lon(indgap(nn-1)+1:indgap(nn)); xgap(2:end-1)'];
             tnew = [tnew; data.time(indgap(nn-1)+1:indgap(nn)); tgap(2:end-1)'];
             ynew = [ynew; data.lat(indgap(nn-1)+1:indgap(nn)); ygap(2:end-1)'];
             uvnew = [uvnew data.uv(:, indgap(nn-1)+1:indgap(nn)) ...
                  nan(length(data.z), length(xgap(2:end-1))) + sqrt(-1)*...
                  nan(length(data.z), length(xgap(2:end-1)))];
        end
        
    end
                
    xnew = [xnew; data.lon(indgap(nn)+1:end)];
    tnew = [tnew; data.time(indgap(nn)+1:end)];
    ynew = [ynew; data.lat(indgap(nn)+1:end)];
    uvnew = [uvnew data.uv(:, indgap(nn)+1:end)];
            
    data.lon = xnew;
    data.lat = ynew;
    data.time = tnew;
    data.uv = uvnew;
end

