function [trendfull confintfull sigtrend_ls trendfull_thse sigtrend_mk] = calctrends2D(data, time, alpha, season_opt);

% [trendfull confintfull sigtrend_ls trendfull_theinseil sigtrend_mk] = calctrends2D(data, time, alpha, season_opt);
% 
% Function to calculate trends and their signifance using least-squares fit
% (parametric) and Thein-seil estimate (non parametric). Significance for
% the least-squares fit is estimated using Fay and Mackinley (2013)
% approach, which uses the two-tailed t-student distribution at 1-alpha/2
% confidence. The trends estimated with the Thein-seil estimate are
% assessed using the modified Mann-Kendall test, which accounts for the
% autocorrelation function, therefore, uses the effective degrees of
% freedom.
%
% Inputs:
%
% + data: contains data in a matrix of size [N x M x O],
% where the M is the number of observations in y, N is the number of
% stations in x, and O is the number of observations in time  
% + time: time vector of size [1 x M].
% + alpha: confidence [0 < alpha < 1].
% + season_opt: option to remove seasonal cycle [annual + semiannual]
% 
% Outputs:
%
% + trendfull: rate of change estiamted from the least-squares fit. Size is
% [O x N].
% + confintfull: confidence interval for the trendfull. Size is equal to
% size(trendfull).
% + sigtrend_ls: 1 for trend being statistically different from 0.
% + trendfull_thse: trends estimated using the Thein-Seil estimator.
% Seasonal cycle (annual + semiannual harmonic) removed before computing
% the trend. Size equal to size(trendfull).
% + sigtrend_mk: 1 for non-parametric trend being statistically different
% from 0. Size equal to size(trendfull_theinseil).
%
%
% Manuel O. Gutierrez-Villanueva 2022/06/13
%
% 2022/07/18 - Option to remove seasonal cycle before performing
% least-squares fit

% checks if alpha is included 
if nargin < 3;
    alpha = 0.05; %standard alpha 
elseif nargin < 4;
    season_opt = 0;
end


% checks size of data and time
[mt, nt] = size(time);
[md, nd, od] = size(data);

if nt~=od
        error('check: data & time must have same dimensions')
else
        % Pre-allocates memory
        trendfull = nan(md, nd);
        confintfull = trendfull;
        sigtrend_ls = nan(size(trendfull));
        sigtrend_mk = sigtrend_ls;
        trendfull_thse = trendfull;

end

% Calculates trends and significances
for m = 1:size(trendfull, 1);
    for n = 1:size(trendfull, 2);
        yy = squeeze(data(m, n, :)); % ptemp
        
        if sum(~isnan(yy))/length(yy)>=0.5; %Minimum number of transects per box
%             xx = squeeze(geosbaroc.time); % time vector
            yy = yy(~isnan(yy)); % remove nans
%             xx = geosbaroc.time(~isnan(yy)); %remove nans
            xx1 = time(~isnan(yy));

            % Ordinary least-squares
            phianual = 2*pi*xx1;
            phisemi = 2*phianual;

            Xfull = [ones(length(xx1), 1) xx1(:) ...
                cos(phianual(:)) sin(phianual(:)) cos(phisemi(:)) ...
                sin(phisemi(:))];

            % Remove seasonal cycle before doing trend
            if season_opt == 1;
                [~, yy] = fitseasoncycle(yy(:), xx1(:));
            end

            
            bregfull = regress(yy(:), Xfull);
            trendfull(m, n) = bregfull(2);

            % Confidence intervals (Fay and McKinley 2013) 
            rmsfull = sqrt(sum((mean(yy) - Xfull*bregfull).^2)/...
                [length(xx1) - 6]);
            tstufull = tinv(alpha/2, length(xx1) - 6 );
            meanstdfull = sqrt(1./sum((xx1(:) -...
                mean(xx1(:))).^2));
            confintfull(m, n) = tstufull.*rmsfull*meanstdfull;

            if trendfull(m, n) - abs(confintfull(m, n)) > 0
                sigtrend_ls(m, n) = 1;
            else
                sigtrend_ls(m, n) = 0;
            end

            % Performs Mann - Kendall test for time series with the
            % seasonal cycle removed
            yy_noseas = yy - Xfull(:, 3:end)*bregfull(3:end);
            [sigtrend_mk(m, n), p_value] = ...
                Mann_Kendall_Modified(yy_noseas(:), alpha);

            % Trends per Thein-Seil test
            trendfull_thse(m, n) = Theil_Sen_Regress(xx1(:), ...
                yy_noseas(:));

        end
    end
end



