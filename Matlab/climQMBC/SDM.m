function SDM_series = SDM(obs,mod,SDM_var,frq,pp_threshold,pp_factor,day_win)
%% SDM_series
%  This function performs bias correction of modeled series based on
%  observed data by the Scaled Distribution Mapping (SDM) method, as
%  described by Switanek et al. (2017). Correction is performed to daily,
%  monthly or annual precipitation or temperature data in a single location.
%  The historical period of the modeled series is bias corrected as a
%  scenario where the complete series is replaced by the bias corrected
%  series. On the other hand, for each year after the historical period the
%  method considers a projected period as a window whose length is equal
%  to the number of years of the historical period and ends in the analyzed
%  year. From the bias corrected series, only the last year is saved and
%  all the others are discarded. This approach aims to build a continuum
%  for the bias corrected series.
%
% Input:
%   obs = A column vector of daily, monthly or annual observed data, without
%         considering leap days. If daily frequency is specified, the length
%         of the column vector should by a multiple of 365 and for monthly
%         frequency, it should be a multiple of 12. [ndata_obs, 1]
%
%   mod = A column vector of daily, monthly or annual modeled or GCM data,
%         without considering leap days. If daily frequency is specified,
%         the length of the column vector should by a multiple of 365 and
%         for monthly frequency, it should be a multiple of 12.
%         [ndata_mod, 1]
%
%   SDM_var = A flag that identifies if data are temperature or precipitation.
%             Temperature:   SDM_var = 0
%             Precipitation: SDM_var = 1
%
%   NOTE: This routine considers that obs and mod series start in the same
%         day/month/year and are continuous until the end day/month/year.
%
% Optional inputs:
%   frq = A string specifying if the input frequency is daily, monthly or
%         annual.
%          Daily:     frq = 'D'
%          Monthly:   frq = 'M'
%          Annual:    frq = 'A' (default)
% 
%   pp_threshold = A float indicating the threshold to consider no-rain
%                  values. Default is 1.
%
%   pp_factor = A float which multiplied to pp_threshold indicates the
%               maximum value of the random values that replace no-rain
%               values. Default is 1/100.
%
%   day_win = (only for frq='D') An integer indicating how many days to
%             consider backwards and forward to get the statistics of each
%             calendar day.The length of the window will be (2*win_day-1).
%             For example, day_win=15 -> window of 29. Default: win = 1
%
% Output:
%   SDM_series = A column vector of data bias corrected with the QM method.
%                [ndata_mod, 1]
%
% References:
%   Switanek, B. M., Troch, P., Castro, L. C., Leuprecht, A., Chang, H. I.,
%   Mukherjee, R., and Demaria, M. C. E. (2017) Scaled distribution
%   mapping: A bias correction method that preserves raw climate model
%   projected changes. Hydrology &amp; Earth System Sciences, 21, 2649-2666,
%   https://doi.org/10.5194/hess-21-2649-2017.
%
% Written by Sebastian Aedo Quililongo (1*)
%            Cristian Chadwick         (2)
%            Fernando Gonzalez-Leiva   (3)
%            Jorge Gironas             (3, 4)
%           
%  (1) Stockholm Environment Institute, Latin America Centre, Bogota,
%      Colombia
%  (2) Faculty of Engineering and Sciences, Universidad Adolfo Ibanez,
%      Santiago, Chile
%  (3) Department of Hydraulics and Environmental Engineering, Pontificia 
%      Universidad Catolica de Chile, Santiago, Chile
%  (4) Centro de Cambio Global UC, Pontificia Universidad Catolica de Chile,
%      Santiago, Chile
%
% *Maintainer contact: sebastian.aedo.q@gmail.com
% Revision: 1, updated Apr 2024
%


%%

% Check for optional arguments
if ~exist('pp_threshold','var')
    pp_threshold = 1;
end

if ~exist('pp_factor','var')
    pp_factor = 1/100;
end

if ~exist('day_win','var')
    day_win = 1;
end

if ~exist('frq','var')
    frq = 'A';
end

% Define the no-rain value and the threshold to truncate the tail of the
% probability ditribution functions
lower_lim = pp_threshold;
CDF_th = 10^-3;

% 0) Format series to matrix with rows as sub-periods and columns as years
%    and, if needed, replace no-rain values with random small values 
[y_obs,obs_series] = formatQM(obs, SDM_var, frq,pp_threshold,pp_factor);
[y_mod,mod_series] = formatQM(mod, SDM_var, frq,pp_threshold,pp_factor);
modh_series = mod_series(:,1:y_obs);

% If frequency is daily, consider a moving window for each day to fit
% the statistics
if frq=='D'
    % Get a 3D array with a moving window centered on each day
    obs_series = day_centered_moving_window(obs_series, day_win);
    mod_series = day_centered_moving_window(mod_series, day_win);
    
    % Get a 4D array with a backward moving window for each projected period
    win_series = projected_backward_moving_window(mod_series, y_obs, frq);
    
    % Reshape to 2D in order to have all the days of the window in each row
    obs_series = reshape(permute(obs_series, [1,3,2]),[365,(2*day_win-1)*y_obs]);
    modh_series = reshape(permute(mod_series(:,:,1:y_obs), [1,3,2]),[365,(2*day_win-1)*y_obs]);
    
    % Reshape to 3D in order to have all the days of the window in each row
    win_series = reshape(permute(win_series,[1,4,2,3]), [365, (2*day_win-1)*y_obs, y_mod-y_obs]);
    
    % Add the historical period to the moving window
    win_series_ = cat(3,modh_series, win_series);
else
    % For non-daily frequency, make sure that day_win=1
    day_win = 1;
    
    % Get a 3D array with a backward moving window for each projected period
    win_series = projected_backward_moving_window(mod_series, y_obs, frq);
    
    % Add the historical period to the moving window
    win_series_ = cat(3,modh_series, permute(win_series, [1,3,2]));
end

SDM  = zeros(size(mod_series,1),y_mod-y_obs);
SDM_h  = zeros(size(obs_series,1),y_obs);
for m = 1:size(mod_series,1)
    % 2) Historical period:
    
    % a) [Switanek et al. (2017), step 1)]
    %     For temperature, get the detrended modeled and observed series in
    %     the historical period.
    %     For precipitation, get the rainday values and its frequency for
    %     the modeled and observed series in the historical period.
    if SDM_var == 0
        D_obs = detrend(obs_series(m,:));
        D_mod = detrend(modh_series(m,:));

        mu_obs = mean(obs_series(m,:));
        mu_mod = mean(modh_series(m,:));
    else
        D_obs = sort(obs_series(m,obs_series(m,:)>lower_lim));
        D_mod = sort(modh_series(m,modh_series(m,:)>lower_lim));

        freq_obs = size(D_obs,2)/size(obs_series,2);
        freq_mod = size(D_mod,2)/size(modh_series(m,:),2);

        if freq_mod==0
            freq_mod = 1/(365*y_obs);
        end
    end

    % b) [Switanek et al. (2017), step 2)]
    %    For temperature, fit normal probability distribution to the
    %    detrended observed and modeled series in the historical period.
    %    For precipitation, fit gamma distribution parameters by maximum
    %    likelihood to the rainday values of the observed and modeled
    %    series in the historical period.
    %    Using these fitted distributions, get the corresponding cumulative
    %    distribution values. Tail probabilities are set to (0 + threshold)
    %    for temperature and to (1 - threshold) for temperature and
    %    precipitation. Threshold is set in the first lines of this
    %    function. Default is CDF_th = 10^-3.
    if SDM_var == 0
        mu_obsD = mean(D_obs);
        mu_modD = mean(D_mod);
        sigma_obsD = std(D_obs,0);
        sigma_modD = std(D_mod,0);

        CDF_obs = normcdf(sort(D_obs),mu_obsD,sigma_obsD);
        CDF_mod = normcdf(sort(D_mod),mu_modD,sigma_modD);
        CDF_obs(CDF_obs<CDF_th) = CDF_th;
        CDF_mod(CDF_mod<CDF_th) = CDF_th;
        CDF_obs(CDF_obs>1-CDF_th) = 1-CDF_th;
        CDF_mod(CDF_mod>1-CDF_th) = 1-CDF_th;
    else
        fit_obs = gamfit(D_obs);
        fit_mod = gamfit(D_mod);
        CDF_obs = gamcdf(D_obs,fit_obs(1),fit_obs(2));
        CDF_mod = gamcdf(D_mod,fit_mod(1),fit_mod(2));
        CDF_obs(CDF_obs>1-CDF_th) = 1-CDF_th;
        CDF_mod(CDF_mod>1-CDF_th) = 1-CDF_th;
    end
    

    % 3) Projected periods:
    for j = 0:(y_mod-y_obs)

        % c) Initialize correction array.
        corr_temp = zeros(1,(2*day_win-1)*y_obs);
        
        % d) Define projected window.
        win_series  = win_series_(m,:,j+1);
        
        % e) [Switanek et al. (2017), step 1)]
        %    For temperature, get the detrended series of the projected
        %    period.
        %    For precipitation, get the rainday values, its frequency, and
        %    expected raindays for the projected period.
        %    Get the index of the sorted detrended or rainday values.
        if SDM_var == 0
            D_win = detrend(win_series);
            exp_D = size(win_series,2);
            [~,win_argsort] = sort(D_win);
            
            mu_win = nanmean(win_series);
        else
            D_win = sort(win_series(win_series>lower_lim));
            exp_D = min(round(size(D_win,2)*freq_obs/freq_mod),size(win_series,2)); 
            [~,win_argsort] = sort(win_series);            
        end
        
        % f) [Switanek et al. (2017), step 2)]
        %    For temperature, fit normal probability distribution to the
        %    detrended values of the projected period.
        %    For precipitation, fit gamma distribution parameters by
        %    maximum likelihood to the rainday values of the projected
        %    period.
        %    Using these fitted distributions, get the corresponding
        %    cumulative distribution values. Tail probabilities are set to
        %    (0 + threshold) for temperature and to (1 - threshold) for
        %    temperature and precipitation. Threshold is set in the first
        %    lines of this function. Default is CDF_th = 10^-3.
        if SDM_var == 0
            mu_winD = nanmean(D_win);
            sigma_winD = nanstd(D_win,0);
            
            diff_win = win_series - D_win;

            CDF_win = normcdf(sort(D_win),mu_winD,sigma_winD);
            CDF_win(CDF_win<CDF_th) = CDF_th;
            CDF_win(CDF_win>1-CDF_th) = 1-CDF_th;
        else
            fit_win = gamfit(D_win);
            CDF_win = gamcdf(D_win,fit_win(1),fit_win(2));    
            CDF_win(CDF_win>1-CDF_th) = 1-CDF_th;
        end

        % g) [Switanek et al. (2017), step 3)]
        %    Get the scaling between the model projected period and
        %    historical period distributions.
        if SDM_var == 0
            SF = (sigma_obsD/sigma_modD)*(norminv(CDF_win,mu_winD,sigma_winD) - norminv(CDF_win,mu_modD,sigma_modD));
        else
            SF = gaminv(CDF_win,fit_win(1),fit_win(2))./gaminv(CDF_win,fit_mod(1),fit_mod(2));
        end
        
        % h) Interpolate observed and modeled CDF of the historical period
        %    to the length of the projected period.
        if length(D_obs)>0
            obs_cdf_intpol = interp1(linspace(1,length(D_obs),length(D_obs)),CDF_obs,linspace(1,length(D_obs),length(D_win)));
        else
            obs_cdf_intpol = CDF_win*0;
        end

        if length(D_mod)>0
            mod_cdf_intpol = interp1(linspace(1,length(D_mod),length(D_mod)),CDF_mod,linspace(1,length(D_mod),length(D_win)));
        else
            mod_cdf_intpol = CDF_win*0;
        end
        
        % i) [Switanek et al. (2017), steps 4 and 5)]
        %    Get recurrence intervals and its scaled version for the
        %    observed, historical modeled and projected period modeled
        %    CDFs.
        if SDM_var == 0
            RI_obs = 1./(0.5 - abs(obs_cdf_intpol - 0.5));
            RI_mod = 1./(0.5 - abs(mod_cdf_intpol - 0.5));
            RI_win = 1./(0.5 - abs(CDF_win - 0.5));
            RI_scaled = RI_obs.*RI_win./RI_mod;
            
            CDF_scaled = sign(obs_cdf_intpol - 0.5).*(1 - 1./RI_scaled);
            CDF_scaled(CDF_scaled<0) = CDF_scaled(CDF_scaled<0) + 1;
            
            CDF_scaled(CDF_scaled<CDF_th) = CDF_th;
            CDF_scaled(CDF_scaled>1-CDF_th) = 1-CDF_th;
        else
            RI_obs = 1./(1-obs_cdf_intpol);
            RI_mod = 1./(1-mod_cdf_intpol);
            RI_win = 1./(1-CDF_win);
            RI_scaled = RI_obs.*RI_win./RI_mod;
            RI_scaled(RI_scaled<1) = 1;

            CDF_scaled = sort(1 - 1./(RI_scaled));
        end
        
        % j) [Switanek et al. (2017), step 6)]
        %    Get the initial bias corrected values. For precipitation,
        %    these values are interpolated to the length of the expected
        %    raindays.
        if SDM_var == 0
            xvals = norminv(sort(CDF_scaled),mu_obsD,sigma_obsD) + SF;
            xvals = xvals - mean(xvals) + mu_obs + (mu_win - mu_mod);
        else
            xvals = gaminv(CDF_scaled,fit_obs(1),fit_obs(2)).*SF;
            if size(D_win,2) > exp_D
                xvals = interp1(linspace(1,length(D_win),length(D_win)),xvals,linspace(1,length(D_win),exp_D));
            else
                xvals = [zeros(1,exp_D - size(D_win,2)) xvals];
            end
        end
        
        % k) [Switanek et al. (2017), step 7)]
        %    Bias corrected values are placed back matching the higher bias
        %    corrected values with the higher rainday or detrended values.
        %    For temperature, the trend of the projected period is added
        %    back.
        corr_temp(win_argsort(end-exp_D+1:end)) = xvals;
        if SDM_var == 0
            corr_temp = corr_temp + diff_win - mu_win;
        end
        
        % l) If the projected period is the historical period (period 0,
        %    j=0) save the complete bias corrected series.
        %    If the projected period is not the historical period (j>0),
        %    save the value of the last year.
        if j == 0
            SDM_h(m,:) = corr_temp(day_win:(2*day_win-1):end);
        else
            SDM(m,j) = corr_temp(end);
        end
    end
end

% Reshape to a column vector
SDM = SDM(:);
SDM_h = SDM_h(:);

% Concat historical and future
SDM_series = [SDM_h' SDM']';

% If precipitation, replace no-rain values with 0
if SDM_var==1
    SDM_series(SDM_series<pp_threshold) = 0;
end
end