function SDM_series = SDM(obs,mod,SDM_var,frq,pp_threshold,pp_factor)
%% SDM_series
%   This function performs bias correction of modeled series based on
%   observed data by the Scaled Distribution Mapping (SDM) method, as
%   described by Switanek et al. (2017). Correction is performed to monthly
%   or annual precipitation or temperature data in a single location. The
%   historical period of the modeled series is bias corrected as a scenario
%   where the complete series is replaced by the bias corrected series. On
%   the other hand, for each year after the historical period the method
%   considers a projected period as a window whose length is equal to the
%   number of years of the historical period and ends in the analyzed year.
%   From the bias corrected series, only the last year is saved and all the
%   others are discarded. This approach aims to build a continuum for the
%   bias corrected series.
%
% Description:
%   0) Check if annual or monthly data is specified.
%
%   1) Format inputs and get statistics of the observed and modeled series 
%      of the historical period.
%
%   2) Historical period:
%      a) [Switanek et al. (2017), step 1)]
%          For temperature, get the detrended modeled and observed series
%          in the historical period.
%          For precipitation, get the rainday values and its frequency for
%          the modeled and observed series in the historical period.
%
%      b) [Switanek et al. (2017), step 2)]
%         For temperature, fit normal probability distribution to the 
%         detrended observed and modeled series in the historical period.
%         For precipitation, fit gamma distribution parameters by maximum
%         likelihood to the rainday values of the observed and modeled 
%         series in the historical period.
%         Using these fitted distributions, get the corresponding 
%         cumulative distribution values. Tail probabilities are set to 
%         (0 + threshold) for temperature and to (1 - threshold) for
%         temperature and precipitation. Threshold is set in the first 
%         lines of this function. Default is CDF_th = 10^-3.
%
%   3) Projected periods:
%      c) Initialize correction array.
%      d) Define projected window.
%      e) [Switanek et al. (2017), step 1)]
%         For temperature, get the detrended series of the projected 
%         period.
%         For precipitation, get the rainday values, its frequency, and 
%         expected raindays for the projected period.
%         Get the index of the sorted detrended or rainday values.
%      f) [Switanek et al. (2017), step 2)]
%         For temperature, fit normal probability distribution to the
%         detrended values of the projected period.
%         For precipitation, fit gamma distribution parameters by maximum
%         likelihood to the rainday values of the projected period.
%         Using these fitted distributions, get the corresponding 
%         cumulative distribution values. Tail probabilities are set to
%         (0 + threshold) for temperature and to (1 - threshold) for 
%         temperature and precipitation. Threshold is set in the first 
%         lines of this function. Default is CDF_th = 10^-3.
%      g) [Switanek et al. (2017), step 3)]
%         Get the scaling between the model projected period and historical
%         period distributions.
%      h) Interpolate observed and modeled CDF of the historical period to
%         the length of the projected period.
%      i) [Switanek et al. (2017), steps 4 and 5)]
%         Get recurrence intervals and its scaled version for the observed,
%         historical modeled and projected period modeled CDFs.
%      j) [Switanek et al. (2017), step 6)]
%         Get the initial bias corrected values. For precipitation, these
%         values are interpolated to the length of the expected raindays.
%      k) [Switanek et al. (2017), step 7)]
%         Bias corrected values are placed back matching the higher bias
%         corrected values with the higher rainday or detrended values.
%         For temperature, the trend of the projected period is added back.
%      l) If the projected period is the historical period save the
%         complete bias corrected series.
%         If the projected period is not the historical period, save the
%         value of the last year. 
%
% Input:
%   obs = A column vector of monthly or annual observed data (temperature
%         or precipitation). If monthly frequency is specified, the length
%         of this vector is 12 times the number of observed years
%         [12 x y_obs, 1]. If annual frequency is specified, the length 
%         of this vector is equal to the number of observed years
%         [y_obs, 1].
%
%   mod = A column vector of monthly or annual modeled data (temperature
%         or precipitation). If monthly frequency is specified, the length
%         of this vector is 12 times the number of observed years
%         [12 x y_mod, 1]. If annual frequency is specified, the length
%         of this vector is equal to the number of observed years
%         [y_mod, 1].
%
%   var = A flag that identifies if data are temperature or precipitation.
%         This flag tells the getDist function if it has to discard
%         distribution functions that allow negative numbers, and if the 
%         terms in the correction equations are multiplied/divided or
%         added/subtracted.
%         Temperature:   var = 0
%         Precipitation: var = 1
%
%   NOTE: This routine considers that obs and mod series start in january
%         of the first year and ends in december of the last year.
%
% Optional inputs:
%   frq = A string specifying if the input is annual or monthly data. If
%         not specified, it is set monthly as default.
%         Monthly:    frq = 'M'
%         annual:     frq = 'A'
% 
%   pp_threshold = A float indicating the threshold to consider physically 
%                  null precipitation values.
%
%   pp_factor = A float indicating the maximum value of the random values
%               that replace physically null precipitation values.
%
% Output:
%   SDM_series = A column vector of monthly or annual modeled data 
%                (temperature or precipitation) corrected by the SDM
%                 method. If monthly frequency is specified, the length of 
%                this  vector is 12 times the number of observed years 
%                [12 x y_mod, 1]. If annual frequency is specified, the 
%                length of this vector is equal to the number of observed
%                years [y_mod, 1].
%
% References:
%   Switanek, B. M., Troch, P., Castro, L. C., Leuprecht, A., Chang, H. I.,
%   Mukherjee, R., and Demaria, M. C. E. (2017) Scaled distribution
%   mapping: A bias correction method that preserves raw climate model
%   projected changes. Hydrology &amp; Earth System Sciences, 21,
%   2649-2666, https://doi.org/10.5194/hess-21-2649-2017.
%

% Written by Sebastian Aedo Quililongo (1*)
%            Cristian Chadwick         (2)
%            Fernando Gonzalez-Leiva   (3)
%            Jorge Gironas             (3)
%            
%   (1) Centro de CAmbio Global UC, Pontificia Universidad Catolica de 
%       Chile, Santiago, Chile
%   (2) Faculty of Engineering and Sciences, Universidad Adolfo Ibanez,
%       Santiago, Chile
%   (3) Department of Hydraulics and Environmental Engineering, Pontificia
%       Universidad Catolica de Chile, Santiago, Chile
%       
%   *Maintainer contact: slaedo@uc.cl
% Revision: 1, updated Jul 2022


%%

% Define optional arguments
if ~exist('pp_threshold','var')
    pp_threshold = 1;
end

if ~exist('pp_factor','var')
    pp_factor = 1/100;
end

lower_lim = pp_threshold;
CDF_th = 10^-3;

% 0) Check if annual or monthly data is specified.
if ~exist('frq','var')
    frq = 'M';
end

% 1) Format inputs and get statistics of the observed and modeled series of
% the historical period (formatQM function of the climQMBC package).
[y_obs,obs_series] = formatQM(obs, SDM_var, frq,pp_threshold,pp_factor);
[y_mod,mod_series] = formatQM(mod, SDM_var, frq,pp_threshold,pp_factor);

SDM  = zeros(size(mod_series,1),size(mod_series,2)-y_obs);
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
        D_mod = detrend(mod_series(m,1:y_obs));

        mu_obs = mean(obs_series(m,:));
        mu_mod = mean(mod_series(m,1:y_obs));
    else
        D_obs = sort(obs_series(m,obs_series(m,:)>lower_lim));
        D_mod = sort(mod_series(m,mod_series(m,1:y_obs)>lower_lim));

        freq_obs = size(D_obs,2)/size(obs_series,2);
        freq_mod = size(D_mod,2)/size(mod_series(m,1:y_obs),2);
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
    for j = 0:length(mod_series)-y_obs
        % c) Initialize correction array.
        corr_temp = zeros(1,y_obs);
        
        % d) Define projected window.
        win_series  = mod_series(m,1+j:y_obs+j);
        
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
        obs_cdf_intpol = interp1(linspace(1,length(D_obs),length(D_obs)),CDF_obs,linspace(1,length(D_obs),length(D_win)));
        mod_cdf_intpol = interp1(linspace(1,length(D_mod),length(D_mod)),CDF_mod,linspace(1,length(D_mod),length(D_win)));
        
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
            SDM_h(m,:) = corr_temp;
        else
            SDM(m,j) = corr_temp(end);
        end
    end
end

SDM = SDM(:);
SDM_h = SDM_h(:);
SDM_series = [SDM_h' SDM']';
if SDM_var==1
    SDM_series(SDM_series<pp_threshold) = 0;
end
end