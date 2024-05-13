function QM_series = QM(obs, mod, allow_negatives, frq,pp_threshold,pp_factor, day_win, user_pdf, pdf_obs, pdf_mod)
%% QM_series:
%   This function performs bias correction of modeled series based on 
%   observed data by the Quantile Mapping (QM) method, as described by
%   Cannon et al. (2015) . Correction is performed to monthly or annual
%   precipitation or temperature data in a single location. An independent
%   probability distribution function is assigned to each month and to each
%   projected period based on the Kolmogorov-Smirnov test. If annual 
%   frequency is specified, this is applied to the complete period. 
%   Each projected period considers a window whose length is equal to the 
%   number of years of the historical period and ends in the analyzed year.
%
% Description:
%   0) Check if annual or monthly data is specified.
%
%   1) Format inputs and get statistics of the observed and modeled series 
%      of the historical period.
%
%   2) Assign a probability distribution function to each month for the
%      observed and modeled data in the historical period. If annual
%      frequency is specified, this is applied to the complete historical
%      period.
%
%   3) Apply the cumulative distribution function of the modeled data, 
%      evaluated with the statistics of the modeled data in the historical 
%      period, to the modeled data.
%
%   4) Apply the inverse cumulative distribution function of the observed
%      data, evaluated with the statistics of the observed data in the
%      historical period, to the probabilities obtained from 3). 
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
%   QM_series = A column vector of monthly or annual modeled data 
%              (temperature or precipitation) corrected by the QM method.
%              If monthly frequency is specified, the length of this vector
%              is 12 times the number of observed years [12 x y_mod, 1]. If
%              annual frequency is specified, the length of this vector
%              is equal to the number of observed years [y_mod, 1].
%
% References:
%   Cannon, A. J., S. R. Sobie, and T. Q. Murdock, 2015: Bias correction of
%   GCM precipitation by quantile mapping: How well do methods preserve
%   changes in quantiles and extremes? J. Climate, 28(17), 6938-6959,
%   https://doi.org/10.1175/JCLI-D-14-00754.1
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

% 0) Check if annual or monthly data is specified.
if ~exist('allow_negatives','var')
    allow_negatives = 1;
end

if ~exist('frq','var')
    frq = 'A';
end

% Define optional arguments
if ~exist('pp_threshold','var')
    pp_threshold = 1;
end

if ~exist('pp_factor','var')
    pp_factor = 1/100;
end

if ~exist('day_win','var')
    day_win = 1;
end

if ~exist('user_pdf','var')
    user_pdf = false;
    pdf_obs = false;
    pdf_mod = false;
end

if (frq=='D' & allow_negatives==false)
    pp_threshold_mod = get_pp_threshold_mod(obs, mod, pp_threshold);
else
pp_threshold_mod = pp_threshold;
end

% 1) Format inputs and get statistics of the observed and modeled series of
%    the historical period (formatQM function of the climQMBC package).
[y_obs,obs_series] = formatQM(obs, allow_negatives, frq,pp_threshold,pp_factor);
[y_mod,mod_series] = formatQM(mod, allow_negatives, frq,pp_threshold_mod,pp_factor);

if frq=='D'
    obs_series_moving = day_centered_moving_window(obs_series, day_win);
    mod_series_moving = day_centered_moving_window(mod_series, day_win);

    obs_series_moving = reshape(permute(obs_series_moving, [1,3,2]),[365,(2*day_win-1)*y_obs]);
    modh_series_moving = reshape(permute(mod_series_moving(:,:,1:y_obs), [1,3,2]),[365,(2*day_win-1)*y_obs]);
    
    if allow_negatives==false
        obs_series_moving = set_norain_to_nan(obs_series_moving, pp_threshold, pp_factor);
        modh_series_moving = set_norain_to_nan(modh_series_moving, pp_threshold_mod, pp_factor);
        mod_series(mod_series<pp_threshold_mod) = NaN;
    end

    [mu_obs, std_obs, skew_obs, skewy_obs] = getStats(obs_series_moving, frq);
    [mu_mod, std_mod, skew_mod, skewy_mod] = getStats(modh_series_moving, frq);
    
else
    [mu_obs, std_obs, skew_obs, skewy_obs] = getStats(obs_series, frq);
    [mu_mod, std_mod, skew_mod, skewy_mod] = getStats(mod_series(:,1:y_obs), frq);
end



% 2) Assign a probability distribution function to each month for the
%    observed and modeled data in the historical period. If annual
%    frequency is specified, this is applied to the complete historical
%    period (getDist function of the climQMBC package).
if user_pdf==false
    if frq=='D'
        [pdf_obs, ks_fail_obs] = getDist(obs_series_moving,allow_negatives,mu_obs,std_obs,skew_obs,skewy_obs);
        [pdf_mod, ks_fail_mod] = getDist(modh_series_moving,allow_negatives,mu_mod,std_mod,skew_mod,skewy_mod);
    else
        [pdf_obs, ks_fail_obs] = getDist(obs_series,allow_negatives,mu_obs,std_obs,skew_obs,skewy_obs);
        [pdf_mod, ks_fail_mod] = getDist(mod_series(:,1:y_obs),allow_negatives,mu_mod,std_mod,skew_mod,skewy_mod);
    end
else
    
    pdf_obs = zeros(size(obs_series,1),1) + pdf_obs;
    pdf_mod = zeros(size(mod_series,1),1) + pdf_mod;
end



% 3) Apply the cumulative distribution function of the modeled data,
%    evaluated with the statistics of the modeled data in the historical
%    period, to the modeled data (getCDF function of the climQMBC package).
%    Equation 1 of Cannon et al. (2015).
prob = getCDF(pdf_mod,mod_series,mu_mod,std_mod,skew_mod,skewy_mod);

% 4) Apply the inverse cumulative distribution function of the observed
%    data, evaluated with the statistics of the observed data in the
%    historical period, to the probabilities obtained from 3) (getCDFinv
%    function of the climQMBC package). Equation 1 of Cannon et al. (2015).
QM_series = getCDFinv(pdf_obs,prob,mu_obs,std_obs,skew_obs,skewy_obs);

QM_series=QM_series(:);
if allow_negatives==0
    QM_series(isnan(QM_series)) = 0;
    QM_series(QM_series<pp_threshold) = 0;
end

if user_pdf==false
    ks_fail = ks_fail_obs + ks_fail_mod;
    if ks_fail>0
        disp('QM: Some of the probability distribution functions did not pass the KS-Test')
    end
end
end