function UQM_series = UQM(obs, mod, mult_change, allow_negatives, frq,pp_threshold,pp_factor,day_win, user_pdf, pdf_obs, pdf_mod)
%% UQM_series:
%   This function performs bias correction of modeled series based on
%   observed data by the Unbiased Quantile Mapping (UQM) method, as 
%   described by Chadwick et al. (2022) . Correction is performed to  
%   monthly or annual precipitation or temperature data in a single 
%   location. An independent probability distribution function is assigned 
%   to each month and to each projected period based on the  
%   Kolmogorov-Smirnov test. If annual frequency is specified, this is 
%   applied to the complete period. Each projected period considers a 
%   window whose length is equal to the  number of years of the historical 
%   period and ends in the analyzed year.
%
%   Correction of the historical period is performed by the Quantile 
%   Mapping (QM) method, as described by Cannon et al. (2015). 
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
%   3) For each projected period, get the delta factor (delta) and time 
%      dependent (aster) statistics (mean, standard deviation, skewness, 
%      and log skewness).
%
%   4) For each projected period:
%      a) Assign a probability distribution function to each month. If 
%         annual frequency is specified, this is applied to the complete 
%         period.
%      b) Apply the cumulative distribution function of the projected 
%         period, evaluated with the statistics of this period, to the last
%         data of the period.
%
%   5) Apply the inverse cumulative distribution function of the observed 
%      data, evaluated with the time dependent statistics, to the values 
%      obtained in 4b).
%
%   6) Perform QM for the historical period.
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
%         annual:   frq = 'A'
% 
%   pp_threshold = A float indicating the threshold to consider physically 
%                  null precipitation values.
%
%   pp_factor = A float indicating the maximum value of the random values
%               that replace physically null precipitation values.
%
% Output:
%   UQM_series = A column vector of monthly or annual modeled data 
%              (temperature or precipitation) corrected by the UQM method.
%              If monthly frequency is specified, the length of this vector
%              is 12 times the number of observed years [12 x y_mod, 1]. If
%              annual frequency is specified, the length of this vector
%              is equal to the number of observed years [y_mod, 1].
%
% References:
%   Chadwick et al. (2022) [under revision]
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
if ~exist('mult_change','var')
    mult_change = 1;
end

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

    win_series = projected_backward_moving_window(mod_series_moving, y_obs, frq);

    obs_series_moving = reshape(permute(obs_series_moving, [1,3,2]),[365,(2*day_win-1)*y_obs]);
    modh_series_moving = reshape(permute(mod_series_moving(:,:,1:y_obs), [1,3,2]),[365,(2*day_win-1)*y_obs]);

    win_series_moving = reshape(permute(win_series,[1,4,2,3]), [365, (2*day_win-1)*y_obs, y_mod-y_obs]);
    
    if allow_negatives==false
        obs_series_moving = set_norain_to_nan(obs_series_moving, pp_threshold, pp_factor);
        modh_series_moving = set_norain_to_nan(modh_series_moving, pp_threshold_mod, pp_factor);
        mod_series(mod_series<pp_threshold_mod) = NaN;
        win_series_moving = set_norain_to_nan(win_series_moving, pp_threshold_mod, pp_factor);
    end

    [mu_obs, std_obs, skew_obs, skewy_obs] = getStats(obs_series_moving, frq);
    [mu_mod, std_mod, skew_mod, skewy_mod] = getStats(modh_series_moving, frq);
    [mu_win, std_win, skew_win, skewy_win] = getStats(win_series_moving, frq);

    mu_win = reshape(mu_win, [size(mu_win,1),size(mu_win,3)]);
    std_win = reshape(std_win, [size(std_win,1),size(std_win,3)]);
    skew_win = reshape(skew_win, [size(skew_win,1),size(skew_win,3)]);
    skewy_win = reshape(skewy_win, [size(skewy_win,1),size(skewy_win,3)]);
    
else
    win_series = projected_backward_moving_window(mod_series, y_obs, frq);

    [mu_obs, std_obs, skew_obs, skewy_obs] = getStats(obs_series, frq);
    [mu_mod, std_mod, skew_mod, skewy_mod] = getStats(mod_series(:,1:y_obs), frq);
    [mu_win, std_win, skew_win, skewy_win] = getStats(win_series, frq);
end

% 2) Assign a probability distribution function to each month of the 
%    historical period (getDist function of the climQMBC package).
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



mu_mod_repeated = repmat(mu_mod,1,(y_mod-y_obs));
std_mod_repeated = repmat(std_mod,1,(y_mod-y_obs));
skew_mod_repeated = repmat(skew_mod,1,(y_mod-y_obs));
skewy_mod_repeated = repmat(skewy_mod,1,(y_mod-y_obs));

mu_obs_repeated = repmat(mu_obs,1,(y_mod-y_obs));
std_obs_repeated = repmat(std_obs,1,(y_mod-y_obs));
skew_obs_repeated = repmat(skew_obs,1,(y_mod-y_obs));
skewy_obs_repeated = repmat(skewy_obs,1,(y_mod-y_obs));


if mult_change==1 % Precipitation
    delta_mu = mu_win./mu_mod_repeated;
    delta_sigma = std_win./std_mod_repeated;
    delta_skew = skew_win./skew_mod_repeated;
    delta_skewy = skewy_win./skewy_mod_repeated;
    
    mu_projected = mu_obs_repeated.*delta_mu;
    sigma_projected = std_obs_repeated.*delta_sigma;
    skew_projected = skew_obs_repeated.*delta_skew;
    skewy_projected = skewy_obs_repeated.*delta_skewy;

else % Temperature
    delta_mu = mu_win - mu_mod_repeated;
    delta_sigma = std_win - std_mod_repeated;
    delta_skew = skew_win - skew_mod_repeated;
    delta_skewy = skewy_win - skewy_mod_repeated;
    
    mu_projected = mu_obs_repeated + delta_mu;
    sigma_projected = std_obs_repeated + delta_sigma;
    skew_projected = skew_obs_repeated + delta_skew;
    skewy_projected = skewy_obs_repeated + delta_skewy;
end

% 4) For each projected period:
ks_fail_win = 0;
pdf_win = zeros(size(mod_series,1),size(mod_series,2)-y_obs);
prob = zeros(size(mod_series,1),size(mod_series,2)-y_obs);
UQM  = zeros(size(mod_series,1),size(mod_series,2)-y_obs);
for j = 1:size(prob,2)
    
    % a) Assign a probability distribution function to each month. If
    %    annual frequency is specified, this is applied to the complete 
    %    period (getDist function of the climQMBC package).
    if user_pdf==false
        if frq=='D'
            [pdf_win(:,j), ks_fail_wtemp] = getDist(win_series_moving(:,:,j),allow_negatives,mu_win(:,j),std_win(:,j),skew_win(:,j),skewy_win(:,j));
        else
            [pdf_win(:,j), ks_fail_wtemp] = getDist(reshape(win_series(:,j,:),size(win_series,[1,3])),allow_negatives,mu_win(:,j),std_win(:,j),skew_win(:,j),skewy_win(:,j));
        end
        ks_fail_win = ks_fail_win + ks_fail_wtemp;

    else
        pdf_win(:,j) = pdf_mod;
    end

    
    
    % b) Apply the cumulative distribution function of the projected
    %    period, evaluated with the statistics of this period, to the last
    %    data of the period (getCDF function of the climQMBC package).
    %    Equation 3 of Cannon et al. (2015).
    if frq=='D'
        prob(:,j) = getCDF(pdf_win(:,j),mod_series(:,y_obs+j),mu_win(:,j),std_win(:,j),skew_win(:,j),skewy_win(:,j));
    else
        prob(:,j) = getCDF(pdf_win(:,j),win_series(:,j,end),mu_win(:,j),std_win(:,j),skew_win(:,j),skewy_win(:,j));
    end

    % 5) Apply the inverse cumulative distribution function of the observed
    %    data, evaluated with the time dependent statistics, to the values
    %    obtained in 4b) (getCDFinv function of the climQMBC package). Equation
    %    X of Chadwick et al. (2021).

    UQM(:,j) = getCDFinv(pdf_obs,prob(:,j),mu_projected(:,j),sigma_projected(:,j),skew_projected(:,j),skewy_projected(:,j));
end

UQM = UQM(:);

% 6) Perform QM for the historical period
mod_h = mod(1:length(obs));
mod_h = mod_h(:);
QM_series = QM(obs,mod_h,allow_negatives,frq,pp_threshold,pp_factor,day_win,user_pdf,pdf_obs,pdf_mod);
UQM_series = [QM_series' UQM']';
if allow_negatives==0
    UQM_series(isnan(UQM_series)) = 0;
    UQM_series(UQM_series<pp_threshold) = 0;
end

if user_pdf==false
    ks_fail = ks_fail_obs + ks_fail_win;
    if ks_fail>0
        disp('UQM: Some of the probability distribution functions did not pass the KS-Test')
    end
end
end