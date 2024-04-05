function DQM_series = DQM(obs, mod, mult_change, allow_negatives, frq,pp_threshold,pp_factor)
%% DQM_series:
%   This function performs bias correction of modeled series based on
%   observed data by the Detrended Quantile Mapping (DQM) method, as 
%   described by Cannon et al. (2015) . Correction is performed to monthly 
%   or annual precipitation or temperature data in a single location. An
%   independent probability distribution function is assigned to each month
%   and to each projected period based on the Kolmogorov-Smirnov test. If 
%   annual frequency is specified, this is applied to the complete period. 
%   Each projected period considers a window whose length is equal to the 
%   number of years of the historical period and ends in the analyzed year.
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
%   3) Extract the long-term trend from the modeled data:
%      a) Get the monthly mean of the historical period. If annual
%         frequency is specified, this is applied to the complete period).
%      b) Get the monthly mean of each projected period. If annual
%         frequency is specified, this is applied to the complete period).
%
%   4)Compute the linear scaled values.
%
%   5) Apply the cumulative distribution function of the modeled data, 
%      evaluated with the statistics of the modeled period, to the future
%      modeled data.
%
%   6) Apply the inverse cumulative distribution function of the observed 
%      data, evaluated with the statistics of the observed data in the
%      historical period, to the probabilities obtained from 5).
%
%   7) Reimpose the trend to the values obtained in 6).
%
%   8) Perform QM for the historical period.
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
%   pp_threshold = A float indicating the threshold to consider physically 
%                  null precipitation values.
%
%   pp_factor = A float indicating the maximum value of the random values
%               that replace physically null precipitation values.
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
% Output:
%   DQM_series = A column vector of monthly or annual modeled data 
%              (temperature or precipitation) corrected by the DQM method.
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


% 1) Format inputs and get statistics of the observed and modeled series of
%    the historical period (formatQM function of the climQMBC package).
[y_obs,obs_series] = formatQM(obs, allow_negatives, frq,pp_threshold,pp_factor);
[y_mod,mod_series] = formatQM(mod, allow_negatives, frq,pp_threshold,pp_factor);

[mu_obs, std_obs, skew_obs, skewy_obs] = getStats(obs_series, frq);
[mu_mod, std_mod, skew_mod, skewy_mod] = getStats(mod_series(:,1:y_obs), frq);

% 2) Assign a probability distribution function to each month for the
%    observed and modeled data in the historical period. If annual
%    frequency is specified, this is applied to the complete historical
%    period (getDist function of the climQMBC package).
pdf_obs = getDist(obs_series,allow_negatives,mu_obs,std_obs,skew_obs,skewy_obs);
pdf_mod = getDist(mod_series(:,1:y_obs),allow_negatives,mu_mod,std_mod,skew_mod,skewy_mod);

% 3) Extract the long-term trend from the modeled data:
%    a) Get the monthly mean of the historical period. If annual
%       frequency is specified, this is applied to the complete period).
win_series = [repmat(mod_series,1,y_obs), zeros(size(mod_series,1),y_obs)];
win_series = reshape(win_series,[size(mod_series,1),y_mod+1,y_obs]);
win_series = win_series(:,2:end-y_obs,:);

mu_win = mean(win_series,3);

mu_mod_repeated = repmat(mu_mod,1,(y_mod-y_obs));


% 4)Compute the linear scaled values (value in square brackets in equation
%   2 of Cannon et al. (2015)).
if mult_change == 1   % Precipitation
    scaling_factor = mu_win./mu_mod_repeated;
    detrended_series = mod_series(:,y_obs+1:end)./scaling_factor;
else          % Temperature
    scaling_factor = mu_win - mu_mod_repeated;
    detrended_series = mod_series(:,y_obs+1:end) - scaling_factor;
end

% 5) Apply the cumulative distribution function of the modeled data,
%    evaluated with the statistics of the modeled period, to the future
%    modeled data (getCDF function of the climQMBC package). Equation 2 of
%    Cannon et al. (2015).
prob = getCDF(pdf_mod,detrended_series,mu_mod,std_mod,skew_mod,skewy_mod);

% 6) Apply the inverse cumulative distribution function of the observed
%    data, evaluated with the statistics of the observed data in the
%    historical period, to the probabilities obtained from 5) (getCDFinv
%    function of the climQMBC package). Equation 2 of Cannon et al. (2015).
DQM_LS = getCDFinv(pdf_obs,prob,mu_obs,std_obs,skew_obs,skewy_obs);

% 7) Reimpose the trend to the values obtained in 6). Equation 2 of Cannon
%    et al. (2015).
if mult_change == 1 % Precipitation
    DQM = DQM_LS.*scaling_factor;
else % Temperature
    DQM = DQM_LS + scaling_factor;
end

DQM = DQM(:);

% 8) Perform QM for the historical period.
mod_h = mod_series(:,1:y_obs);
mod_h = mod_h(:);
QM_series = QM(obs,mod_h,allow_negatives,frq);
DQM_series = [QM_series' DQM']';
if allow_negatives==0
    DQM_series(DQM_series<pp_threshold) = 0;
end
end
