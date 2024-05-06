function QM_series = QM(obs, mod, allow_negatives, frq,pp_threshold,pp_factor, win, user_pdf, pdf_obs, pdf_mod)
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

if ~exist('win','var')
    win = 1;
end

if ~exist('user_pdf','var')
    user_pdf = false;
    user_obs = false;
    user_mod = false;
end


% 1) Format inputs and get statistics of the observed and modeled series of
%    the historical period (formatQM function of the climQMBC package).
[y_obs,obs_series] = formatQM(obs, allow_negatives, frq,pp_threshold,pp_factor);
[y_mod,mod_series] = formatQM(mod, allow_negatives, frq,pp_threshold,pp_factor);

if frq=='D'
    obs_series_moving = cat(1,obs_series(end-win+1:end,:),repmat(obs_series, [win*2,1]),obs_series(1:win,:));
    obs_series_moving = reshape(obs_series_moving,[size(obs_series,1)+1, win*2, size(obs_series,2)]);
    obs_series_moving = obs_series_moving(1:end-1,2:end,:);

    mod_series_moving = cat(1,mod_series(end-win+1:end,:),repmat(mod_series, [win*2,1]),mod_series(1:win,:));
    mod_series_moving = reshape(mod_series_moving,[size(mod_series,1)+1, win*2, size(mod_series,2)]);
    mod_series_moving = mod_series_moving(1:end-1,2:end,:);

    [mu_obs, std_obs, skew_obs, skewy_obs] = getStats(obs_series_moving, frq);
    [mu_mod, std_mod, skew_mod, skewy_mod] = getStats(mod_series_moving(:,:,1:y_obs), frq);
    
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
        pdf_obs = getDist(reshape(obs_series_moving, 365,[]),allow_negatives,mu_obs,std_obs,skew_obs,skewy_obs);
        pdf_mod = getDist(reshape(mod_series_moving(:,:,1:y_obs), 365,[]),allow_negatives,mu_mod,std_mod,skew_mod,skewy_mod);
    else
        pdf_obs = getDist(obs_series,allow_negatives,mu_obs,std_obs,skew_obs,skewy_obs);
        pdf_mod = getDist(mod_series(:,1:y_obs),allow_negatives,mu_mod,std_mod,skew_mod,skewy_mod);
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
    QM_series(QM_series<pp_threshold) = 0;
end
end