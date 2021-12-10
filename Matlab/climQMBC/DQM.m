function DQM_series = DQM(obs,mod,var,frq)
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

% Written by Sebastian Aedo Quililongo (1)
%            Cristian Chadwick         (2)
%            Fernando Gonzalez-Leiva   (1)
%            Jorge Gironas             (1)
%            
%   (1) Pontificia Universidad Catolica de Chile, Santiago, Chile
%       Department of Environmental and Hydraulic Engineering
%   (2) Universidad Adolfo Ibanez, Santiago, Chile
%       Faculty of Engineering and Sciences
% Maintainer contact: slaedo@uc.cl
% Revision: 0, updated Dec 2021


%%
% 0) Check if annual or monthly data is specified.
if ~exist('frq','var')
    frq = 'M';
end

% 1) Format inputs and get statistics of the observed and modeled series of
%    the historical period (formatQM function of the climQMBC package).
[y_obs,obs_series,mod_series,mu_obs,std_obs,skew_obs,skewy_obs,mu_mod,std_mod,skew_mod,skewy_mod] = formatQM(obs,mod,frq);

% 2) Assign a probability distribution function to each month for the
%    observed and modeled data in the historical period. If annual
%    frequency is specified, this is applied to the complete historical
%    period (getDist function of the climQMBC package).
PDF_obs = getDist(obs_series,mu_obs,std_obs,skew_obs,skewy_obs,var);
PDF_mod = getDist(mod_series(:,1:y_obs),mu_mod,std_mod,skew_mod,skewy_mod,var);

% 3) Extract the long-term trend from the modeled data:
%    a) Get the monthly mean of the historical period. If annual
%       frequency is specified, this is applied to the complete period).
xbarmh = mean(mod_series(:,1:y_obs),2);

%    b) Get the monthly mean of each projected period. If annual
%       frequency is specified, this is applied to the complete period).
y_mod   = size(mod_series,2);
xbarmp = zeros(size(mod_series,1),(y_mod-y_obs));
for j=1:(y_mod-y_obs)
    xbarmp(:,j) = mean(mod_series(:,(j+1):(y_obs+j)),2);
end

% 4)Compute the linear scaled values (value in square brackets in equation
%   2 of Cannon et al. (2015)).
LS = zeros(size(xbarmp));
if var == 1   % Precipitation
    for m = 1:size(mod_series,1)
        LS(m,:) = (mod_series(m,y_obs+1:end).* xbarmh(m,1))./ xbarmp(m,:);
    end    
else          % Temperature
    for m = 1:size(mod_series,1)
        LS(m,:) = (mod_series(m,y_obs+1:end)+ xbarmh(m))- xbarmp(m,:);
    end     
end

% 5) Apply the cumulative distribution function of the modeled data,
%    evaluated with the statistics of the modeled period, to the future
%    modeled data (getCDF function of the climQMBC package). Equation 2 of
%    Cannon et al. (2015).
Taot = getCDF(PDF_mod,LS,mu_mod,std_mod,skew_mod,skewy_mod);

% 6) Apply the inverse cumulative distribution function of the observed
%    data, evaluated with the statistics of the observed data in the
%    historical period, to the probabilities obtained from 5) (getCDFinv
%    function of the climQMBC package). Equation 2 of Cannon et al. (2015).
DQM_LS = getCDFinv(PDF_obs,Taot,mu_obs,std_obs,skew_obs,skewy_obs);

% 7) Reimpose the trend to the values obtained in 6). Equation 2 of Cannon
%    et al. (2015).
if var == 1 % Precipitation
    DQM = DQM_LS .* (xbarmp ./ repmat(xbarmh,1,y_mod-y_obs)); 
else % Temperature
    DQM = DQM_LS + (xbarmp - repmat(xbarmh,1,y_mod-y_obs)); 
end

DQM = DQM(:);

% 8) Perform QM for the historical period.
mod_h = mod_series(:,1:y_obs);
mod_h = mod_h(:);
QM_series = QM(obs,mod_h,var,frq);
DQM_series = [QM_series' DQM']';
end  
