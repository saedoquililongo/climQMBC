function [y_obs,obs_series,mod_series,mu_obs,sigma_obs,skew_obs,skewy_obs,mu_mod,sigma_mod,skew_mod,skewy_mod] = formatQM(obs,mod,frq)
%% formatQM:
% This function formats the inputs and gets basic statistics for the
% different Quantile Mapping (QM, DQM, QDM, UQM and SDM) methods available
% in the climQMBC package. If monthly data is specified, the input series
% will be reshaped to a matrix of 12 rows and several columns equal to the
% number of years of each series. If annual data is specified, the input
% is reshaped to a row vector with same entries as the input series.
%
% Description:
%   0) Check if annual or monthly data is specified.
%
%   1) Get number of years of the observed period.
%
%   2) If monthly data is specified, reshape the input series to a matrix
%      of 12 rows and several columns equal to the number of years of each
%      series. If annual data is specified, reshape the input to a row
%      vector with same entries as the input series.
%
%   3) If monthly data is specified, get monthly mean, standard deviation,
%      skewness, and log-skewness for the historical period of the observed
%      and modeled series. If annual data is specified, get monthly mean,
%      standard deviation, skewness, and log-skewness for the historical
%      period of the observed and modeled series.
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
%   frq = A string specifying if the input is annual or monthly data. If
%         not specified, it is set monthly as default.
%         Monthly:   var = 'M'
%
%   NOTE: This routine considers that obs and mod series start in january
%         of the first year and ends in december of the last year.
%
% Output:
%   y_obs = Number of observed years
%
%   obs_series = A column vector of monthly or annual observed data 
%               (temperature or precipitation). If monthly frequency is
%               specified, the length of this vector is 12 times the number
%               of observed years [12, y_obs]. If annual frequency is
%               specified, the length of this vector is equal to the number
%               of observed years [1, y_obs].
%
%   mod_series = A column vector of monthly or annual modeled data 
%               (temperature or precipitation). If monthly frequency is
%               specified, the length of this vector is 12 times the number
%               of observed years [12, y_mod]. If annual frequency is
%               specified, the length of this vector is equal to the number
%               of observed years [1, y_mod].
%
%   mu_obs = If monthly frequency is specified, a column vector of monthly 
%            mean of observed data [12,1]. If annual frequency is
%            specified, the mean of the observed data (float).
%
%   mu_mod = If monthly frequency is specified, a column vector of monthly 
%            mean of modeled data of the historical period [12,1]. If
%            annual frequency is specified, the mean of the modeled data
%            of the historical period(float).
%
%   sigma_obs = If monthly frequency is specified, a column vector of 
%               monthly standard deviation of observed data [12,1]. If
%               annual frequency is specified, the standard deviation of
%               the observed data (float).
%
%   sigma_mod = If monthly frequency is specified, a column vector of 
%               monthly standard deviation of modeled data of the
%               historical period [12,1]. If annual frequency is
%               specified, the standard deviation of the modeled data of
%               the historical period(float).
%
%   skew_obs = If monthly frequency is specified, a column vector of 
%               monthly skewness of observed data [12,1]. If annual
%               frequency is specified, the skewness of the observed data
%               (float).
%
%   skew_mod = If monthly frequency is specified, a column vector of 
%               monthly skewness of modeled data of the historical period
%               [12,1]. If annual frequency is specified, the skewness of
%               the modeled data of the historical period(float).
%
%   skewy_obs = If monthly frequency is specified, a column vector of 
%               monthly skewness of the logarithm of observed data [12,1].
%               If annual frequency is specified, the skewness of the
%               logarithm of the observed data (float).
%
%   skewy_mod = If monthly frequency is specified, a column vector of 
%               monthly skewness of the logarithm of modeled data of the
%               historical period [12,1]. If annual frequency is
%               specified, the skewness of the logarithm of the modeled
%               data of the historical period(float).
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
%   Maintainer contact: slaedo@uc.cl
% Revision: 0, updated Dec 2021

%%
% 0) Check if annual or monthly data is specified.
if frq == 'A'
    I = 1;
else
    I = 12;
end

% 1) Get number of years of the observed period.
y_obs=length(obs)/I;

% 2) If monthly data is specified, reshape the input series to a matrix of
%    12 rows and several columns equal to the number of years of each 
%    series. If annual data is specified, reshape the input to a row 
%    vector with same entries as the input series.
obs_series = reshape(obs,I,[]);
mod_series = reshape(mod,I,[]);

% 3) If monthly data is specified, get monthly mean, standard deviation, 
%    skewness, and log-skewness for the historical period of the observed
%    and modeled series. If annual data is specified, get monthly mean,
%    standard deviation, skewness, and log-skewness for the historical
%    period of the observed and modeled series.
mu_obs  = nanmean(obs_series,2);     % Mean
sigma_obs = nanstd(obs_series,0,2);  % Standard deviation
skew_obs = skewness(obs_series,0,2); % Skewness
Ln_obs  = log(obs_series);
Ln_obs(imag(Ln_obs)~=0) = 0;
Ln_obs(isinf(Ln_obs)) = log(0.01);
skewy_obs = skewness(Ln_obs,0,2);    % Log-skewness

mu_mod  = nanmean(mod_series(:,1:y_obs),2);     % Mean
sigma_mod = nanstd(mod_series(:,1:y_obs),0,2);  % Standard deviation
skew_mod = skewness(mod_series(:,1:y_obs),0,2); % Skewness
Ln_mod  = log(mod_series(:,1:y_obs));
Ln_mod(imag(Ln_mod)~=0) = 0;
Ln_mod(isinf(Ln_mod)) = log(0.01);
skewy_mod = skewness(Ln_mod,0,2);               % Log-skewness
end