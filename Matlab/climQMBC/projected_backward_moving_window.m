function win_series = projected_backward_moving_window(series, projected_win, frq)
%% formatQM:
% This function formats the inputs and gets basic statistics for the
% different Quantile Mapping (QM, DQM, QDM, UQM and SDM) methods available
% in the climQMBC package. If monthly data is specified, the input series
% will be reshaped to a matrix of 12 rows and several columns equal to the
% number of years of each series. If annual data is specified, the input
% is reshaped to a row vector with same entries as the input series. For
% precipitation, physically null values (values below pp_threshold) are
% replaced by random positive values below pp_factor.
%
% Description:
%   0) Check if annual or monthly data is specified.
%
%   1) If variable is precipitation, replace low values with random values.
%
%   2) Get number of years of the observed period.
%
%   3) If monthly data is specified, reshape the input series to a matrix
%      of 12 rows and several columns equal to the number of years of each
%      series. If annual data is specified, reshape the input to a row
%      vector with same entries as the input series.
%
%   4) If monthly data is specified, get monthly mean, standard deviation,
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
%   var = A flag that identifies if data are temperature or precipitation.
%         This flag tells the getDist function if it has to discard
%         distribution functions that allow negative numbers, and if the 
%         terms in the correction equations are multiplied/divided or
%         added/subtracted.
%         Temperature:   var = 0
%         Precipitation: var = 1
%
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

if frq == 'D'
    y_mod = size(series,3);
    day_win = (size(series,2)+1)/2;

    win_series = cat(3, repmat(series, [1,1,projected_win]),zeros(size(series,1), day_win*2-1, projected_win));
    win_series = reshape(win_series, [size(series,1), day_win*2-1, y_mod+1, projected_win]);
    win_series = win_series(:,:,2:(y_mod-projected_win+1),:);
else
    y_mod = size(series,2);

    win_series = [repmat(series,1,projected_win), zeros(size(series,1),projected_win)];
    win_series = reshape(win_series,[size(series,1),y_mod+1,projected_win]);
    win_series = win_series(:,2:end-projected_win,:);

end