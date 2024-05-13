function  series_moving = day_centered_moving_window(series, day_win)
%% day_centered_moving_window:
%  This function transform a 2D matrix of dimension [days in year, years]
%  to a 3D array, where the new dimension is a centered moving window for
%  each day of the year.
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
series_moving = cat(1,series(end-day_win+1:end,:),repmat(series, [day_win*2,1]),series(1:day_win,:));
series_moving = reshape(series_moving,[size(series,1)+1, day_win*2, size(series,2)]);
series_moving = series_moving(1:end-1,2:end,:);

end