function [mu, sigma, skew, skewy] = getStats(series, frq)
%% getStats:
%  This function computes the mean, standard deviation, skewness and
%  skewness of the logarithmic values, for each sub-period within the year,
%  according to the frequency initially defined. Statistics are computed
%  along axis 1. 
%
% Input:
%   series = An array of daily, monthly or annual data, without considering
%            leap days. Possible input dimensions:
%            2D array: [sub-periods, years]
%                      [sub-periods + days window, years]
%            3D array: [sub_periods, projected periods window, years]
%                      [sub-periods + days window, projected periods window, years]
%
% Output:
%   mu = Mean values of each sub-period (and projected period if input is a
%        3D array).
%
%   sigma = Standar deviation values of each sub-period (and projected 
%           period if input is a 3D array).
%   
%   skew = Skewness values of each sub-period (and projected period if input
%          is a 3D array).
%
%   skewy = Skewness values of the logarithm of theseries of each sub-period
%           (and projected period if input is a 3D array).
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
if frq == 'D'
    if length(size(series)) == 4
        dim_stats = [2 4];
    else
        dim_stats = 2;
    end
else
    if length(size(series)) == 3
        dim_stats = 3;
    else
        dim_stats = 2;
    end
end

% Get the mean, standard deviation, skewness and skewness of the
% logarithmic values of each year sub-period of the series.
mu  = nanmean(series,dim_stats);     % Mean
sigma = nanstd(series,0,dim_stats);  % Standard deviation
skew = skewness(series,0,dim_stats); % Skewness
series_log = log(series);
series_log(isinf(series_log)) = log(0.01);
skewy = skewness(series_log,0,dim_stats);    % Log-skewness

end