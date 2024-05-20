function [years, series] = formatQM(series_, allow_negatives, frq, pp_threshold, pp_factor)
%% formatQM:
%  This function formats a time series from a column vector to a matrix.
%  The number of rows depends on the frequency of the data (1, 12 and 365
%  for annual, monthyl and daily data, respectively). If negative values
%  are not allowed, determined by allow_negatives=1, no-rain values (values
%  below pp_threshold) are replaced by random positive values below the
%  product pp_threshold*pp_factor.
%
% Input:
%   series_ = A column vector of daily data, without considering leap days.
%             The length of the column vector should by a multiple of 365.
%             [ndata_obs, 1]
%
%   allow_negatives = A flag that identifies if data allows negative values
%                     and also to replace no-rain values with random small
%                     values (Chadwick et al., 2023) to avoid numerical
%                     problems with the probability distribution functions.
%                     allow_negatives = 1 or True: Allow negatives
%                     allow_negatives = 0 or False: Do not allow negative
%
%   frq = A string specifying if the input frequency is daily, monthly or
%         annual.
%         Daily:     frq = 'D'
%         Monthly:   frq = 'M'
%         Annual:    frq = 'A'
%
%   pp_threshold = A float indicating the threshold to consider no-rain
%                  values.
%
%   pp_factor = A float which multiplied to pp_threshold indicates the
%               maximum value of the random values that replace no-rain
%               values.
%
% Output:
%   years = Number of years of the series.
%
%   series_matrix = A matrix of daily, monthly or annual observed data,
%                   without considering leap days. For daily, monthly and
%                   annual frequency, the number of rows will be 365, 12
%                   and 1, respectively. The number of columns will be
%                   equal to the number of years. [1, 12 or 365, years]
%
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

% This prevents modyfing the original input passed by reference
% Must look for a more elegant way to skip this step
series = series_;

% Check frequency and define the number of rows
if frq=='D'
    I = 365;
elseif frq=='M'
    I = 12;
elseif frq=='A'
    I = 1;
else
    I = 1;
end

% If variable is precipitation, determined by allow_negatives=1, replace 
% low values with random values.
if allow_negatives == 0
    bool_low = series<pp_threshold;
    series(bool_low) = rand(size(series(bool_low)))*pp_threshold*pp_factor;
end

% Get number of years of the series.
years = length(series)/I;

% Reshape to get a matrix of shape [I , years]
series = reshape(series,I,[]);

end