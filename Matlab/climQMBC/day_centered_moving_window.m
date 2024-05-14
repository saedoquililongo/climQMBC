function  series_moving = day_centered_moving_window(series, day_win)
%% day_centered_moving_window:
%  This function transform a 2D matrix of dimension [days in year, years]
%  to a 3D array, where the new dimension is a centered moving window for
%  each day of the year.
%
% Input:
%   series_matrix = A matrix of daily data, without considering leap days.
%                   The number of rows is 365 and number of columns is the
%                   number of years of the series. [365, years]
%
%   day_win = An integer indicating how many days to consider backwards and
%             forward to get the statistics of each calendar day. The
%             length of the window will be (2*win_day-1). For example,
%             day_win=15 -> window of 29.
%
% Output:
%   series_moving = A 3D array of daily data, without considering leap days,
%                   with a dimension that considers a centered moving window
%                   for each day of the year. [365, 2*win-1, years]
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
series_moving = cat(1,series(end-day_win+1:end,:),repmat(series, [day_win*2,1]),series(1:day_win,:));
series_moving = reshape(series_moving,[size(series,1)+1, day_win*2, size(series,2)]);
series_moving = series_moving(1:end-1,2:end,:);

end