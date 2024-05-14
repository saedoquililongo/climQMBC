function win_series = projected_backward_moving_window(series, projected_win, frq)
%% projected_backward_moving_window:
%  This function transform a 2D or 3D array of dimensions [days in year, year]
%  or [days in year, centered window for each day, years] to a 3D or 4D
%  array adding a new dimension to represent a backward moving window for
%  each projected period.
%
% Input:
%   series = A 2D or 3D array of daily, monthly or annual observed data,
%            without considering leap days. For daily, monthly and annual
%            frequency, the number of rows will be 365, 12 and 1,
%            respectively. If daily data, it is a 3D array that has a 
%            dimension with a centered moving window for each day of the
%            year. [1, 12 or 365, years] or [365, day centered window, years]
%
%   projected_win = An integer the length of the backward moving window.
%
%   frq = A string specifying if the input frequency is daily, monthly or
%         annual.
%         Daily:    frq = 'D'
%         Monthly:  frq = 'M'
%         Annual:   frq = 'A'
%
% Output:
%   win_series = An array of daily, monthly or annual future data, without
%                considering leap days, with a dimension associated to a
%                backward moving window for each  projected period.
%                Possible dimensions:
%                3D array: [sub-periods, projected periods window, years]
%                4D array: [sub-periods, days window, projected periods window, years]
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