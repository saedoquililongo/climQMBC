function  series_moving = set_norain_to_nan(series_moving,pp_threshold,pp_factor,min_rainday)
%% set_norain_to_nan:
%  This function replace no-rain values with nans, leaving a minimum amout
%  of values (min_rainday) to fit a distribution. If fewer than min_rainday
%  values are above pp_threshold, random small values will be added.
%
% Input:
%   series_moving = A 2D or 3D array of daily, monthly or annual observed
%                   data, without considering leap days. For daily, monthly
%                   and annual frequency, the number of rows will be 365,
%                   12 and 1, respectively. If daily data, it is a 3D array
%                   that has a  dimension with a centered moving window
%                   for each day of the year. [1, 12 or 365, years] or
%                   [365, day centered window, years]
%
%   pp_factor = A float which multiplied to pp_threshold indicates the
%               maximum value of the random values that replace no-rain
%               values.
%
%   pp_threshold = A float indicating the threshold to consider no-rain
%                  values.
%
% Optional input
%   min_rainday = Minimum amount of values to keep for each sub-period and
%                 projected period, to ensure a minimum amount to fit a
%                 distribution. Default is 30
%
% Output:
%   series_moving = The input, but with no-rain values replaced with nans.
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

if ~exist('min_rainday','var')
    min_rainday = 30;
end

% Falta evaluar qué ocurre con la información de períodos proyectados
rainday_count = sum(series_moving>pp_threshold,2);

if length(size(series_moving))==2
    for per = 1:size(series_moving,1)
        bool_low = series_moving(per,:)<pp_threshold;
        replace_values_nans = rand(1,sum(bool_low))*pp_threshold*pp_factor;
        replace_values_nans(max(1, (min_rainday-rainday_count(per))):end) = nan;
        series_moving(per, bool_low) = replace_values_nans;
    end
else
    for per1 = 1:size(rainday_count,1)
        for per2 = 1:size(rainday_count,3)
            bool_low = series_moving(per1,:,per2)<pp_threshold;
            replace_values_nans = rand(1,sum(bool_low))*pp_threshold*pp_factor;
            replace_values_nans(max(1, (min_rainday-rainday_count(per1,1,per2))):end) = nan;
            series_moving(per1, bool_low, per2) = replace_values_nans;
        end
    end
end


end