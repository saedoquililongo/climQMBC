function pp_threshold_mod = get_pp_threshold_mod(obs,mod,pp_threshold)
%% get_pp_threshold_mod:
%  This functions gets a no-rain value threshold for the modeled series by 
%  mathcing the no-rain days of the modeled series in the historical period
%  with the observed series.
%
% Input:
%   obs = A column vector of daily observed data, without considering leap 
%         days. The length of the column vector should by a multiple of 365.
%         [ndata_obs, 1]
%
%   mod = A column vector of daily modeled data, without considering leap
%         days. The length of the column vector should by a multiple of 365.
%         [ndata_mod, 1]
%
%   pp_threshold = A float indicating the threshold to consider no-rain 
%                  values in the observed data.
%
% Output:
%   pp_threshold_mod = A float indicating the threshold to consider no-rain
%                      values in the modeled data.
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

%%
% Get the rain days of the observed series
obs_rainday_hist = sum(obs>pp_threshold);

% Define the pp_threshold of the modeled series by matching the rain days
% of the observed series
mod_sort_descending = sort(mod(1:length(obs)), 'descend');
days_kept = min(obs_rainday_hist, length(obs));

if days_kept~=length(obs)
    pp_threshold_mod = mod_sort_descending(days_kept+1);
else
    pp_threshold_mod = 0;
end
end