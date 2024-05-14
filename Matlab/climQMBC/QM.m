function QM_series = QM(obs, mod, allow_negatives, frq,pp_threshold,pp_factor, day_win, user_pdf, pdf_obs, pdf_mod)
%% QM_series:
%  This function performs bias correction of modeled series based on
%  observed data by the Quantile Mapping (QM) method, as described by
%  Cannon et al. (2015). Correction is performed to daily, monthly or annual
%  data in a single location. An independent probability distribution
%  function is assigned to  each sub-annual period based on the
%  Kolmogorov-Smirnov test.
%
% Input:
%   obs = A column vector of daily, monthly or annual observed data,
%         without considering leap days. If daily frequency is specified,
%         the length of the column vector should by a multiple of 365 and
%         for monthly frequency, it should be a multiple of 12.
%         [ndata_obs, 1]
%
%   mod = A column vector of daily, monthly or annual modeled or GCM data,
%         without considering leap days. If daily frequency is specified,
%         the length of the column vector should by a multiple of 365 and
%         for monthly frequency, it should be a multiple of 12.
%         [ndata_mod, 1]
%
%   NOTE: This routine considers that obs and mod series start in the same
%         day/month/year and are continuous until the end day/month/year.
%
% Optional inputs:
%   allow_negatives = A flag that identifies if data allows negative values
%                     and also to replace no-rain values with random small 
%                     values (Chadwick et al., 2023) to avoid numerical
%                     problems with the probability distribution functions.
%                     allow_negatives = 1 or True: Allow negatives (default)
%                     allow_negatives = 0 or False: Do not allow negative
% 
%   frq = A string specifying if the input frequency is daily, monthly or
%         annual.
%         Daily:     frq = 'D'
%         Monthly:   frq = 'M'
%         Annual:    frq = 'A' (default)
% 
%   pp_threshold = A float indicating the threshold to consider no-rain
%                  values. Default is 1.
%
%   pp_factor = A float which multiplied to pp_threshold indicates the
%               maximum value of the random values that replace no-rain
%               values. Default is 1/100.
%
%   day_win = (only for frq='D') An integer indicating how many days to
%              consider backwards and forward to get the statistics
%              of each calendar day.The length of the window will be 
%              (2*win_day-1). For example, day_win=15 -> window of 29.
%              Default: win = 1
%
%   user_pdf = A flag indicating if the user will define the probability
%              distribution functions (pdf) for the observed and modeled
%              series. The distributions will be the same for all periods
%              and sub-periods.
%              user_pdf = 1 or True: User defines the pdf
%              user_pdf = 0 or False: pdf defined by the Kolmogorov
%                                     -Smirnov test (default)
%
%   pdf_obs = An integer indicating the probability distribution function
%             (pdf) to be used for the observed data. The pdf will be the
%             same for all periods and sub-periods. Default: None
%
%   pdf_mod = An integer indicating the probability distribution function
%             (pdf) to be used for the modeled data. The pdf will be the
%             same for all periods and sub-periods. Default: None
%
%   NOTE: The available distributions for pdf_obs and pdf_mod are:
%         0) Normal
%         1) Log-Normal
%         2) Gamma 2 parameters
%         3) Gamma 3 parameters
%            (Pearson 3 parameters)
%         4) Log-Gamma 3 parameters
%            (Log-Pearson 3 parameters)
%         5) Gumbel
%         6) Exponential
%
% Output:
%   QM_series = A column vector of data bias corrected with the QM method.
%               [ndata_mod, 1]
%
% References:
%   Cannon, A. J., S. R. Sobie, and T. Q. Murdock. (2015). Bias correction
%   of GCM precipitation by quantile mapping: How well do methods preserve
%   changes in quantiles and extremes? J. Climate, 28(17), 6938-6959, 
%   https://doi.org/10.1175/JCLI-D-14-00754.1
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

% Check for optional arguments
if ~exist('allow_negatives','var')
    allow_negatives = 1;
end

if ~exist('frq','var')
    frq = 'A';
end

if ~exist('pp_threshold','var')
    pp_threshold = 1;
end

if ~exist('pp_factor','var')
    pp_factor = 1/100;
end

if ~exist('day_win','var')
    day_win = 1;
end

if ~exist('user_pdf','var')
    user_pdf = false;
    pdf_obs = false;
    pdf_mod = false;
end

% Only for daily frequency data that does not allow negatives (as 
% precipitacion). Get a pp_threshold for the modeled series that result
% in the same amount of rain-day values in the historical period as in the
% observed series
if (frq=='D' & allow_negatives==false)
    pp_threshold_mod = get_pp_threshold_mod(obs, mod, pp_threshold);
else
    pp_threshold_mod = pp_threshold;
end

% 1) Format series to matrix with rows as sub-periods and columns as years
% and, if needed, replace no-rain values with random small values 
[y_obs,obs_series] = formatQM(obs, allow_negatives, frq,pp_threshold,pp_factor);
[y_mod,mod_series] = formatQM(mod, allow_negatives, frq,pp_threshold_mod,pp_factor);

% 2) Get statistics of the observed and modeled series in the historical
% period.
if frq=='D'
    % Get a 3D array with a moving window centered on each day
    obs_series_moving = day_centered_moving_window(obs_series, day_win);
    mod_series_moving = day_centered_moving_window(mod_series, day_win);

    % Reshape to 2D in order to have all the days of the window in each row
    obs_series_moving = reshape(permute(obs_series_moving, [1,3,2]),[365,(2*day_win-1)*y_obs]);
    modh_series_moving = reshape(permute(mod_series_moving(:,:,1:y_obs), [1,3,2]),[365,(2*day_win-1)*y_obs]);
    
    % Replace no-rain values with nans but keep at least 30 values to get
    % the statistics. If fewer than 30 values are rain values, the missing
    % ones will be previously random values replaced
    if allow_negatives==false
        obs_series_moving = set_norain_to_nan(obs_series_moving, pp_threshold, pp_factor);
        modh_series_moving = set_norain_to_nan(modh_series_moving, pp_threshold_mod, pp_factor);
        mod_series(mod_series<pp_threshold_mod) = NaN;
    end

    % Get the statistics for each row
    [mu_obs, std_obs, skew_obs, skewy_obs] = getStats(obs_series_moving, frq);
    [mu_mod, std_mod, skew_mod, skewy_mod] = getStats(modh_series_moving, frq);
    
else
    % Get the statistics for each row
    [mu_obs, std_obs, skew_obs, skewy_obs] = getStats(obs_series, frq);
    [mu_mod, std_mod, skew_mod, skewy_mod] = getStats(mod_series(:,1:y_obs), frq);
end


% 3) Assign a probability distribution function to each sub-period for the
%    observed and modeled data in the historical period.
if user_pdf==false
    if frq=='D'
        [pdf_obs, ks_fail_obs] = getDist(obs_series_moving,allow_negatives,mu_obs,std_obs,skew_obs,skewy_obs);
        [pdf_mod, ks_fail_mod] = getDist(modh_series_moving,allow_negatives,mu_mod,std_mod,skew_mod,skewy_mod);
    else
        [pdf_obs, ks_fail_obs] = getDist(obs_series,allow_negatives,mu_obs,std_obs,skew_obs,skewy_obs);
        [pdf_mod, ks_fail_mod] = getDist(mod_series(:,1:y_obs),allow_negatives,mu_mod,std_mod,skew_mod,skewy_mod);
    end
else
    pdf_obs = zeros(size(obs_series,1),1) + pdf_obs;
    pdf_mod = zeros(size(mod_series,1),1) + pdf_mod;
end

% 4) Apply the cumulative distribution function of the modeled data,
%    evaluated with the statistics of the modeled data in the historical
%    period and the modeled data. Equation 1 of Cannon et al. (2015).
prob = getCDF(pdf_mod,mod_series,mu_mod,std_mod,skew_mod,skewy_mod);

% 5) Apply the inverse cumulative distribution function of the observed
%    data, evaluated with the statistics of the observed data in the
%    historical period, to the probabilities obtained from 2). Equation 1
%    of Cannon et al. (2015).
QM_series = getCDFinv(pdf_obs,prob,mu_obs,std_obs,skew_obs,skewy_obs);

% 6) Reshape to a column vector
QM_series=QM_series(:);

% 7) If does not allow negative values, replace nans (should be the result 
%    of correcting removed no-rain values) and replace no-rain values with 0
if allow_negatives==0
    QM_series(isnan(QM_series)) = 0;
    QM_series(QM_series<pp_threshold) = 0;
end

if user_pdf==false
    % Check if ks-test failures
    ks_fail = ks_fail_obs + ks_fail_mod;
    if ks_fail>0
        disp('QM: Some of the probability distribution functions did not pass the KS-Test')
    end
end
end