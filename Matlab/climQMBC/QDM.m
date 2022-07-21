function QDM_series = QDM(obs,mod,var,frq,pp_threshold,pp_factor, rel_change_th,inv_mod_th)
%% QDM_series:
%   This function performs bias correction of modeled series based on
%   observed data by the Quantile Delta Mapping (QDM) method, as described
%   by Cannon et al. (2015) . Correction is performed to monthly or
%   annual precipitation or temperature data in a single location. An
%   independent probability distribution function is assigned to each month
%   and to each projected period based on the Kolmogorov-Smirnov test. If 
%   annual frequency is specified, this is applied to the complete period. 
%   Each projected period considers a window whose length is equal to the 
%   number of years of the historical period and ends in the analyzed year.
%
%   Correction of the historical period is performed by the Quantile 
%   Mapping (QM) method, as described by Cannon et al. (2015). 
%
% Description:
%   0) Check if annual or monthly data is specified.
%
%   1) Format inputs and get statistics of the observed and modeled series 
%      of the historical period.
%
%   2) Assign a probability distribution function to each month for the
%      observed and modeled data in the historical period. If annual
%      frequency is specified, this is applied to the complete historical
%      period.
%
%   3) For each projected period:
%      a) Assign a probability distribution function to each month. If 
%         annual frequency is specified, this is applied to the complete 
%         period.
%      b) Apply the cumulative distribution function of the projected 
%         period, evaluated with the statistics of this period, to the last
%         data of the period.
%
%   4) Apply the inverse cumulative distribution function:
%      a) Of the observed data, evaluated with the statistics of the 
%         observed data in the historical period, to the probabilities 
%         obtained from 3b).
%      b) Of the modeled data, evaluated with the statistics of the 
%         observed data in the historical period, to the probabilities 
%         obtained from 3b).
%
%   5) Get the delta factor or relative change and apply it to the value 
%      obtained in 4b).
%
%   6) Perform QM for the historical period.
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
%   NOTE: This routine considers that obs and mod series start in january
%         of the first year and ends in december of the last year.
%
% Optional inputs:
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
% Output:
%   QDM_series = A column vector of monthly or annual modeled data 
%              (temperature or precipitation) corrected by the QDM method.
%              If monthly frequency is specified, the length of this vector
%              is 12 times the number of observed years [12 x y_mod, 1]. If
%              annual frequency is specified, the length of this vector
%              is equal to the number of observed years [y_mod, 1].
%
% References:
%   Cannon, A. J., S. R. Sobie, and T. Q. Murdock, 2015: Bias correction of
%   GCM precipitation by quantile mapping: How well do methods preserve
%   changes in quantiles and extremes? J. Climate, 28(17), 6938-6959,
%   https://doi.org/10.1175/JCLI-D-14-00754.1
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

% Define optional arguments
if ~exist('pp_threshold','var')
    pp_threshold = 1;
end

if ~exist('pp_factor','var')
    pp_factor = 1/100;
end

if ~exist('rel_change_th','var')
    rel_change_th = 2;
end

if ~exist('inv_mod_th','var')
    inv_mod_th = pp_threshold;
end

% 0) Check if annual or monthly data is specified.
if ~exist('frq','var')
    frq = 'M';
end

% 1) Format inputs and get statistics of the observed and modeled series of
%    the historical period (formatQM function of the climQMBC package).
[y_obs,obs_series,mod_series,mu_obs,std_obs,skew_obs,skewy_obs,mu_mod,std_mod,skew_mod,skewy_mod] = formatQM(obs,mod,var,frq,pp_threshold,pp_factor);

% 2) Assign a probability distribution function to each month of the 
%    historical period (getDist function of the climQMBC package).
PDF_obs = getDist(obs_series,mu_obs,std_obs,skew_obs,skewy_obs,var);
PDF_mod = getDist(mod_series(:,1:y_obs),mu_mod,std_mod,skew_mod,skewy_mod,var);

% 3) For each projected period:
PDF_win = zeros(size(mod_series,1),size(mod_series,2)-y_obs);
Taot = zeros(size(mod_series,1),size(mod_series,2)-y_obs);
for j = 1:size(Taot,2)
    win_series  = mod_series(:,1+j:y_obs+j);
    
    mux  = nanmean(win_series,2);
    sigmax = nanstd(win_series,0,2);
    skewx = skewness(win_series,0,2);
    Ln_win  = log(win_series);
    Ln_win(imag(Ln_win)~=0)= 0;
    skewy = skewness(Ln_win,0,2);
    
    % a) Assign a probability distribution function to each month. If
    %    annual frequency is specified, this is applied to the complete 
    %    period (getDist function of the climQMBC package).
    PDF_win(:,j) = getDist(win_series,mux,sigmax,skewx,skewy,var);
    
    % b) Apply the cumulative distribution function of the projected
    %    period, evaluated with the statistics of this period, to the last
    %    data of the period (getCDF function of the climQMBC package).
    %    Equation 3 of Cannon et al. (2015).
    Taot(:,j) = getCDF(PDF_win(:,j),mod_series(:,y_obs+j),mux,sigmax,skewx,skewy);
end

% 4) Apply the inverse cumulative distribution function:
%    a) Of the observed data, evaluated with the statistics of the observed
%       data in the historical period, to the probabilities obtained from 
%       3b) (getCDFinv function of the climQMBC package). Equation 5 of 
%       Cannon et al. (2015).
inv_obs = getCDFinv(PDF_obs,Taot,mu_obs,std_obs,skew_obs,skewy_obs);

%    b) Of the modeled data, evaluated with the statistics of the observed
%       data in the historical period, to the probabilities obtained from 
%       3b) (getCDFinv function of the climQMBC package). Equations 4 of 
%       Cannon et al. (2015).
inv_mod = getCDFinv(PDF_mod,Taot,mu_mod,std_mod,skew_mod,skewy_mod);

% 5) Get the delta factor or relative change and apply it to the value
%    obtained in 4b). Equation 4 and 6 of Cannon et al. (2015).
if var == 1
    DM = mod_series(:,y_obs+1:end)./inv_mod;
    DM((DM>rel_change_th) & (inv_mod<inv_mod_th)) = rel_change_th;
    QDM = inv_obs.*DM;
else
    DM = mod_series(:,y_obs+1:end) - inv_mod;
    QDM = inv_obs + DM;
end

QDM = QDM(:);

% 6) Perform QM for the historical period.
mod_h = mod_series(:,1:y_obs);
mod_h = mod_h(:);
QM_series = QM(obs,mod_h,var,frq);
QDM_series = [QM_series' QDM']';
QDM_series(QDM_series<pp_threshold) = 0;
end