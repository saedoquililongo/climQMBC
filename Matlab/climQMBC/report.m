function [QM_series,DQM_series,QDM_series,UQM_series,SDM_series] = report(obs,mod,SDM_var,mult_change,allow_negatives,fun,y_init,y_wind)
%% report:
%   This function generates two report of the performance of the different
%   methods (QM, DQM, QDM, UQM and SDM) available in the climQMBC package.
%   These reports are based on the mean and standard deviation of the
%   series in the historical and future projected periods.
%
%   The first report is a summary table with the overall performance of the
%   QM method in the historical period and of the different methods in
%   future projected periods. The methods performance is addressed by
%   comparing its variations in the mean and standard deviation with the
%   ones of the modeled data.
%
%   The second report consist of three figures. Figure 1 shows the
%   cumulative distribution of the observed, modeled, and corrected series
%   in the historical and future period. This figure also shows the
%   complete time series. Figure 2 and 3 shows the monthly mean and
%   standard deviation, respectively, of each series in the historical and
%   future projected periods. In these two figures, in the future projected
%   periods, the observed series is replaced by an objective series which
%   is computed as the observed series (monthly mean or standard
%   deviation), scaled by the variation between the projected and the
%   historical period of the modeled series. Each projected period is
%   centered in a moving window whose length is equal to the length of the
%   historical period.
%
%   NOTE: 
%      1) Even if a set of methods are specified in the optional inputs,
%         this function returns all five methods implemented in the
%         climQMBC package.
%      2) This report function was built for monthly data.
%
% Description:
%   0) Get the number of observed and modeled years.
%
%   1) Set non-declared arguments.
%
%   2) Apply bias correction methods.
%
%   3) Get observed, modeled and bias corrected statistics of the
%      historical and complete future period.
%
%   4) Get observed, modeled and bias corrected statistics of the projected
%      periods.
%
%   5) Get delta mean and delta standard deviation.
%
%   6) Display report:
%          a) Report description.
%          b) Table.
%          c) Figures.
%
% Input:
%   obs = A column vector of monthly observed data (temperature or
%         precipitation) [12 x y_obs, 1].
%
%   mod = A column vector of monthly modeled data (temperature or
%         precipitation) [12 x y_mod, 1].
%
%   var = A flag that identifies if data are temperature or precipitation.
%         This flag tells the getDist function if it has to discard
%         distribution functions that allow negative numbers, and if the 
%         terms in the correction equations are multiplied/divided or
%         added/subtracted.
%         Temperature:   var = 0
%         Precipitation: var = 1
%
% Optional inputs:
%   fun = A cell array of strings with the desired bias correction methods
%         to be reported. If this input is not recieved by the function,
%         all bias correction methods available in the climQMBC package
%         will be reported. The methods supported are:
%          a) 'QM' : Quantile Mapping
%          b) 'DQM': Detrended Quantile Mapping
%          c) 'QDM': Quantile Delta Mapping
%          d) 'UQM': Unbiased Quantile Mapping
%          e) 'SDM': Scaled Distribution Mapping
%         Default: fun = {'QM','DQM','QDM','UQM','SDM'}
%
%   y_init = First year of the observed and modeled series (integer).
%            Default: y_init = 0
%
%   y_wind = A row array of integers with the year of the center of the
%            projected periods to be reported.
%            Default: y_wind = [floor(y_obs+y_obs/2) floor(y_mod-y_obs/2)]
%                     This value sets a first projected period just after
%                     the end of historical period, and a second projected 
%                     period just before the end of the modeled series.
%
%   NOTE: As MATLAB reads the inputs in sequence, each optional inputs must
%         be added as defined by the function. If y_init or y_wind is a
%         desired argument, fun must be specified.
%
% Output:
%   QM_series = A column vector of monthly modeled data (temperature or 
%               precipitation) corrected by the QM method [12 x y_mod, 1].
%
%   DQM_series = A column vector of monthly modeled data (temperature or 
%                precipitation) corrected by the DQM method 
%                [12 x y_mod, 1].
%
%   QDM_series = A column vector of monthly modeled data (temperature or 
%                precipitation) corrected by the QdM method 
%                [12 x y_mod, 1].
%
%   UQM_series = A column vector of monthly modeled data (temperature or 
%                precipitation) corrected by the UQM method 
%                [12 x y_mod, 1].
%
%   SDM_series = A column vector of monthly modeled data (temperature or 
%                precipitation) corrected by the SDM method 
%                [12 x y_mod, 1].
%
%   NOTE: This function returns all five bias correction methods, 
%         independently of which methods are specified for this report.
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
% Revision: 0, updated Dec 2021
% 

%%
% 0) Get the number of observed and modeled years
n_obs = length(obs);
y_obs = n_obs/12;

n_mod = length(mod);
y_mod = n_mod/12;

% 1) Set non-declared arguments
if ~exist('mult_change','var')
  mult_change=1;
end

if ~exist('allow_negatives','var')
  allow_negatives=1;
end

if ~exist('fun','var')
  fun = {'QM','DQM','QDM','UQM','SDM'};
end
if ~exist('y_init','var')
  y_init = 0;
end
if ~exist('y_wind','var')
  y_wind = [floor(y_obs+y_obs/2) floor(y_mod-y_obs/2)];
  w_label = {'First period','Last period'};
  rep_head = ['Variation           Hist.     First     Last'];
else
    rep_head = ['Variation           Hist'];
    w_label = {};
    for w = 1:length(y_wind)
        rep_head = [rep_head '      ' num2str(y_wind(w))];
        w_label{end+1} = num2str(y_wind(w));
    end
end

y_wind = y_wind - y_init;

% 2) Apply QM methods
QM_series = QM(obs,mod,allow_negatives,'M');
DQM_series = DQM(obs,mod,mult_change,allow_negatives,'M');
QDM_series = QDM(obs,mod,mult_change,allow_negatives,'M');
UQM_series = UQM(obs,mod,mult_change,allow_negatives,'M');
SDM_series = SDM(obs,mod,SDM_var,'M');

% 3) Get observed, modeled and bias corrected statistics
% a) Get obs, mod, and QMs as [12 x n] matrix
obs_M = reshape(obs,12,[]);
mod_M = reshape(mod,12,[]);
QM_M = reshape(QM_series,12,[]);
DQM_M = reshape(DQM_series,12,[]);
QDM_M = reshape(QDM_series,12,[]);
UQM_M = reshape(UQM_series,12,[]);
SDM_M = reshape(SDM_series,12,[]);

% b) Mean
mu_obs = mean(mean(obs_M));
mu_mod_h = mean(mean(mod_M(:,1:y_obs)));
mu_mod_f = mean(mean(mod_M(:,y_obs+1:y_mod)));
mu_QM_h = mean(mean(QM_M(:,1:y_obs)));
mu_QM_f = mean(mean(QM_M(:,y_obs+1:y_mod)));
mu_DQM_f = mean(mean(DQM_M(:,y_obs+1:y_mod)));
mu_QDM_f = mean(mean(QDM_M(:,y_obs+1:y_mod)));
mu_UQM_f = mean(mean(UQM_M(:,y_obs+1:y_mod)));
mu_SDM_f = mean(mean(SDM_M(:,y_obs+1:y_mod)));

% c) Monthly mean
mu_obs_M = mean(obs_M,2);
mu_mod_Mh = mean(mod_M(:,1:y_obs),2);
mu_mod_Mf = mean(mod_M(:,y_obs+1:y_mod),2);
mu_QM_Mh = mean(QM_M(:,1:y_obs),2);
mu_QM_Mf = mean(QM_M(:,y_obs+1:y_mod),2);
mu_DQM_Mf = mean(DQM_M(:,y_obs+1:y_mod),2);
mu_QDM_Mf = mean(QDM_M(:,y_obs+1:y_mod),2);
mu_UQM_Mf = mean(UQM_M(:,y_obs+1:y_mod),2);
mu_SDM_Mf = mean(SDM_M(:,y_obs+1:y_mod),2);

% d) Standard deviation
s_obs = std(reshape(obs_M,[],1),0);
s_mod_h = std(reshape(mod_M(:,1:y_obs),[],1),0);
s_mod_f = std(reshape(mod_M(:,y_obs+1:y_mod),[],1),0);
s_QM_h = std(reshape(QM_M(:,1:y_obs),[],1),0);
s_QM_f = std(reshape(QM_M(:,y_obs+1:y_mod),[],1),0);
s_DQM_f = std(reshape(DQM_M(:,y_obs+1:y_mod),[],1),0);
s_QDM_f = std(reshape(QDM_M(:,y_obs+1:y_mod),[],1),0);
s_UQM_f = std(reshape(UQM_M(:,y_obs+1:y_mod),[],1),0);
s_SDM_f = std(reshape(SDM_M(:,y_obs+1:y_mod),[],1),0);

% e) Monthly standard deviation
s_obs_M = std(obs_M,0,2);
s_mod_Mh = std(mod_M(:,1:y_obs),0,2);
s_mod_Mf = std(mod_M(:,y_obs+1:y_mod),0,2);
s_QM_Mh = std(QM_M(:,1:y_obs),0,2);
s_QM_Mf = std(QM_M(:,y_obs+1:y_mod),0,2);
s_DQM_Mf = std(DQM_M(:,y_obs+1:y_mod),0,2);
s_QDM_Mf = std(QDM_M(:,y_obs+1:y_mod),0,2);
s_UQM_Mf = std(UQM_M(:,y_obs+1:y_mod),0,2);
s_SDM_Mf = std(SDM_M(:,y_obs+1:y_mod),0,2);

% 4) Get observed, modeled and bias corrected statistics of the projected
%    periods
% a) Projected period mean
mu_mod_w = zeros(length(y_wind),1);
mu_QM_w = zeros(length(y_wind),1);
mu_DQM_w = zeros(length(y_wind),1);
mu_QDM_w = zeros(length(y_wind),1);
mu_UQM_w = zeros(length(y_wind),1);
mu_SDM_w = zeros(length(y_wind),1);

% b) Projected period monthly mean
mu_mod_Mw = zeros(12,length(y_wind));
mu_QM_Mw = zeros(12,length(y_wind));
mu_DQM_Mw = zeros(12,length(y_wind));
mu_QDM_Mw = zeros(12,length(y_wind));
mu_UQM_Mw = zeros(12,length(y_wind));
mu_SDM_Mw = zeros(12,length(y_wind));

% c) Projected period standard deviation
s_mod_w = zeros(length(y_wind),1);
s_QM_w = zeros(length(y_wind),1);
s_DQM_w = zeros(length(y_wind),1);
s_QDM_w = zeros(length(y_wind),1);
s_UQM_w = zeros(length(y_wind),1);
s_SDM_w = zeros(length(y_wind),1);

% d) Projected period monthly standard deviation
s_mod_Mw = zeros(12,length(y_wind));
s_QM_Mw = zeros(12,length(y_wind));
s_DQM_Mw = zeros(12,length(y_wind));
s_QDM_Mw = zeros(12,length(y_wind));
s_UQM_Mw = zeros(12,length(y_wind));
s_SDM_Mw = zeros(12,length(y_wind));

% Compute projected period statistics
for w = 1:length(y_wind)
    % Projected period mean
    mu_mod_w(w) = mean(mean(mod_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2))));
    mu_QM_w(w) = mean(mean(QM_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2))));
    mu_DQM_w(w) = mean(mean(DQM_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2))));
    mu_QDM_w(w) = mean(mean(QDM_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2))));
    mu_UQM_w(w) = mean(mean(UQM_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2))));
    mu_SDM_w(w) = mean(mean(SDM_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2))));
    
    % Projected period monthly mean
    mu_mod_Mw(:,w) = mean(mod_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2)),2);
    mu_QM_Mw(:,w) = mean(QM_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2)),2);
    mu_DQM_Mw(:,w) = mean(DQM_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2)),2);
    mu_QDM_Mw(:,w) = mean(QDM_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2)),2);
    mu_UQM_Mw(:,w) = mean(UQM_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2)),2);
    mu_SDM_Mw(:,w) = mean(SDM_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2)),2);
    
    % Projected period standard deviation
    s_mod_w(w) = std(reshape(mod_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2)),[],1),0);
    s_QM_w(w) = std(reshape(QM_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2)),[],1),0);
    s_DQM_w(w) = std(reshape(DQM_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2)),[],1),0);
    s_QDM_w(w) = std(reshape(QDM_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2)),[],1),0);
    s_UQM_w(w) = std(reshape(UQM_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2)),[],1),0);
    s_SDM_w(w) = std(reshape(SDM_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2)),[],1),0);
    
    % Projected period monthly standard deviation
    s_mod_Mw(:,w) = std(mod_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2)),0,2);
    s_QM_Mw(:,w) = std(QM_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2)),0,2);
    s_DQM_Mw(:,w) = std(DQM_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2)),0,2);
    s_QDM_Mw(:,w) = std(QDM_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2)),0,2);
    s_UQM_Mw(:,w) = std(UQM_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2)),0,2);
    s_SDM_Mw(:,w) = std(SDM_M(:,floor(y_wind(w)-y_obs/2):floor(y_wind(w)+y_obs/2)),0,2);
end

% 5) Get delta mean and delta standard deviation
dm_mod_w = zeros(length(y_wind),1);
dm_QM_w = zeros(length(y_wind),1);
dm_DQM_w = zeros(length(y_wind),1);
dm_QDM_w = zeros(length(y_wind),1);
dm_UQM_w = zeros(length(y_wind),1);
dm_SDM_w = zeros(length(y_wind),1);

dm_mod_Mw = zeros(12,length(y_wind));
dm_QM_Mw = zeros(12,length(y_wind));
dm_DQM_Mw = zeros(12,length(y_wind));
dm_QDM_Mw = zeros(12,length(y_wind));
dm_UQM_Mw = zeros(12,length(y_wind));
dm_SDM_Mw = zeros(12,length(y_wind));

ds_mod_w = zeros(length(y_wind),1);
ds_QM_w = zeros(length(y_wind),1);
ds_DQM_w = zeros(length(y_wind),1);
ds_QDM_w = zeros(length(y_wind),1);
ds_UQM_w = zeros(length(y_wind),1);
ds_SDM_w = zeros(length(y_wind),1);

ds_mod_Mw = zeros(12,length(y_wind));
ds_QM_Mw = zeros(12,length(y_wind));
ds_DQM_Mw = zeros(12,length(y_wind));
ds_QDM_Mw = zeros(12,length(y_wind));
ds_UQM_Mw = zeros(12,length(y_wind));
ds_SDM_Mw = zeros(12,length(y_wind));

if mult_change == 1
    dm_obs = mu_QM_h/mu_obs;
    dm_mod = mu_mod_f/mu_mod_h;
    dm_QM = mu_QM_f/mu_QM_h;
    dm_DQM = mu_DQM_f/mu_QM_h;
    dm_QDM = mu_QDM_f/mu_QM_h;
    dm_UQM = mu_UQM_f/mu_QM_h;
    dm_SDM = mu_SDM_f/mu_QM_h;
    
    ds_obs = s_QM_h/s_obs;
    ds_mod = s_mod_f/s_mod_h;
    ds_QM = s_QM_f/s_QM_h;
    ds_DQM = s_DQM_f/s_QM_h;
    ds_QDM = s_QDM_f/s_QM_h;
    ds_UQM = s_UQM_f/s_QM_h;
    ds_SDM = s_SDM_f/s_QM_h;
    
    % Scale factor for the objective series in the future projected series.
    dm_mod_M = mu_mod_Mf./mu_mod_Mh;
    ds_mod_M = s_mod_Mf./s_mod_Mh;
    
    for w = 1:length(y_wind)
        dm_mod_w(w) = mu_mod_w(w)/mu_mod_h;
        dm_QM_w(w) = mu_QM_w(w)/mu_QM_h;
        dm_DQM_w(w) = mu_DQM_w(w)/mu_QM_h;
        dm_QDM_w(w) = mu_QDM_w(w)/mu_QM_h;
        dm_UQM_w(w) = mu_UQM_w(w)/mu_QM_h;
        dm_SDM_w(w) = mu_SDM_w(w)/mu_QM_h;
        
        ds_mod_w(w) = s_mod_w(w)/s_mod_h;
        ds_QM_w(w) = s_QM_w(w)/s_QM_h;
        ds_DQM_w(w) = s_DQM_w(w)/s_QM_h;
        ds_QDM_w(w) = s_QDM_w(w)/s_QM_h;
        ds_UQM_w(w) = s_UQM_w(w)/s_QM_h;
        ds_SDM_w(w) = s_SDM_w(w)/s_QM_h;
        
        dm_mod_Mw(:,w) = mu_mod_Mw(:,w)./mu_mod_Mh;
        dm_QM_Mw(:,w) = mu_QM_Mw(:,w)./mu_QM_Mh;
        dm_DQM_Mw(:,w) = mu_DQM_Mw(:,w)./mu_QM_Mh;
        dm_QDM_Mw(:,w) = mu_QDM_Mw(:,w)./mu_QM_Mh;
        dm_UQM_Mw(:,w) = mu_UQM_Mw(:,w)./mu_QM_Mh;
        dm_SDM_Mw(:,w) = mu_SDM_Mw(:,w)./mu_QM_Mh;
        
        ds_mod_Mw(:,w) = s_mod_Mw(:,w)./s_mod_Mh;
        ds_QM_Mw(:,w) = s_QM_Mw(:,w)./s_QM_Mh;
        ds_DQM_Mw(:,w) = s_DQM_Mw(:,w)./s_QM_Mh;
        ds_QDM_Mw(:,w) = s_QDM_Mw(:,w)./s_QM_Mh;
        ds_UQM_Mw(:,w) = s_UQM_Mw(:,w)./s_QM_Mh;
        ds_SDM_Mw(:,w) = s_SDM_Mw(:,w)./s_QM_Mh;
    end
else
    dm_obs = mu_QM_h-mu_obs;
    dm_mod = mu_mod_f-mu_mod_h;
    dm_QM = mu_QM_f-mu_QM_h;
    dm_DQM = mu_DQM_f-mu_QM_h;
    dm_QDM = mu_QDM_f-mu_QM_h;
    dm_UQM = mu_UQM_f-mu_QM_h;
    dm_SDM = mu_SDM_f-mu_QM_h;
    
    ds_obs = s_QM_h-s_obs;
    ds_mod = s_mod_f-s_mod_h;
    ds_QM = s_QM_f-s_QM_h;
    ds_DQM = s_DQM_f-s_QM_h;
    ds_QDM = s_QDM_f-s_QM_h;
    ds_UQM = s_UQM_f-s_QM_h;
    ds_SDM = s_SDM_f-s_QM_h;
    
    % Scale factor for the objective series in the future projected series.
    dm_mod_M = mu_mod_Mf - mu_mod_Mh;
    ds_mod_M = s_mod_Mf - s_mod_Mh;
    
    for w = 1:length(y_wind)
        dm_mod_w(w) = mu_mod_w(w)-mu_mod_h;
        dm_QM_w(w) = mu_QM_w(w)-mu_QM_h;
        dm_DQM_w(w) = mu_DQM_w(w)-mu_QM_h;
        dm_QDM_w(w) = mu_QDM_w(w)-mu_QM_h;
        dm_UQM_w(w) = mu_UQM_w(w)-mu_QM_h;
        dm_SDM_w(w) = mu_SDM_w(w)-mu_QM_h;
        
        ds_mod_w(w) = s_mod_w(w)-s_mod_h;
        ds_QM_w(w) = s_QM_w(w)-s_QM_h;
        ds_DQM_w(w) = s_DQM_w(w)-s_QM_h;
        ds_QDM_w(w) = s_QDM_w(w)-s_QM_h;
        ds_UQM_w(w) = s_UQM_w(w)-s_QM_h;
        ds_SDM_w(w) = s_SDM_w(w)-s_QM_h;
        
        dm_mod_Mw(:,w) = mu_mod_Mw(:,w)-mu_mod_Mh;
        dm_QM_Mw(:,w) = mu_QM_Mw(:,w)-mu_QM_Mh;
        dm_DQM_Mw(:,w) = mu_DQM_Mw(:,w)-mu_QM_Mh;
        dm_QDM_Mw(:,w) = mu_QDM_Mw(:,w)-mu_QM_Mh;
        dm_UQM_Mw(:,w) = mu_UQM_Mw(:,w)-mu_QM_Mh;
        dm_SDM_Mw(:,w) = mu_SDM_Mw(:,w)-mu_QM_Mh;
        
        ds_mod_Mw(:,w) = s_mod_Mw(:,w)-s_mod_Mh;
        ds_QM_Mw(:,w) = s_QM_Mw(:,w)-s_QM_Mh;
        ds_DQM_Mw(:,w) = s_DQM_Mw(:,w)-s_QM_Mh;
        ds_QDM_Mw(:,w) = s_QDM_Mw(:,w)-s_QM_Mh;
        ds_UQM_Mw(:,w) = s_UQM_Mw(:,w)-s_QM_Mh;
        ds_SDM_Mw(:,w) = s_SDM_Mw(:,w)-s_QM_Mh;
    end
end

% 6) Display report
% a) Report description
disp('<strong>Description of this report:</strong>')
disp(' ')
disp('Description of the table')
disp ('The first line shows the difference (relative or absolute, according')
disp('to the variable analyzed) between the QM method and the observed data')
disp('in the historical period. For precipitation, a value of 1 means that')
disp('the mean precipitation of the observed data and QM series are the same')
disp('for the historical period. For temperature, this is achieved with a')
disp('value of 0. The second line shows the headers of the periods reported.')
disp('Remember that the moving window is centered in the period and its')
disp('length is equal to the length of the historical period. The first ')
disp('column (Hist.) shows the difference between the future and the')
disp('historical period. The next columns show the difference between the')
disp('corresponding period and the historical period. The third and the')
disp('following lines show the difference in each period for the')
disp('corresponding series. As a general rule, it is desirable that the')
disp('difference between the future periods and the historical period of the')
disp('bias corrected series match the differences of the modeled series.')
disp(' ')
disp('Description of the figures')
disp('Figure 1: (left) shows the cumulative distribution functions of each')
disp('series split between historical and future period. (right) shows the')
disp('observed, modeled, and corrected series.')
disp('Figure 2: shows the monthly mean of the historical and future period.')
disp('Additionally, monthly mean of the projected periods is displayed.')
disp('Figure 3: shows the monthly standard deviation of the historical and')
disp('future period. Additionally, monthly mean of the projected periods is')
disp('displayed.')
disp('* In Fig. 2 and 3, the obj (red) line represents objective values.')
disp('This objective series is computed as the observed series (monthly mean')
disp('or standard deviation), scaled by the variation between the projected')
disp('and the historical period of the modeled series. Each projected period')
disp('is centered in  a moving window whose length is equal to the length of')
disp('the historical period.')
disp(' ')

% Format each line of the table report
n = 3; % Number of decimals to which round the values
sp = ['-',' ',' ']; % A small trick to prevent negative values to use more
                    % space than positive values

rep_mu_M = ['Modeled          :' ' ' sp(sign(dm_mod)+2) num2str(abs(round(dm_mod,n)))];
rep_mu_QM = ['QM               :' ' ' sp(sign(dm_QM)+2) num2str(abs(round(dm_QM,n)))];
rep_mu_DQM = ['DQM              :' ' ' sp(sign(dm_DQM)+2) num2str(abs(round(dm_DQM,n)))];
rep_mu_QDM = ['QDM              :' ' ' sp(sign(dm_QDM)+2) num2str(abs(round(dm_QDM,n)))];
rep_mu_UQM = ['UQM              :' ' ' sp(sign(dm_UQM)+2) num2str(abs(round(dm_UQM,n)))];
rep_mu_SDM = ['SDM              :' ' ' sp(sign(dm_SDM)+2) num2str(abs(round(dm_SDM,n)))];

rep_s_M = ['Modeled          :' ' ' sp(sign(ds_mod)+2) num2str(abs(round(ds_mod,n)))];
rep_s_QM = ['QM               :' ' ' sp(sign(ds_QM)+2) num2str(abs(round(ds_QM,n)))];
rep_s_DQM = ['DQM              :' ' ' sp(sign(ds_DQM)+2) num2str(abs(round(ds_DQM,n)))];
rep_s_QDM = ['QDM              :' ' ' sp(sign(ds_QDM)+2) num2str(abs(round(ds_QDM,n)))];
rep_s_UQM = ['UQM              :' ' ' sp(sign(ds_UQM)+2) num2str(abs(round(ds_UQM,n)))];
rep_s_SDM = ['SDM              :' ' ' sp(sign(ds_SDM)+2) num2str(abs(round(ds_SDM,n)))];

for w = 1:length(y_wind)
    rep_mu_M = [rep_mu_M '  ; ' sp(sign(dm_mod_w(w))+2) num2str(abs(round(dm_mod_w(w),n)))];
    rep_mu_QM = [rep_mu_QM '  ; ' sp(sign(dm_QM_w(w))+2) num2str(abs(round(dm_QM_w(w),n)))];
    rep_mu_DQM = [rep_mu_DQM '  ; ' sp(sign(dm_DQM_w(w))+2) num2str(abs(round(dm_DQM_w(w),n)))];
    rep_mu_QDM = [rep_mu_QDM '  ; ' sp(sign(dm_QDM_w(w))+2) num2str(abs(round(dm_QDM_w(w),n)))];
    rep_mu_UQM = [rep_mu_UQM '  ; ' sp(sign(dm_UQM_w(w))+2) num2str(abs(round(dm_UQM_w(w),n)))];
    rep_mu_SDM = [rep_mu_SDM '  ; ' sp(sign(dm_SDM_w(w))+2) num2str(abs(round(dm_SDM_w(w),n)))];
    
    rep_s_M = [rep_s_M '  ; ' sp(sign(ds_mod_w(w))+2) num2str(abs(round(ds_mod_w(w),n)))];
    rep_s_QM = [rep_s_QM '  ; ' sp(sign(ds_QM_w(w))+2) num2str(abs(round(ds_QM_w(w),n)))];
    rep_s_DQM = [rep_s_DQM '  ; ' sp(sign(ds_DQM_w(w))+2) num2str(abs(round(ds_DQM_w(w),n)))];
    rep_s_QDM = [rep_s_QDM '  ; ' sp(sign(ds_QDM_w(w))+2) num2str(abs(round(ds_QDM_w(w),n)))];
    rep_s_UQM = [rep_s_UQM '  ; ' sp(sign(ds_UQM_w(w))+2) num2str(abs(round(ds_UQM_w(w),n)))];
    rep_s_SDM = [rep_s_SDM '  ; ' sp(sign(ds_SDM_w(w))+2) num2str(abs(round(ds_SDM_w(w),n)))];
end

% b) Table
disp('<strong>Mean</strong>')
disp(['QM performance in historical period  :' ' ' sp(sign(dm_obs)+2) num2str(abs(round(dm_obs,n)))])
disp(rep_head)
disp(rep_mu_M)
if any(strcmp(fun,'QM'))
    disp(rep_mu_QM)
end
if any(strcmp(fun,'DQM'))
    disp(rep_mu_DQM)
end
if any(strcmp(fun,'QDM'))
    disp(rep_mu_QDM)
end
if any(strcmp(fun,'UQM'))
    disp(rep_mu_UQM)
end
if any(strcmp(fun,'SDM'))
    disp(rep_mu_SDM)
end
disp(' ')

disp('<strong>Standard deviation</strong>')
disp(['QM performance in historical period  :' ' ' sp(sign(ds_obs)+2) num2str(abs(round(ds_obs,n)))])
disp(rep_head)
disp(rep_s_M)
if any(strcmp(fun,'QM'))
    disp(rep_s_QM)
end
if any(strcmp(fun,'DQM'))
    disp(rep_s_DQM)
end
if any(strcmp(fun,'QDM'))
    disp(rep_s_QDM)
end
if any(strcmp(fun,'UQM'))
    disp(rep_s_UQM)
end
if any(strcmp(fun,'SDM'))
    disp(rep_s_SDM)
end
disp(' ')

% c) Figures
% Figure 1: Cumulative probabilities and time series
% Get empirical cumulative distribution functions
[fo,xo] = ecdf(obs);
[fm,xm] = ecdf(mod(1:length(obs)));
[fm_,xm_] = ecdf(mod(length(obs)+1:end));
[fq,xq] = ecdf(QM_series(1:length(obs)));
[fq1,xq1] = ecdf(QM_series(length(obs)+1:end));
[fq2,xq2] = ecdf(DQM_series(length(obs)+1:end));
[fq3,xq3] = ecdf(QDM_series(length(obs)+1:end));
[fq4,xq4] = ecdf(UQM_series(length(obs)+1:end));
[fq5,xq5] = ecdf(SDM_series(length(obs)+1:end));

lgnd = {'Obs_h','Mod_h','Mod_f','QM_h'};
% Plot empirical cumulative distribution function
figure();
subplot(1,2,1)
plot(xo,fo,'r')
grid on
hold on
plot(xm,fm,'b')
plot(xm_,fm_,'b--')
plot(xq,fq,'k')
if any(strcmp(fun,'QM'))
    plot(xq1,fq1,'k--')
    lgnd{end+1} = 'QM_f';
end
if any(strcmp(fun,'DQM'))
    plot(xq2,fq2,'g--')
    lgnd{end+1} = 'DQM_f';
end
if any(strcmp(fun,'QDM'))
    plot(xq3,fq3,'y--')
    lgnd{end+1} = 'QDM_f';
end
if any(strcmp(fun,'UQM'))
    plot(xq4,fq4,'cyan--')
    lgnd{end+1} = 'UQM_f';
end
if any(strcmp(fun,'SDM'))
    plot(xq5,fq5,'m--')
    lgnd{end+1} = 'SDM_f';
end
hold off
title('Empirical cumulative ditribution functions')
if mult_change ==1
    xlabel('Precipitation (mm)') 
else
    xlabel('Temperature (C)') 
end
ylabel('Probability') 
legend(lgnd,'location','southeast')


lgnd = {'Mod'};
% Plot series
subplot(1,2,2)
plot(mod,'b')
grid on
hold on
if any(strcmp(fun,'QM'))
    plot(QM_series,'k')
    lgnd{end+1} = 'QM';
end
if any(strcmp(fun,'DQM'))
    plot(DQM_series,'g')
    lgnd{end+1} = 'DQM';
end
if any(strcmp(fun,'QDM'))
    plot(QDM_series,'y')
    lgnd{end+1} = 'QDM';
end
if any(strcmp(fun,'UQM'))
    plot(UQM_series,'cyan')
    lgnd{end+1} = 'UQM';
end
if any(strcmp(fun,'SDM'))
    plot(SDM_series,'m')
    lgnd{end+1} = 'SDM';
end
plot(obs,'r')
lgnd{end+1} = 'Obs';
hold off
if mult_change ==1
    title('Precipitation time series')
    ylabel('Precipitation (mm)') 
else
    title('Temperature time series')
    ylabel('Temperature (C)')
end
xlabel('Months since starting date')
legend(lgnd,'orientation','horizontal','location','North')
figure;

% Figure 2: Monthly means
% Historical period
subplot(1,length(y_wind)+2,1)
plot(mu_mod_Mh,'b')
grid on
hold on
plot(mu_QM_Mh,'k')
hold on
plot(mu_obs_M,'r')
hold off
title('Monthly mean of the historical period')
xlabel('Month') 
legend({'Mod_h','QM_h','Obs'},'orientation','horizontal','location','North')
xlim([1,12])

% Future period
lgnd = {'Mod'};
subplot(1,length(y_wind)+2,2)
plot(mu_mod_Mf,'b')
grid on
hold on
if any(strcmp(fun,'QM'))
    plot(mu_QM_Mf,'k')
    lgnd{end+1} = 'QM';
end
if any(strcmp(fun,'DQM'))
    plot(mu_DQM_Mf,'g')
    lgnd{end+1} = 'DQM';
end
if any(strcmp(fun,'QDM'))
    plot(mu_QDM_Mf,'y')
    lgnd{end+1} = 'QDM';
end
if any(strcmp(fun,'UQM'))
    plot(mu_UQM_Mf,'cyan')
    lgnd{end+1} = 'UQM';
end
if any(strcmp(fun,'SDM'))
    plot(mu_SDM_Mf,'m')
    lgnd{end+1} = 'SDM';
end
if mult_change==1
    plot(mu_obs_M.*dm_mod_M ,'r')
elseif mult_change ==0
    plot(mu_obs_M+dm_mod_M ,'r')
end
lgnd{end+1} = 'Obj';
hold off
xlim([1,12])
title('Monthly mean of the future period')
xlabel('Month')

% Projected periods
for w=1:length(y_wind)
    subplot(1,length(y_wind)+2,2+w)
    plot(mu_mod_Mw(:,w),'b')
    grid on
    hold on 
    if any(strcmp(fun,'QM'))
        plot(mu_QM_Mw(:,w),'k')
    end
    if any(strcmp(fun,'DQM'))
        plot(mu_DQM_Mw(:,w),'g')
    end
    if any(strcmp(fun,'QDM'))
        plot(mu_QDM_Mw(:,w),'y')
    end
    if any(strcmp(fun,'UQM'))
        plot(mu_UQM_Mw(:,w),'cyan')
    end
    if any(strcmp(fun,'SDM'))
        plot(mu_SDM_Mw(:,w),'m')
    end
    if mult_change==1
        plot(mu_obs_M.*dm_mod_Mw(:,w) ,'r')
    elseif mult_change ==0
        plot(mu_obs_M+dm_mod_Mw(:,w) ,'r')
    end
    hold off
    title(w_label(w))
    xlabel('Months')
    xlim([1,12])
end
legend(lgnd,'orientation','horizontal','location','North')
figure;

% Figure 3: Monthly standard deviation
% Historical period
lgnd = {'Mod'};
subplot(1,length(y_wind)+2,1)
plot(s_mod_Mh,'b')
grid on
hold on
plot(s_QM_Mh,'k')
hold on
plot(s_obs_M,'r')
hold off
title('Monthly std. dev. of the historical period')
xlabel('Months') 
legend({'Mod_h','QM_h','Obs'},'orientation','horizontal','location','North')
xlim([1,12])

% Future period
subplot(1,length(y_wind)+2,2)
plot(s_mod_Mf,'b')
grid on
hold on
if any(strcmp(fun,'QM'))
    plot(s_QM_Mf,'k')
    lgnd{end+1} = 'QM';
end
if any(strcmp(fun,'DQM'))
    plot(s_DQM_Mf,'g')
    lgnd{end+1} = 'DQM';
end
if any(strcmp(fun,'QDM'))
    plot(s_QDM_Mf,'y')
    lgnd{end+1} = 'QDM';
end
if any(strcmp(fun,'UQM'))
    plot(s_UQM_Mf,'cyan')
    lgnd{end+1} = 'UQM';
end
if any(strcmp(fun,'SDM'))
    plot(s_SDM_Mf,'m')
    lgnd{end+1} = 'SDM';
end
if mult_change==1
    plot(s_obs_M.*ds_mod_M ,'r')
elseif mult_change ==0
    plot(s_obs_M+ds_mod_M ,'r')
end
lgnd{end+1} = 'Obj';
hold off
title('Monthly std. dev. of the future period')
xlabel('Months')
xlim([1,12])

% Projected periods
for w=1:length(y_wind)
    subplot(1,length(y_wind)+2,2+w)
    plot(s_mod_Mw(:,w),'b')
    grid on
    hold on
    if any(strcmp(fun,'QM'))
        plot(s_QM_Mw(:,w),'k')
    end
    if any(strcmp(fun,'DQM'))
        plot(s_DQM_Mw(:,w),'g')
    end
    if any(strcmp(fun,'QDM'))
        plot(s_QDM_Mw(:,w),'y')
    end
    if any(strcmp(fun,'UQM'))
        plot(s_UQM_Mw(:,w),'cyan')
    end
    if any(strcmp(fun,'SDM'))
        plot(s_SDM_Mw(:,w),'m')
    end
    if mult_change==1
        plot(s_obs_M.*ds_mod_Mw(:,w) ,'r')
    elseif mult_change ==0
        plot(s_obs_M+ds_mod_Mw(:,w) ,'r')
    end
    hold off
    title(w_label(w))
    xlabel('Months')
    xlim([1,12])
end
legend(lgnd,'orientation','horizontal','location','North')
end