%%Description
% We hope this script will help the user understand how to use the main
% functionalities of the climQMBC package. Alongside this script, in a
% folder called Sample_data, the user will find sample datasets for monthly
% precipitation (mm) and mean temperature (C). The datasets include
% modeled and observed data. The files name structure is xxx_yyy.csv, where
% xxx can be mod or obs for modeled or observed data, respectively, and yyy
% can be tmp or pp for temperature or precipitation data, respectively. The
% historical period of the sample datasets begin in 1979 and has a length
% of 36 years (1979 to 2014, including both years), and the modeled period
% begins in 1979 and has a length of 122 years (1979 to 2100, including
% both years).
% 
% This script is divided into three examples. The first two examples can be
% used to evaluate the performance of each bias correction method available
% in the climQMBC package. Example 1 shows how to use the report function
% with the minimum number of inputs. The five methods available in the
% climQMBC package will be reported, and the future projected windows
% displayed are the first period after the historical period and the last
% period before the end of the modeled period. Remember that the projected
% periods length is equal to the length of the historical period. Example 2
% shows how to use the report function for specific bias correction methods
% and projected periods. The Quantile Delta Mapping (QDM), Unbiased
% Quantile Mapping (UQM), and Scaled Distribution Mapping (SDM) methods
% will be reported. The report will analyze the projected periods centered
% in 2035 2060 and 2080. Example 3 shows how each bias correction method
% available in the climQMBC package should be called. The outputs of each
% function are columns vector with monthly corrected data.
% 
% Feel free to uncomment each example, modify the periods and try your own
% datasets.
%  
% 
% Written by Sebastian Aedo Quililongo (1*)
%            Cristian Chadwick         (2)
%            Fernando Gonzalez-Leiva   (3)
%            Jorge Gironas             (3, 4)
%            
%   (1) Stockholm Environment Institute, Latin America Centre, Bogota,
%       Colombia
%   (2) Faculty of Engineering and Sciences, Universidad Adolfo Ibanez,
%       Santiago, Chile
%   (3) Department of Hydraulics and Environmental Engineering, Pontificia
%       Universidad Catolica de Chile, Santiago, Chile
%   (4) Centro de CAmbio Global UC, Pontificia Universidad Catolica de 
%       Chile, Santiago, Chile
%       
%   *Maintainer contact: sebastian.aedo.q@gmail.com
% Revision: 1, updated Jul 2022

%% Set path to the climQMBC package and import datasets
clc
clear all

addpath('./climQMBC')  

%  =============================================================================
%  I) Montly and annual data
%  =============================================================================
%  variable:
%     - pr  (precipitation)
%     - tas (temperature)
%  allow_negatives:
%     - 0 (variables like precipitation)
%     - 1 (variables like temperature)
%  mult_change:
%     - 0 (additive change: fut = hist + delta) 
%     - 1 (multiplicative change: fut = hist*delta)
%  SDM_var: (for Scaled Distribution Mapping only)
%     - 0 (temperature: normal distribution and additive changes) 
%     - 1 (precipitation: gamma distribution and multiplicative changes)
%  frq:
%     - 'D': Daily data (use in section II. Section I works for 'M' and 'A')
%     - 'M': Monthly data (report function works only with 'M')
%     - 'A': Anual data

variable = 'pr';
mult_change = 1;
allow_negatives = 0;
SDM_var = 1;

% Load observed and model data.
obs = csvread(strcat('../Sample_data/obs_',variable,'_M.csv'),1,3);
mod = csvread(strcat('../Sample_data/mod_',variable,'_M.csv'),1,3);


%% Example 1
%   Example 1 shows how to use the report function with the minimum number
%   of inputs. The five methods available in the climQMBC package will be
%   reported, and the future projected windows displayed are the first
%   period after the historical period and the last period before the end
%   of the modeled period. Remember that the projected periods length is
%   equal to the length of the historical period.

[QM_series,DQM_series,QDM_series,UQM_series,SDM_series] = report(obs,mod,SDM_var,mult_change,allow_negatives);


%% Example 2
%   Example 2 shows how to use the report function for specific bias
%   correction methods and projected periods. The Quantile Delta Mapping
%   (QDM), Unbiased Quantile Mapping (UQM), and Scaled Distribution Mapping
%   (SDM) methods will be reported. The report will analyze the projected
%   periods centered in 2035 2060 and 2080.

% [QM_series,DQM_series,QDM_series,UQM_series,SDM_series] = report(obs,mod,SDM_var,mult_change,allow_negatives,{'QDM','UQM','SDM'},1979,[2035 2060 2080]);


%% Example 3
%   Example 3 shows how each bias correction method available in the
%   climQMBC package should be called. The outputs of each function are
%   columns vector with monthly corrected data.

% frq = 'A'; % 'M' for monthly data; 'A' for annual data
% QM_series = QM(obs,mod,allow_negatives,frq);
% DQM_series = DQM(obs,mod,mult_change,allow_negatives,frq);
% QDM_series = QDM(obs,mod,mult_change,allow_negatives,frq);
% UQM_series = UQM(obs,mod,mult_change,allow_negatives,frq);
% sdm_series = SDM(obs,mod,SDM_var,frq);


% %% AA
obs = csvread('../Sample_data/obs_tas_D.csv', 1,4);
mod = csvread('../Sample_data/mod_tas_D.csv', 1,4);

frq = 'D';
allow_negatives=0;
mult_change = 1;
pp_threshold=1;
pp_factor=1/(100*100);
day_win=15;
SDM_var = 0;


% qm_series = QM(obs,mod,allow_negatives,frq,pp_threshold, pp_factor, day_win,true,3,3);
% dqm_series = DQM(obs,mod,mult_change,allow_negatives,frq,pp_threshold, pp_factor, day_win,true,1,1);
% qdm_series = QDM(obs,mod,mult_change,allow_negatives,frq,pp_threshold,pp_factor,2,pp_threshold,day_win);
% uqm_series = UQM(obs,mod,mult_change,allow_negatives,frq,pp_threshold, pp_factor, day_win);
% sdm_series = SDM(obs,mod,SDM_var,frq,pp_threshold,pp_factor,day_win);

% plot(uqm_series,'o')

% csvwrite('../../SDM.txt',sdm_series);