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

%% Set path to the climQMBC package and import datasets
clc
clear all

addpath('./climQMBC')  
var_name = {'tmax','pp'};


% Load observed and model data. Remember that for temperature, var = 0, and
% for precipitation, var = 1.
var = 1;
obs = csvread(strcat('../Sample_data/obs_',char(var_name(var+1)),'.csv'));
mod = csvread(strcat('../Sample_data/mod_',char(var_name(var+1)),'.csv'));

%% Example 1
%   Example 1 shows how to use the report function with the minimum number
%   of inputs. The five methods available in the climQMBC package will be
%   reported, and the future projected windows displayed are the first
%   period after the historical period and the last period before the end
%   of the modeled period. Remember that the projected periods length is
%   equal to the length of the historical period.

[QM_series,DQM_series,QDM_series,UQM_series,SDM_series] = report(obs,mod,var);


%% Example 2
%   Example 2 shows how to use the report function for specific bias
%   correction methods and projected periods. The Quantile Delta Mapping
%   (QDM), Unbiased Quantile Mapping (UQM), and Scaled Distribution Mapping
%   (SDM) methods will be reported. The report will analyze the projected
%   periods centered in 2035 2060 and 2080.

% [QM_series,DQM_series,QDM_series,UQM_series,SDM_series] = report(obs,mod,var,{'QDM','UQM','SDM'},1979,[2035 2060 2080]);


%% Example 3
%   Example 3 shows how each bias correction method available in the
%   climQMBC package should be called. The outputs of each function are
%   columns vector with monthly corrected data.

% frq = 'M'; % 'M' for monthly data; 'A' for annual data
% QM_series = QM(obs,mod,var,frq);
% DQM_series = DQM(obs,mod,var,frq);
% QDM_series = QDM(obs,mod,var,frq);
% UQM_series = UQM(obs,mod,var,frq);
% SDM_series = SDM(obs,mod,var,frq);
