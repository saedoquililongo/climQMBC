# -*- coding: utf-8 -*-
"""
We hope this script will help the user understand how to use the main 
functionalities of the climQMBC package. Alongside this script, in a folder 
called Sample_data, the user will find sample datasets for daily and monthly 
precipitation (mm) and mean temperature (C). The datasets include modeled and
observed data. The files name structure is xxx_yyy_pp.csv, where xxx can be mod 
or obs for modeled or observed data, respectively, yyy can be tmp or pp 
for temperature or precipitation data, respectively, and pp can be D or M for
daily and monthly data, respectively. The historical period of the sample
datasets begin in 1985 and has a length of 30 years (1985 to 2014, including
both years), and the modeled period begins in 1985 and has a length 116 years 
(1985 to 2100, including both years).

This script is divided into four examples (three with monthly data and one with
daily data. The first two examples can be used to evaluate the performance of 
each bias correction method available in the climQMBC package. Example 1 shows
how to use the report function with the minimum number of inputs. The five
methods available in the climQMBC package will be reported, and the future
projected windows displayed are the first period after the historical period
and the last period before the end of the modeled period. Remember that the
projected periods length is equal to the length of the historical period.
Example 2 shows how to use the report function for specific bias correction
methods and projected periods. The report will analyze the projected periods
centered in 2035 2060 and 2080. Example 3 and 4 shows how each bias correction
method available in the climQMBC package should be called for monthly and daily 
frequency. The outputs of each function are columns vector with daily or monthly 
corrected data.

Feel free to uncomment each example, modify the periods and try your own 
datasets.


Written by Sebastian Aedo Quililongo (1*)
           Cristian Chadwick         (2)
           Fernando Gonzalez-Leiva   (3)
           Jorge Gironas             (3, 4)
           
  (1) Stockholm Environment Institute, Latin America Centre, Bogota, Colombia
  (2) Faculty of Engineering and Sciences, Universidad Adolfo Ibanez, Santiago,
      Chile
  (3) Department of Hydraulics and Environmental Engineering, Pontificia
      Universidad Catolica de Chile, Santiago, Chile
  (4) Centro de Cambio Global UC, Pontificia Universidad Catolica de Chile,
      Santiago, Chile
      
*Maintainer contact: sebastian.aedo.q@gmail.com
Revision: 1, updated Apr 2024
"""

from climQMBC.methods import QM, DQM, QDM, UQM, SDM
from climQMBC.report import report
import pandas as pd


# =============================================================================
# I) Montly and annual data
# =============================================================================
# variable:
#    - pr  (precipitation)
#    - tas (temperature)
# allow_negatives:
#    - 0 (variables like precipitation)
#    - 1 (variables like temperature)
# mult_change:
#    - 0 (additive change: fut = hist + delta) 
#    - 1 (multiplicative change: fut = hist*delta)
# SDM_var: (for Scaled Distribution Mapping only)
#    - 0 (temperature: normal distribution and additive changes) 
#    - 1 (precipitation: gamma distribution and multiplicative changes)
# frq:
#    - 'D': Daily data (use in section II. Section I works for 'M' and 'A')
#    - 'M': Monthly data (report function works only with 'M')
#    - 'A': Anual data

variable = 'pr'
mult_change = 0
allow_negatives = 1
SDM_var = 1

# Load observed and model data.
obs = pd.read_csv(f'../Sample_data/obs_{variable}_M.csv')[[variable]].values
mod = pd.read_csv(f'../Sample_data/mod_{variable}_M.csv')[[variable]].values

# Example 1
# Example 1 shows how to use the report function with the minimum number
# of inputs. The five methods available in the climQMBC package will be
# reported, and the future projected windows displayed are the first
# period after the historical period and the last period before the end
# of the modeled period. Remember that the projected periods length is
# equal to the length of the historical period.

qm_series,dqm_series,qdm_series,uqm_series,sdm_series = report(obs, mod, SDM_var=SDM_var, mult_change=mult_change, allow_negatives=allow_negatives)


# Example 2
# Example 2 shows how to use the report function for specific bias
# correction methods and projected periods. The Quantile Delta Mapping
# (QDM), Unbiased Quantile Mapping (UQM), and Scaled Distribution Mapping
# (SDM) methods will be reported. The report will analyze the projected periods
# centered in 2035 2060 and 2080.

# qm_series,dqm_series,qdm_series,uqm_series,sdm_series = report(obs, mod, SDM_var=SDM_var, mult_change=mult_change, allow_negatives=allow_negatives, y_init = 1979,y_wind = [2035,2060,2080])


# Example 3
# Example 3 shows how each bias correction method available in the
# climQMBC package should be called. The outputs of each function are
# columns vector with monthly corrected data.

# frq = 'A' # 'M' or 'A''
# qm_series = QM(obs, mod, allow_negatives=allow_negatives, frq=frq)
# dqm_series = DQM(obs, mod, mult_change=mult_change, allow_negatives=allow_negatives, frq=frq)
# qdm_series = QDM(obs, mod, mult_change=mult_change, allow_negatives=allow_negatives, frq=frq)
# uqm_series = UQM(obs, mod, mult_change=mult_change, allow_negatives=allow_negatives, frq=frq)
# sdm_series = SDM(obs, mod, SDM_var=SDM_var, frq=frq)



# %%
# =============================================================================
# II) Daily data
# =============================================================================
# variable:
#    - pr  (precipitation)
#    - tas (temperature)
# allow_negatives:
#    - 0 (variables like precipitation)
#    - 1 (variables like temperature)
# mult_change:
#    - 0 (additive change: fut = hist + delta) 
#    - 1 (multiplicative change: fut = hist*delta)
# frq:
#    - 'D': Daily data
#    - 'M': Monthly data
#    - 'A': Anual data
# SDM_var: (for Scaled Distribution Mapping only)
#    - 0 (temperature: normal distribution and additive changes) 
#    - 1 (precipitation: gamma distribution and multiplicative changes)
# day_win: An integer to define a moving window for each day of the year and 
#          compute the statistics for each probability distribution function
#          and projected change. The lenght of the window is computed as 2*win-1
# pp_threshold: A float to define the threshold to consider rain or no-rain values
# pp_factor: A float to scale pp_threshold and set as limit of the random low
#            values to replace no-rain values
frq = 'D'

variable = 'tas'
mult_change = 1
allow_negatives = 0
SDM_var = 0

day_win = 15
pp_threshold=1
pp_factor=1/10000

# Load observed and model data.
obs = pd.read_csv(f'../Sample_data/obs_{variable}_D.csv')[[variable]].values
mod = pd.read_csv(f'../Sample_data/mod_{variable}_D.csv')[[variable]].values

# Example 4
# Example 4 shows how each bias correction method available in the
# climQMBC package should be called. The outputs of each function are
# columns vector with daily corrected data.

# qm_series = QM(obs, mod, allow_negatives=allow_negatives, frq=frq, day_win=day_win, pp_threshold=pp_threshold, pp_factor=pp_factor)
# dqm_series = DQM(obs, mod, allow_negatives=allow_negatives, frq=frq, mult_change=mult_change, day_win=day_win, pp_threshold=pp_threshold, pp_factor=pp_factor)
# qdm_series = QDM(obs, mod, allow_negatives=allow_negatives, frq=frq, mult_change=mult_change, day_win=day_win, pp_threshold=pp_threshold, pp_factor=pp_factor)
# uqm_series = UQM(obs, mod, allow_negatives=allow_negatives, frq=frq, mult_change=mult_change, day_win=day_win, pp_threshold=pp_threshold, pp_factor=pp_factor)
# sdm_series = SDM(obs, mod, SDM_var, frq=frq, pp_threshold=pp_threshold, pp_factor=pp_factor, day_win=day_win)