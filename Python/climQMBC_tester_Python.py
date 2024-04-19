# -*- coding: utf-8 -*-
"""
We hope this script will help the user understand how to use the main 
functionalities of the climQMBC package. Alongside this script, in a folder 
called Sample_data, the user will find sample datasets for monthly 
precipitation (mm) and mean temperature (C). The datasets include modeled and
observed data. The files name structure is xxx_yyy.csv, where xxx can be mod 
or obs for modeled or observed data, respectively, and yyy can be tmp or pp 
for temperature or precipitation data, respectively. The historical period of 
the sample datasets begin in 1979 and has a length of 36 years (1979 to 2014,
including both years), and the modeled period begins in 1979 and has a length
of 122 years (1979 to 2100, including both years).

This script is divided into three examples. The first two examples can be used
to evaluate the performance of each bias correction method available in the 
climQMBC package.Example 1 shows how to use the report function with the 
minimum number of inputs. The five methods available in the climQMBC package 
will be reported, and the future projected windows displayed are the first 
period after the historical period and the last period before the end of the 
modeled period. Remember that the projected periods length is equal to the 
length of the historical period. Example 2 shows how to use the report 
function for specific bias correction methods and projected periods. The 
Quantile Delta Mapping (QDM), Unbiased Quantile Mapping (UQM), and Scaled 
Distribution Mapping (SDM) methods will be reported. The report will analyze 
the projected periods centered in 2035 2060 and 2080. Example 3 shows how each 
bias correction method available in the climQMBC package should be called. The 
outputs of each function are columns vector with monthly corrected data.

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
import matplotlib.pylab as plt
import pandas as pd
import numpy as np


# =============================================================================
# I) Montly and annual data
# =============================================================================
# variable:
#    - pp
#    - tmax
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

variable = 'pp'
mult_change = 1
allow_negatives = 0
SDM_var = 1

# Load observed and model data.
obs = np.array(pd.read_csv('../Sample_data/obs_'+variable+'.csv',header=None))
mod = np.array(pd.read_csv('../Sample_data/mod_'+variable+'.csv',header=None))

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
# (SDM) methods will be reported. The report will analyze the projected

# qm_series,dqm_series,qdm_series,uqm_series,sdm_series = report(obs, mod, SDM_var=SDM_var, mult_change=mult_change, allow_negatives=allow_negatives, y_init = 1979,y_wind = [2035,2060,2080])

# Example 3
# Example 3 shows how each bias correction method available in the
# climQMBC package should be called. The outputs of each function are
# columns vector with monthly corrected data.

frq = 'M' # 'M' or 'A''
# qm_series = QM(obs, mod, allow_negatives=allow_negatives, frq=frq)
# dqm_series = DQM(obs, mod, mult_change=mult_change, allow_negatives=allow_negatives, frq=frq)
# qdm_series = QDM(obs, mod, mult_change=mult_change, allow_negatives=allow_negatives, frq=frq)
# uqm_series = UQM(obs, mod, mult_change=mult_change, allow_negatives=allow_negatives, frq=frq)
# sdm_series = SDM(obs, mod, SDM_var=1, frq=frq)



# =============================================================================
# II) Daily data
# =============================================================================
# variable:
#    - pp
#    - tmax
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
# day_win: An integer to define a moving window for each day of the year and 
#          compute the statistics for each probability distribution function
#          and projected change. The lenght of the window is computed as 2*win-1
# pp_threshold: A float to define the threshold to consider rain or no-rain values
# pp_factor: A float to scale pp_threshold and set as limit of the random low
#            values to replace no-rain values
variable = 'pp'
mult_change = 1
allow_negatives = 0
frq = 'D'

win = 5
pp_threshold=1
pp_factor=1/10000

# Load observed and model data.
obs = pd.read_csv(f'../Sample_data/obs_{variable}_D.csv')[[variable]].values
mod = pd.read_csv(f'../Sample_data/mod_{variable}_D.csv')[[variable]].values


qm_series = QM(obs, mod, allow_negatives=allow_negatives, frq=frq, win=win, pp_threshold=pp_threshold, pp_factor=pp_factor)
dqm_series = DQM(obs, mod, allow_negatives=allow_negatives, frq=frq, mult_change=mult_change, win=win, pp_threshold=pp_threshold, pp_factor=pp_factor)
qdm_series = QDM(obs, mod, allow_negatives=allow_negatives, frq=frq, mult_change=mult_change, win=win, pp_threshold=pp_threshold, pp_factor=pp_factor)
uqm_series = UQM(obs, mod, allow_negatives=allow_negatives, frq=frq, mult_change=mult_change, win=win, pp_threshold=pp_threshold, pp_factor=pp_factor)
# sdm_series = SDM(obs, mod, SDM_var, frq='A', pp_threshold=1, pp_factor=1/100)

# %%
plt.figure()

plt.title('QM')
plt.subplot(2,2,1)
plt.plot(qm_series)
plt.plot(obs[:,0])
plt.ylim(0,100)

plt.title('DQM')
plt.subplot(2,2,2)
plt.plot(dqm_series)
plt.plot(obs[:,0])
plt.ylim(0,100)

plt.title('QDM')
plt.subplot(2,2,3)
plt.plot(qdm_series)
plt.plot(obs[:,0])
plt.ylim(0,100)

plt.title('UQM')
plt.subplot(2,2,4)
plt.plot(uqm_series)
plt.plot(obs[:,0])
plt.ylim(0,100)

plt.tight_layout()
plt.show()

# %%

plt.figure(dpi=300, figsize=(8,8))
plt.suptitle('Serie diaria', fontweight='bold')

plt.subplot(3,2,1)
plt.title('QM')
plt.plot(mod[:,0], label='Mod');plt.plot(qm_series, label='Cor');plt.plot(obs[:,0], label='Obs')
plt.legend(ncol=3, loc='upper center')

plt.subplot(3,2,2)
plt.title('DQM')
plt.plot(mod[:,0]);plt.plot(dqm_series);plt.plot(obs[:,0])

plt.subplot(3,2,3)
plt.title('QDM')
plt.plot(mod[:,0]);plt.plot(qdm_series);plt.plot(obs[:,0])

plt.subplot(3,2,4)
plt.title('UQM')
plt.plot(mod[:,0]);plt.plot(uqm_series);plt.plot(obs[:,0])

plt.subplot(3,2,(5,6))
plt.title('Serie anual', fontweight='bold')
plt.plot(mod.reshape((365,80), order='F').sum(0), label='Mod')
plt.plot(qm_series.reshape((365,80), order='F').sum(0), label='QM')
plt.plot(dqm_series.reshape((365,80), order='F').sum(0), label='DQM')
plt.plot(qdm_series.reshape((365,80), order='F').sum(0), label='QDM')
plt.plot(uqm_series.reshape((365,80), order='F').sum(0), label='UQM')
plt.plot(obs.reshape((365,30), order='F').sum(0), label='Obs')
plt.plot(mod.reshape((365,80), order='F').sum(0)*obs.reshape((365,30), order='F').sum(0).mean()/mod.reshape((365,80), order='F').sum(0)[:30].mean(), label='Obj')
plt.legend()

plt.tight_layout()
plt.show()

# %%
win = 15
mod_series_d = np.vstack([mod.reshape((365,80), order='F')[-win:],np.tile(mod.reshape((365,80), order='F'),(win*2,1)),mod.reshape((365,80), order='F')[:win]])
mod_series_d = mod_series_d.reshape(mod.reshape((365,80), order='F').shape[0]+1,win*2,mod.reshape((365,80), order='F').shape[1], order='F')[:-1,1:]

obs_series_d = np.vstack([obs.reshape((365,30), order='F')[-win:],np.tile(obs.reshape((365,30), order='F'),(win*2,1)),obs.reshape((365,30), order='F')[:win]])
obs_series_d = obs_series_d.reshape(obs.reshape((365,30), order='F').shape[0]+1,win*2,obs.reshape((365,30), order='F').shape[1], order='F')[:-1,1:]

qm_series_d = np.vstack([qm_series.reshape((365,80), order='F')[-win:],np.tile(qm_series.reshape((365,80), order='F'),(win*2,1)),qm_series.reshape((365,80), order='F')[:win]])
qm_series_d = qm_series_d.reshape(qm_series.reshape((365,80), order='F').shape[0]+1,win*2,qm_series.reshape((365,80), order='F').shape[1], order='F')[:-1,1:]

plt.figure(dpi=300, figsize=(8,6))
plt.suptitle('Estacionalidad diaria histórica', fontweight='bold')

plt.subplot(2,2,1)
plt.title('Promedio diario')
plt.plot(mod.reshape((365,80), order='F')[:,:30].mean(1), label='Mod')
plt.plot(qm_series.reshape((365,80), order='F')[:,:30].mean(1), label='QM_rainday')
plt.plot(obs.reshape((365,30), order='F')[:,:30].mean(1), label='Obs')
plt.legend()

plt.subplot(2,2,2)
plt.title('Desviación estándar diaria')
plt.plot(mod.reshape((365,80), order='F')[:,:30].std(1))
plt.plot(qm_series.reshape((365,80), order='F')[:,:30].std(1))
plt.plot(obs.reshape((365,30), order='F')[:,:30].std(1))

plt.subplot(2,2,3)
plt.title('Promedio diario móvil')
plt.plot(mod_series_d[:,:,:30].mean((1,2)))
plt.plot(qm_series_d[:,:,:30].mean((1,2)))
plt.plot(obs_series_d[:,:,:30].mean((1,2)))

plt.subplot(2,2,4)
plt.title('Desviación estándar diaria móvil')
plt.plot(mod_series_d[:,:,:30].std((1,2)))
plt.plot(qm_series_d[:,:,:30].std((1,2)))
plt.plot(obs_series_d[:,:,:30].std((1,2)))

plt.tight_layout()
plt.show()

# %%
win = 15
qm_series_d = np.vstack([qm_series.reshape((365,80), order='F')[-win:],np.tile(qm_series.reshape((365,80), order='F'),(win*2,1)),qm_series.reshape((365,80), order='F')[:win]])
qm_series_d = qm_series_d.reshape(qm_series.reshape((365,80), order='F').shape[0]+1,win*2,qm_series.reshape((365,80), order='F').shape[1], order='F')[:-1,1:]

dqm_series_d = np.vstack([dqm_series.reshape((365,80), order='F')[-win:],np.tile(dqm_series.reshape((365,80), order='F'),(win*2,1)),dqm_series.reshape((365,80), order='F')[:win]])
dqm_series_d = dqm_series_d.reshape(dqm_series.reshape((365,80), order='F').shape[0]+1,win*2,dqm_series.reshape((365,80), order='F').shape[1], order='F')[:-1,1:]

qdm_series_d = np.vstack([qdm_series.reshape((365,80), order='F')[-win:],np.tile(qdm_series.reshape((365,80), order='F'),(win*2,1)),qdm_series.reshape((365,80), order='F')[:win]])
qdm_series_d = qdm_series_d.reshape(qdm_series.reshape((365,80), order='F').shape[0]+1,win*2,qdm_series.reshape((365,80), order='F').shape[1], order='F')[:-1,1:]

uqm_series_d = np.vstack([uqm_series.reshape((365,80), order='F')[-win:],np.tile(uqm_series.reshape((365,80), order='F'),(win*2,1)),uqm_series.reshape((365,80), order='F')[:win]])
uqm_series_d = uqm_series_d.reshape(uqm_series.reshape((365,80), order='F').shape[0]+1,win*2,uqm_series.reshape((365,80), order='F').shape[1], order='F')[:-1,1:]

plt.figure(dpi=300, figsize=(8,6))
plt.suptitle('Estacionalidad diaria futura', fontweight='bold')

plt.subplot(2,2,1)
plt.title('Promedio diario')
plt.plot(mod.reshape((365,80), order='F')[:,30:].mean(1), label='Mod')
plt.plot(qm_series.reshape((365,80), order='F')[:,30:].mean(1), label='QM')
plt.plot(dqm_series.reshape((365,80), order='F')[:,30:].mean(1), label='DQM')
plt.plot(qdm_series.reshape((365,80), order='F')[:,30:].mean(1), label='QDM')
plt.plot(uqm_series.reshape((365,80), order='F')[:,30:].mean(1), label='UQM')
plt.plot(obs.reshape((365,30), order='F')[:,:30].mean(1)*mod.reshape((365,80), order='F')[:,30:].mean(1)/mod.reshape((365,80), order='F')[:,:30].mean(1), label='Obj')
plt.legend()

plt.subplot(2,2,2)
plt.title('Desviación estándar diaria')
plt.plot(mod.reshape((365,80), order='F')[:,30:].std(1))
plt.plot(qm_series.reshape((365,80), order='F')[:,30:].std(1))
plt.plot(dqm_series.reshape((365,80), order='F')[:,30:].std(1))
plt.plot(qdm_series.reshape((365,80), order='F')[:,30:].std(1))
plt.plot(uqm_series.reshape((365,80), order='F')[:,30:].std(1))
plt.plot(obs.reshape((365,30), order='F')[:,:30].std(1)*mod.reshape((365,80), order='F')[:,30:].std(1)/mod.reshape((365,80), order='F')[:,:30].std(1))

plt.subplot(2,2,3)
plt.title('Promedio diario móvil')
plt.plot(mod_series_d[:,:,30:].mean((1,2)))
plt.plot(qm_series_d[:,:,30:].mean((1,2)))
plt.plot(dqm_series_d[:,:,30:].mean((1,2)))
plt.plot(qdm_series_d[:,:,30:].mean((1,2)))
plt.plot(uqm_series_d[:,:,30:].mean((1,2)))
plt.plot(obs_series_d[:,:,:30].mean((1,2))*mod_series_d[:,:,30:].mean((1,2))/mod_series_d[:,:,:30].mean((1,2)))

plt.subplot(2,2,4)
plt.title('Desviación estándar diaria móvil')
plt.plot(mod_series_d[:,:,30:].std((1,2)))
plt.plot(qm_series_d[:,:,30:].std((1,2)))
plt.plot(dqm_series_d[:,:,30:].std((1,2)))
plt.plot(qdm_series_d[:,:,30:].std((1,2)))
plt.plot(uqm_series_d[:,:,30:].std((1,2)))
plt.plot(obs_series_d[:,:,:30].std((1,2))*mod_series_d[:,:,30:].std((1,2))/mod_series_d[:,:,:30].std((1,2)))

plt.tight_layout()
plt.show()


# %%


plt.figure(dpi=300, figsize=(8,6))
plt.suptitle('Tasas de cambio diaria', fontweight='bold')

plt.subplot(2,2,1)
plt.title('Promedio diario')
plt.plot(mod.reshape((365,80), order='F')[:,30:].mean(1)/mod.reshape((365,80), order='F')[:,:30].mean(1),label='Mod')
plt.plot(qm_series.reshape((365,80), order='F')[:,30:].mean(1)/qm_series.reshape((365,80), order='F')[:,:30].mean(1), label='QM')
plt.plot(dqm_series.reshape((365,80), order='F')[:,30:].mean(1)/dqm_series.reshape((365,80), order='F')[:,:30].mean(1), label='DQM')
plt.plot(qdm_series.reshape((365,80), order='F')[:,30:].mean(1)/qdm_series.reshape((365,80), order='F')[:,:30].mean(1),label='QDM')
plt.plot(uqm_series.reshape((365,80), order='F')[:,30:].mean(1)/uqm_series.reshape((365,80), order='F')[:,:30].mean(1),label='UQM')
plt.legend()

plt.subplot(2,2,2)
plt.title('Desviación estándar diaria')
plt.plot(mod.reshape((365,80), order='F')[:,30:].std(1)/mod.reshape((365,80), order='F')[:,:30].std(1))
plt.plot(qm_series.reshape((365,80), order='F')[:,30:].std(1)/qm_series.reshape((365,80), order='F')[:,:30].std(1))
plt.plot(dqm_series.reshape((365,80), order='F')[:,30:].std(1)/dqm_series.reshape((365,80), order='F')[:,:30].std(1))
plt.plot(qdm_series.reshape((365,80), order='F')[:,30:].std(1)/qdm_series.reshape((365,80), order='F')[:,:30].std(1))
plt.plot(uqm_series.reshape((365,80), order='F')[:,30:].std(1)/uqm_series.reshape((365,80), order='F')[:,:30].std(1))

plt.subplot(2,2,3)
plt.title('Promedio diario móvil')
plt.plot(mod_series_d[:,:,30:].mean((1,2))/mod_series_d[:,:,:30].mean((1,2)))
plt.plot(qm_series_d[:,:,30:].mean((1,2))/qm_series_d[:,:,:30].mean((1,2)))
plt.plot(dqm_series_d[:,:,30:].mean((1,2))/dqm_series_d[:,:,:30].mean((1,2)))
plt.plot(qdm_series_d[:,:,30:].mean((1,2))/qdm_series_d[:,:,:30].mean((1,2)))
plt.plot(uqm_series_d[:,:,30:].mean((1,2))/uqm_series_d[:,:,:30].mean((1,2)))

plt.subplot(2,2,4)
plt.title('Desviación estándar diaria móvil')
plt.plot(mod_series_d[:,:,30:].std((1,2))/mod_series_d[:,:,:30].std((1,2)))
plt.plot(qm_series_d[:,:,30:].std((1,2))/qm_series_d[:,:,:30].std((1,2)))
plt.plot(dqm_series_d[:,:,30:].std((1,2))/dqm_series_d[:,:,:30].std((1,2)))
plt.plot(qdm_series_d[:,:,30:].std((1,2))/qdm_series_d[:,:,:30].std((1,2)))
plt.plot(uqm_series_d[:,:,30:].std((1,2))/uqm_series_d[:,:,:30].std((1,2)))

plt.tight_layout()
plt.show()


