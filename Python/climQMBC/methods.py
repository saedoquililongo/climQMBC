# -*- coding: utf-8 -*-
"""
This script contains the bias correction methods implemented in the climQMBC
package, including the Quantile Mapping (QM), the Detrended Quantile Mapping 
(DQM), the Quantile Delta Mapping (QDM), the Unbiased Quantile Mapping (UQM), 
and the Scaled Distribution Mapping (SDM).

References:
    Cannon, A. J., S. R. Sobie, and T. Q. Murdock, 2015: Bias correction of 
    GCM precipitation by quantile mapping: How well do methods preserve changes
    in quantiles and extremes? J. Climate, 28(17), 6938-6959, 
    https://doi.org/10.1175/JCLI-D-14-00754.1
    
    Switanek, B. M., Troch, P., Castro, L. C., Leuprecht, A., Chang, H. I., 
    Mukherjee, R., and Demaria, M. C. E. (2017) Scaled distribution mapping: A 
    bias correction method that preserves raw climate model projected changes.
    Hydrology &amp; Earth System Sciences, 21, 2649-2666,
    https://doi.org/10.5194/hess-21-2649-2017.
    
    Chadwick, C., Gironás, J., González-Leiva, F., and Aedo, S. (2023). Bias 
    adjustment to preserve changes in variability: the unbiased mapping of GCM 
    changes. Hydrological Sciences Journal,
    https://doi.org/10.1080/02626667.2023.2201450


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
Revision: 1, updated Jul 2022
"""

from .utils import formatQM, getDist, getCDF, getCDFinv
from scipy.signal import detrend
import scipy.stats as stat
import numpy as np


def QM(obs, mod, var, frq='M', pp_threshold=1, pp_factor=1/100):
    """
    This function performs bias correction of modeled series based on observed
    data by the Quantile Mapping (QM) method, as described by Cannon et al. 
    (2015). Correction is performed to monthly or annual precipitation or 
    temperature data in a single location. An independent probability 
    distribution function is assigned to each month and to each projected 
    period based on the Kolmogorov-Smirnov test. If annual frequency is 
    specified, this is applied to the complete period. Each projected period 
    considers a window whose length is equal to the number of years of the 
    historical period and ends in the analyzed year.
    
    Description:
        0) Format inputs and get statistics of the observed and modeled series
           of the historical period.
           
        1) Assign a probability distribution function to each month for the 
           observed and modeled data in the historical period. If annual 
           frequency is specified, this is applied to the complete historical 
           period.
           
        2) Apply the cumulative distribution function of the modeled data, 
           evaluated with the statistics of the modeled data in the historical 
           period, to the modeled data.

        3) Apply the inverse cumulative distribution function of the observed 
           data, evaluated with the statistics of the observed data in the 
           historical period, to the probabilities obtained from 2).

    Inputs:
        obs:    A column vector of monthly or annual observed data (temperature
                or precipitation). If monthly frequency is specified, 
                the length of this vector is 12 times the number of observed 
                years [12 x y_obs, 1]. If annual frequency is specified, the 
                length of this vector is equal to the number of observed years
                [y_obs, 1].

        mod:    A column vector of monthly or annual modeled data (temperature
                or precipitation). If monthly frequency is specified, 
                the length of this vector is 12 times the number of observed
                years [12 x y_mod, 1]. If annual frequency is specified, the 
                length of this vector is equal to the number of observed years
                [y_mod, 1].

        var:    A flag that identifies if data are temperature or 
                precipitation. This flag tells the getDist function if it has 
                to discard distribution functions that allow negative numbers,
                and if the terms in the correction equations are
                multiplied/divided or added/subtracted.
                    Temperature:   var = 0
                    Precipitation: var = 1

        NOTE: This routine considers that obs and mod series start in january
        of the first year and ends in december of the last year.

    Optional inputs:
        frq:    A string specifying if the input is annual or monthly data. If
                not specified, it is set monthly as default.
                    Monthly:   frq = 'M'
                    Annual:    frq = 'A'

        pp_threshold:    A float indicating the threshold to consider 
                         physically null precipitation values.
                
        pp_factor:       A float indicating the maximum value of the random
                         values that replace physically null precipitation 
                         values.

    Output:
        QM_series:  A column vector of monthly or annual modeled data 
                    (temperature or precipitation) corrected by the QM method. 
                    If monthly frequency is specified, the length of this 
                    vector is 12 times the number of observed years 
                    [12 x y_mod, 1]. If annual frequency is specified, the 
                    length of this vector is equal to the number of observed 
                    years [y_mod, 1].
    
    References:
        Cannon, A. J., S. R. Sobie, and T. Q. Murdock, 2015: Bias correction of 
        GCM precipitation by quantile mapping: How well do methods preserve 
        changes in quantiles and extremes? J. Climate, 28(17), 6938-6959, 
        https://doi.org/10.1175/JCLI-D-14-00754.1

    """
    
    # 0) Format inputs and get statistics of the observed and modeled series of
    #    the historical period (formatQM function of the climQMBC package).
    y_obs, obs_series, mod_series, mu_obs, std_obs, skew_obs, skewy_obs, mu_mod, std_mod, skew_mod, skewy_mod = formatQM(obs, mod, var, frq, pp_threshold, pp_factor)
    
    # 1) Assign a probability distribution function to each month for the
    #    observed and modeled data in the historical period. If annual
    #    frequency is specified, this is applied to the complete historical
    #    period (getDist function of the climQMBC package).
    PDF_obs = getDist(obs_series, mu_obs, std_obs, skew_obs, skewy_obs, var)
    PDF_mod = getDist(mod_series[:,:y_obs], mu_mod, std_mod, skew_mod, skewy_mod, var)

    # 2) Apply the cumulative distribution function of the modeled data,
    #    evaluated with the statistics of the modeled data in the historical
    #    period, to the modeled data (getCDF function of the climQMBC package).
    #    Equation 1 of Cannon et al. (2015).
    Taot = getCDF(PDF_mod, mod_series, mu_mod, std_mod, skew_mod, skewy_mod)
    
    # 3) Apply the inverse cumulative distribution function of the observed
    #    data, evaluated with the statistics of the observed data in the
    #    historical period, to the probabilities obtained from 2) (getCDFinv
    #    function of the climQMBC package). Equation 1 of Cannon et al. (2015).
    QM_series = getCDFinv(PDF_obs, Taot, mu_obs, std_obs, skew_obs, skewy_obs)
    QM_series = QM_series.reshape(-1, order='F')
    if var==1:
        QM_series[QM_series<pp_threshold] = 0

    return QM_series


def DQM(obs, mod, var, frq='M', pp_threshold=1, pp_factor=1/100):
    """
    This function performs bias correction of modeled series based on observed
    data by the Detrended Quantile Mapping (DQM) method, as described by Cannon
    et al. (2015). Correction is performed to monthly or annual precipitation 
    or temperature data in a single location. An independent probability 
    distribution function is assigned to each month and to each projected 
    period based on the Kolmogorov-Smirnov test. If annual frequency is 
    specified, this is applied to the complete period. Each projected period 
    considers a window whose length is equal to the number of years of the 
    historical period and ends in the analyzed year.
    
    Correction of the historical period is performed by the Quantile Mapping 
    (QM) method, as described by Cannon et al. (2015).
    
    Description:
        0) Format inputs and get statistics of the observed and modeled series
           of the historical period.
           
        1) Assign a probability distribution function to each month for the 
           observed and modeled data in the historical period. If annual 
           frequency is specified, this is applied to the complete historical 
           period.
           
        2) Extract the long-term trend from the modeled data:
            a) Get the monthly mean of the historical period. If annual 
               frequency is specified, this is applied to the complete period.
            b) Get the monthly mean of each projected period. If annual 
              frequency is specified, this is applied to the complete period.
            
        3) Compute the linear scaled values.
        
        4) Apply the cumulative distribution function of the modeled data, 
           evaluated with the statistics of the modeled period, to the future 
           modeled data.
        
        5) Apply the inverse cumulative distribution function of the observed 
           data, evaluated with the statistics of the observed data in the 
           historical period, to the probabilities obtained from 4).
        
        6) Reimpose the trend to the values obtained in 5).
    
        End) Perform QM for the historical period.

    Inputs:
        obs:    A column vector of monthly or annual observed data (temperature
                or precipitation). If monthly frequency is specified, 
                the length of this vector is 12 times the number of observed 
                years [12 x y_obs, 1]. If annual frequency is specified, the 
                length of this vector is equal to the number of observed years
                [y_obs, 1].

        mod:    A column vector of monthly or annual modeled data (temperature
                or precipitation). If monthly frequency is specified, 
                the length of this vector is 12 times the number of observed
                years [12 x y_mod, 1]. If annual frequency is specified, the 
                length of this vector is equal to the number of observed years
                [y_mod, 1].

        var:    A flag that identifies if data are temperature or 
                precipitation. This flag tells the getDist function if it has 
                to discard distribution functions that allow negative numbers,
                and if the terms in the correction equations are
                multiplied/divided or added/subtracted.
                    Temperature:   var = 0
                    Precipitation: var = 1

        NOTE: This routine considers that obs and mod series start in january
        of the first year and ends in december of the last year.

    Optional inputs:
        frq:    A string specifying if the input is annual or monthly data. If
                not specified, it is set monthly as default.
                    Monthly:   frq = 'M'
                    Annual:    frq = 'A'
                    
        pp_threshold:    A float indicating the threshold to consider 
                         physically null precipitation values.
                
        pp_factor:       A float indicating the maximum value of the random
                         values that replace physically null precipitation 
                         values.

    Output:
        DQM_series:  A column vector of monthly or annual modeled data 
                    (temperature or precipitation) corrected by the DQM method. 
                    If monthly frequency is specified, the length of this 
                    vector is 12 times the number of observed years 
                    [12 x y_mod, 1]. If annual frequency is specified, the 
                    length of this vector is equal to the number of observed 
                    years [y_mod, 1].
    
    References:
        Cannon, A. J., S. R. Sobie, and T. Q. Murdock, 2015: Bias correction of 
        GCM precipitation by quantile mapping: How well do methods preserve
        changes in quantiles and extremes? J. Climate, 28(17), 6938-6959, 
        https://doi.org/10.1175/JCLI-D-14-00754.1
    
    """
    
    # 0) Format inputs and get statistics of the observed and modeled series of
    #    the historical period (formatQM function of the climQMBC package).
    y_obs, obs_series, mod_series, mu_obs, std_obs, skew_obs, skewy_obs, mu_mod, std_mod, skew_mod, skewy_mod = formatQM(obs, mod, var, frq, pp_threshold, pp_factor)
    
    # 1) Assign a probability distribution function to each month for the
    #    observed and modeled data in the historical period. If annual
    #    frequency is specified, this is applied to the complete historical
    #    period (getDist function of the climQMBC package).
    PDF_obs = getDist(obs_series, mu_obs, std_obs, skew_obs, skewy_obs, var)
    PDF_mod = getDist(mod_series[:,:y_obs], mu_mod, std_mod, skew_mod, skewy_mod, var)    
        
    # 2) Extract the long-term trend from the modeled data:
    #    a) Get the monthly mean of the historical period. If annually
    #       frequency is specified, this is applied to the complete period).
    xbarmh = np.nanmean(mod_series[:,:y_obs], 1)
    
    #    b) Get the monthly mean of each projected period. If annually
    #       frequency is specified, this is applied to the complete period).
    y_mod = mod_series.shape[1]
    xbarmp = np.zeros((mod_series.shape[0], y_mod-y_obs))
    for j in range(y_mod-y_obs):
        xbarmp[:,j] = np.nanmean(mod_series[:,j+1:y_obs+j+1], 1)
    
    # 3)Compute the linear scaled values (value in square brackets in equation
    #   2 of Cannon et al. (2015)).
    LS = np.zeros(xbarmp.shape)
    if var==1: # Precipitacion
        for m in range(mod_series.shape[0]):
            LS[m,:] = mod_series[m,y_obs:]*xbarmh[m]/xbarmp[m]
    else:        # Temperature
        for m in range(mod_series.shape[0]):
            LS[m,:] = mod_series[m,y_obs:] + xbarmh[m] - xbarmp[m]
    
    # 4) Apply the cumulative distribution function of the modeled data,
    #    evaluated with the statistics of the modeled period, to the future
    #    modeled data (getCDF function of the climQMBC package). Equation 2 of
    #    Cannon et al. (2015).
    Taot = getCDF(PDF_mod, LS, mu_mod, std_mod, skew_mod, skewy_mod)
    
    # 5) Apply the inverse cumulative distribution function of the observed
    #    data, evaluated with the statistics of the observed data in the
    #    historical period, to the probabilities obtained from 4) (getCDFinv
    #    function of the climQMBC package). Equation 2 of Cannon et al. (2015).
    DQM_LS = getCDFinv(PDF_obs, Taot, mu_obs, std_obs, skew_obs, skewy_obs)
    
    # 6) Reimpose the trend to the values obtained in 5). Equation 2 of Cannon
    #    et al. (2015).
    if var==1: # Precipitation
        DQM = DQM_LS*xbarmp/np.tile(xbarmh,(y_mod-y_obs, 1)).T
    else: # Temperature
        DQM = DQM_LS + xbarmp - np.tile(xbarmh,(y_mod-y_obs, 1)).T
    
    DQM = DQM.reshape(-1, order='F')
    
    # 7) Perform QM for the historical period.
    mod_h = mod_series[:,:y_obs].reshape(-1, order='F')
    QM_series = QM(obs, mod_h, var, frq)
    DQM_series = np.hstack([QM_series, DQM])
    if var==1:
        DQM_series[DQM_series<pp_threshold] = 0
    
    return DQM_series


def QDM(obs, mod, var, frq='M', pp_threshold=1, pp_factor=1/100, rel_change_th=2, inv_mod_th=None):
    """
    This function performs bias correction of modeled series based on observed
    data by the Quantile Delta Mapping (QDM) method, as described by Cannon
    et al. (2015). Correction is performed to monthly or annual precipitation 
    or temperature data in a single location. An independent probability 
    distribution function is assigned to each month and to each projected 
    period based on the Kolmogorov-Smirnov test. If annual frequency is 
    specified, this is applied to the complete period. Each projected period 
    considers a window whose length is equal to the number of years of the 
    historical period and ends in the analyzed year.
    
    Correction of the historical period is performed by the Quantile Mapping 
    (QM) method, as described by Cannon et al. (2015).
    
    Description:
        0) Format inputs and get statistics of the observed and modeled series
           of the historical period.
           
        1) Assign a probability distribution function to each month for the 
           observed and modeled data in the historical period. If annual 
           frequency is specified, this is applied to the complete historical 
           period.
           
        2) For each projected period:
            a) Assign a probability distribution function to each month. 
               If annual frequency is specified, this is applied to the 
               complete period.
            b) Apply the cumulative distribution function of the projected 
               period, evaluated with the statistics of this period, to the 
               last data of the period.

        3) Apply the inverse cumulative distribution function:
            a) Of the observed data, evaluated with the statistics of the 
               observed data in the historical period, to the probabilities
               obtained from 2b).
            b) Of the modeled data, evaluated with the statistics of the 
               observed data in the historical period, to the probabilities
               obtained from 2b).

        4) Get the delta factor or relative change and apply it to the value
           obtained in 3b).
    
        End) Perform QM for the historical period.

    Inputs:
        obs:    A column vector of monthly or annual observed data (temperature
                or precipitation). If monthly frequency is specified, 
                the length of this vector is 12 times the number of observed 
                years [12 x y_obs, 1]. If annual frequency is specified, the 
                length of this vector is equal to the number of observed years
                [y_obs, 1].

        mod:    A column vector of monthly or annual modeled data (temperature
                or precipitation). If monthly frequency is specified, 
                the length of this vector is 12 times the number of observed
                years [12 x y_mod, 1]. If annual frequency is specified, the 
                length of this vector is equal to the number of observed years
                [y_mod, 1].

        var:    A flag that identifies if data are temperature or 
                precipitation. This flag tells the getDist function if it has 
                to discard distribution functions that allow negative numbers,
                and if the terms in the correction equations are
                multiplied/divided or added/subtracted.
                    Temperature:   var = 0
                    Precipitation: var = 1

        NOTE: This routine considers that obs and mod series start in january
        of the first year and ends in december of the last year.

    Optional inputs:
        frq:             A string specifying if the input is annual or monthly
                         data. If not specified, it is set monthly as default.
                             Monthly:   frq = 'M'
                             Annual:    frq = 'A'
                    
        pp_threshold:     A float indicating the threshold to consider 
                          physically null precipitation values.
                
        pp_factor:        A float indicating the maximum value of the random
                          values that replace physically null precipitation 
                          values.
                         
        rel_change_th:    A float indicating the maximum scaling factor 
                          (Equation 4 of Cannon et al. (2015)) when the
                          denominator is below inv_mod_th.
                
        inv_mod_th:       A float indicating the upper threshold of the 
                          denominator of the scaling factor (Equation 4 of
                          Cannon et al. (2015)) to truncate the scaling factor.
                          This parameter is defined as default as the
                          pp_threshold parameter, described above.
                          
    Output:
        QDM_series:  A column vector of monthly or annual modeled data 
                    (temperature or precipitation) corrected by the QDM method. 
                    If monthly frequency is specified, the length of this 
                    vector is 12 times the number of observed years 
                    [12 x y_mod, 1]. If annual frequency is specified, the 
                    length of this vector is equal to the number of observed 
                    years [y_mod, 1].
    
    References:
        Cannon, A. J., S. R. Sobie, and T. Q. Murdock, 2015: Bias correction of 
        GCM precipitation by quantile mapping: How well do methods preserve
        changes in quantiles and extremes? J. Climate, 28(17), 6938-6959, 
        https://doi.org/10.1175/JCLI-D-14-00754.1

    """
    
    if inv_mod_th is None:
        inv_mod_th = pp_threshold
    
    # 0) Format inputs and get statistics of the observed and modeled series of
    #    the historical period (formatQM function of the climQMBC package).
    y_obs, obs_series, mod_series, mu_obs, std_obs, skew_obs, skewy_obs, mu_mod, std_mod, skew_mod, skewy_mod = formatQM(obs, mod, var, frq, pp_threshold, pp_factor)
    
    # 1) Assign a probability distribution function to each month for the
    #    observed and modeled data in the historical period. If annual
    #    frequency is specified, this is applied to the complete historical
    #    period (getDist function of the climQMBC package).
    PDF_obs = getDist(obs_series, mu_obs, std_obs, skew_obs, skewy_obs, var)
    PDF_mod = getDist(mod_series[:,:y_obs], mu_mod, std_mod, skew_mod, skewy_mod, var)  

    # 2) For each projected period:
    PDF_win = np.zeros((mod_series.shape[0],mod_series.shape[1]-y_obs))
    Taot = np.zeros((mod_series.shape[0],mod_series.shape[1]-y_obs))
    for j in range(Taot.shape[1]):
        win_series = mod_series[:,j+1:y_obs+j+1]
        
        mux = np.nanmean(win_series, 1)
        sigmax = np.nanstd(win_series, 1, ddof=1)
        skewx = stat.skew(obs_series, 1, bias=False)
        with np.errstate(divide = 'ignore'):
            Ln_win  = np.log(win_series)
        Ln_win[np.imag(Ln_win)!=0] = 0
        Ln_win[np.isinf(Ln_win)] = np.log(0.01)
        skewy = stat.skew(Ln_win, 1, bias=False)
                
        # a) Assign a probability distribution function to each month. If
        #    annual frequency is specified, this is applied to the complete 
        #    period (getDist function of the climQMBC package).
        PDF_win[:,j] = getDist(win_series, mux, sigmax, skewx, skewy, var)
        
        # b) Apply the cumulative distribution function of the projected
        #    period, evaluated with the statistics of this period, to the last
        #    data of the period (getCDF function of the climQMBC package).
        #    Equation 3 of Cannon et al. (2015).
        Taot[:,j] = getCDF(PDF_win[:,j], mod_series[:,y_obs+j:y_obs+j+1], mux, sigmax, skewx, skewy)[:,0]
    
    # 3) Apply the inverse cumulative distribution function:
    #    a) Of the observed data, evaluated with the statistics of the observed
    #       data in the historical period, to the probabilities obtained from 
    #       2b) (getCDFinv function of the climQMBC package). Equation 5 of 
    #       Cannon et al. (2015).
    inv_obs = getCDFinv(PDF_obs, Taot, mu_obs, std_obs, skew_obs, skewy_obs)
    
    #    b) Of the modeled data, evaluated with the statistics of the observed
    #       data in the historical period, to the probabilities obtained from 
    #       2b) (getCDFinv function of the climQMBC package). Equations 4 of 
    #       Cannon et al. (2015).
    inv_mod = getCDFinv(PDF_mod, Taot, mu_mod, std_mod, skew_mod, skewy_mod)
    
    # 4) Get the delta factor or relative change and apply it to the value
    #    obtained in 3b). Equation 4 and 6 of Cannon et al. (2015).
    if var==1:
        DM = mod_series[:,y_obs:]/inv_mod
        DM[(DM>rel_change_th) & (inv_mod<inv_mod_th)] = rel_change_th
        QDM = inv_obs*DM
    else:
        DM = mod_series[:,y_obs:] - inv_mod
        QDM = inv_obs + DM
    
    QDM = QDM.reshape(-1,order='F')
    
    # 5) Perform QM for the historical period.
    mod_h = mod_series[:,:y_obs].reshape(-1, order='F')
    QM_series = QM(obs, mod_h, var, frq)
    QDM_series = np.hstack([QM_series, QDM])
    if var==1:
        QDM_series[QDM_series<pp_threshold] = 0
    
    return QDM_series


def UQM(obs, mod, var, frq='M', pp_threshold=1, pp_factor=1/100):
    """
    This function performs bias correction of modeled series based on observed
    data by the Unbiased Quantile Mapping (UQM) method, as described by 
    Chadwick et al. (2023). Correction is performed to monthly or annual
    precipitation or temperature data in a single location. An independent
    probability distribution function is assigned to each month and to each
    projected period based on the Kolmogorov-Smirnov test. If annual frequency
    is specified, this is applied to the complete period. Each projected period 
    considers a window whose length is equal to the number of years of the 
    historical period and ends in the analyzed year.
    
    Correction of the historical period is performed by the Quantile Mapping 
    (QM) method, as described by Cannon et al. (2015).
    
    Description:
        0) Format inputs and get statistics of the observed and modeled series
           of the historical period.
           
        1) Assign a probability distribution function to each month for the 
           observed and modeled data in the historical period. If annual 
           frequency is specified, this is applied to the complete historical 
           period.
           
        2) For each projected period, get the delta factor (delta) and time
           dependent (aster) statistics (mean, standard deviation, skewness, 
           and log skewness).

        3) For each projected period:
            a) Assign a probability distribution function to each month.
               If annual frequency is specified, this is applied to the 
               complete period.
            b) Apply the cumulative distribution function of the projected 
               period, evaluated with the statistics of this period, to the 
               last data of the period.

        4) Apply the inverse cumulative distribution function of the observed
           data, evaluated with the time dependent statistics, to the values
           obtained in 3b).
    
        End) Perform QM for the historical period.

    Inputs:
        obs:    A column vector of monthly or annual observed data (temperature
                or precipitation). If monthly frequency is specified, 
                the length of this vector is 12 times the number of observed 
                years [12 x y_obs, 1]. If annual frequency is specified, the 
                length of this vector is equal to the number of observed years
                [y_obs, 1].

        mod:    A column vector of monthly or annual modeled data (temperature
                or precipitation). If monthly frequency is specified, 
                the length of this vector is 12 times the number of observed
                years [12 x y_mod, 1]. If annual frequency is specified, the 
                length of this vector is equal to the number of observed years
                [y_mod, 1].

        var:    A flag that identifies if data are temperature or 
                precipitation. This flag tells the getDist function if it has 
                to discard distribution functions that allow negative numbers,
                and if the terms in the correction equations are
                multiplied/divided or added/subtracted.
                    Temperature:   var = 0
                    Precipitation: var = 1

        NOTE: This routine considers that obs and mod series start in january
        of the first year and ends in december of the last year.

    Optional inputs:
        frq:    A string specifying if the input is annual or monthly data. If
                not specified, it is set monthly as default.
                    Monthly:   frq = 'M'
                    Annual:    frq = 'A'
                    
        pp_threshold:    A float indicating the threshold to consider 
                         physically null precipitation values.
                
        pp_factor:       A float indicating the maximum value of the random
                         values that replace physically null precipitation 
                         values.

    Output:
        UQM_series:  A column vector of monthly or annual modeled data 
                    (temperature or precipitation) corrected by the UQM method. 
                    If monthly frequency is specified, the length of this 
                    vector is 12 times the number of observed years 
                    [12 x y_mod, 1]. If annual frequency is specified, the 
                    length of this vector is equal to the number of observed 
                    years [y_mod, 1].
    
    References:
        Chadwick, C., Gironás, J., González-Leiva, F., and Aedo, S. (2023). Bias 
        adjustment to preserve changes in variability: the unbiased mapping of GCM 
        changes. Hydrological Sciences Journal,
        https://doi.org/10.1080/02626667.2023.2201450
  
    """
    
    # 0) Format inputs and get statistics of the observed and modeled series of
    #    the historical period (formatQM function of the climQMBC package).
    y_obs, obs_series, mod_series, mu_obs, std_obs, skew_obs, skewy_obs, mu_mod, std_mod, skew_mod, skewy_mod = formatQM(obs, mod, var, frq, pp_threshold, pp_factor)
    
    # 1) Assign a probability distribution function to each month for the
    #    observed and modeled data in the historical period. If annual
    #    frequency is specified, this is applied to the complete historical
    #    period (getDist function of the climQMBC package).
    PDF_obs = getDist(obs_series, mu_obs, std_obs, skew_obs, skewy_obs, var)
    
    # 2) For each projected period, get the delta factor (delta) and time
    #    dependent (aster) statistics (mean, standard deviation, skewness, and
    #    log skewness). Equations 13 to 14 in Chadwick et al. (2023).
    xbarmt = np.zeros((mod_series.shape[0], mod_series.shape[1]-y_obs))
    xhatmt = np.zeros(xbarmt.shape)
    skwmt = np.zeros(xbarmt.shape)
    Lskwmt = np.zeros(xbarmt.shape)
    
    Dmu = np.zeros(xbarmt.shape)
    Dsigma = np.zeros(xbarmt.shape)
    Dskw = np.zeros(xbarmt.shape)
    DLskw = np.zeros(xbarmt.shape)
    
    muAster = np.zeros(xbarmt.shape)
    sigmaAster = np.zeros(xbarmt.shape)
    skwAster = np.zeros(xbarmt.shape)
    LskwAster = np.zeros(xbarmt.shape)
    
    if var==1:  # Precipitation
        for i in range(xbarmt.shape[0]):
            for j in range(xbarmt.shape[1]):
                xbarmt[i,j] = np.nanmean(mod_series[i,j+1:y_obs+j+1], 0)
                xhatmt[i,j] = np.nanstd(mod_series[i,j+1:y_obs+j+1], 0, ddof=1)
                
                skwmt[i,j] = stat.skew(mod_series[i,j+1:y_obs+j+1], 0, bias=False)
                Lnskwmt = np.log(mod_series[i,j+1:y_obs+j+1])
                Lnskwmt[np.isinf(Lnskwmt)]= np.log(0.01)
                Lskwmt[i,j] = stat.skew(Lnskwmt,0,bias=False)    
    
                Dmu[i,j] = (xbarmt[i,j]- mu_mod[i])/mu_mod[i]
                Dsigma[i,j] = (xhatmt[i,j]- std_mod[i])/std_mod[i]
                Dskw[i,j] = (skwmt[i,j] - skew_mod[i])/skew_mod[i]
                DLskw[i,j] = (Lskwmt[i,j]- skewy_mod[i])/skewy_mod[i]
                
                muAster[i,j] = mu_obs[i]*(1 + Dmu[i,j])
                sigmaAster[i,j] = std_obs[i]*(1 + Dsigma[i,j])
                skwAster[i,j] = skew_obs[i]*(1 + Dskw[i,j])
                LskwAster[i,j] = skewy_obs[i]*(1 + DLskw[i,j])
    else:  # Temperature
        for i in range(xbarmt.shape[0]):
            for j in range(xbarmt.shape[1]):
                xbarmt[i,j] = np.nanmean(mod_series[i,j+1:y_obs+j+1], 0)
                xhatmt[i,j] = np.nanstd(mod_series[i,j+1:y_obs+j+1], 0, ddof=1)
                
                skwmt[i,j] = stat.skew(mod_series[i,j+1:y_obs+j+1], 0, bias=False)
                Lnskwmt = np.log(mod_series[i,j+1:y_obs+j+1])
                Lnskwmt[np.isinf(Lnskwmt)] = np.log(0.01)
                Lskwmt[i,j] = stat.skew(Lnskwmt, 0, bias=False)    
    
                Dmu[i,j] = xbarmt[i,j] - mu_mod[i]
                Dsigma[i,j] = xhatmt[i,j] - std_mod[i]
                Dskw[i,j] = skwmt[i,j] - skew_mod[i]
                DLskw[i,j] = Lskwmt[i,j] - skewy_mod[i]
                
                muAster[i,j] = mu_obs[i] + Dmu[i,j]
                sigmaAster[i,j] = std_obs[i] + Dsigma[i,j]
                skwAster[i,j] = skew_obs[i] + Dskw[i,j]
                LskwAster[i,j] = skewy_obs[i] + DLskw[i,j]

    # 3) For each projected period:
    PDF_win = np.zeros((mod_series.shape[0], mod_series.shape[1]-y_obs))
    Taot = np.zeros((mod_series.shape[0], mod_series.shape[1]-y_obs))
    UQM = np.zeros((mod_series.shape[0], mod_series.shape[1]-y_obs))
    for j in range(Taot.shape[1]):
        win_series = mod_series[:,j+1:y_obs+j+1]
        
        mux  = np.nanmean(win_series, 1)
        sigmax = np.nanstd(win_series, 1, ddof=1)
        skewx = stat.skew(obs_series,1 , bias=False)
        with np.errstate(divide='ignore'):
            Ln_win  = np.log(win_series)
        Ln_win[np.imag(Ln_win)!=0] = 0
        Ln_win[np.isinf(Ln_win)] = np.log(0.01)
        skewy = stat.skew(Ln_win, 1, bias=False)
                
        # a) Assign a probability distribution function to each month. If
        #    annual frequency is specified, this is applied to the complete 
        #    period (getDist function of the climQMBC package).
        PDF_win[:,j] = getDist(win_series, mux, sigmax, skewx, skewy, var)
        
        # b) Apply the cumulative distribution function of the projected
        #    period, evaluated with the statistics of this period, to the last
        #    data of the period (getCDF function of the climQMBC package).
        #    Equation in Chadwick et al. (2023).
        Taot[:,j] = getCDF(PDF_win[:,j], mod_series[:,y_obs+j:y_obs+j+1], mux, sigmax, skewx, skewy)[:,0]
        
    # 4) Apply the inverse cumulative distribution function of the observed
    #    data, evaluated with the time dependent statistics, to the values
    #    obtained in 3b) (getCDFinv function of the climQMBC package). Equation
    #    15 in Chadwick et al. (2023).
    for yr in range(Taot.shape[1]):
        UQM[:,yr] = getCDFinv(PDF_obs, Taot[:,yr:yr+1], muAster[:,yr], sigmaAster[:,yr], skwAster[:,yr], LskwAster[:,yr])[:,0]
    
    UQM = UQM.reshape(-1, order='F')
    
    # 5) Perform QM for the historical period.
    mod_h = mod_series[:,:y_obs].reshape(-1, order='F')
    QM_series = QM(obs, mod_h, var, frq)
    UQM_series = np.hstack([QM_series, UQM])
    if var==1:
        UQM_series[UQM_series<pp_threshold] = 0
    
    return UQM_series


def SDM(obs, mod, var, frq='M', pp_threshold=1, pp_factor=1/100):
    """
    This function performs bias correction of modeled series based on observed
    data by the Scaled Distribution Mapping (SDM) method, as described by 
    Switanek et al. (2017). Correction is performed to monthly or annual 
    precipitation or temperature data in a single location. The historical 
    period of the modeled series is bias corrected as a scenario where the 
    complete series is replaced by the bias corrected series. On the other 
    hand, for each year after the historical period the method considers a 
    projected period as a window whose length is equal to the number of years 
    of the historical period and ends in the analyzed year. From the bias 
    corrected series, only the last year is saved and all the others are 
    discarded. This approach aims to build a continuum for the bias corrected     
    series.
    
    Description:
        0) Format inputs and get statistics of the observed and modeled series
           of the historical period.
           
        1) Historical period:
            a) [Switanek et al. (2017), step 1)]
                For temperature, get the detrended modeled and observed series in
                the historical period.
                For precipitation, get the rainday values and its frequency for the
                modeled and observed series in the historical period.
            b) [Switanek et al. (2017), step 2)]
                For temperature, fit normal probability distribution to the
                detrended observed and modeled series in the historical period.
                For precipitation, fit gamma distribution parameters by maximum 
                likelihood to the rainday values of the observed and modeled series
                in the historical period.
                Using these fitted distributions, get the corresponding cumulative
                distribution values. Tail probabilities are set to (0 + threshold)
                for temperature and to (1 - threshold) for temperature and
                precipitation. Threshold is set in the first lines of this function.
                Default is CDF_th = 10^-3.
        
        2) Projected periods:
            c) Initialize correction array.
            d) Define projected window.
            e) [Switanek et al. (2017), step 1)]
                For temperature, get the detrended series of the projected
                period.
                For precipitation, get the rainday values, its frequency, and 
                expected raindays for the projected period.
                Get the index of the sorted detrended or rainday values.
            f) [Switanek et al. (2017), step 2)]
                For temperature, fit normal probability distribution to the 
                detrended values of the projected period.
                For precipitation, fit gamma distribution parameters by maximum
                likelihood to the rainday values of the projected period.
                Using these fitted distributions, get the corresponding 
                cumulative distribution values. Tail probabilities are set to 
                (0 + threshold) for temperature and to (1 - threshold) for 
                temperature and precipitation. Threshold is set in the first 
                lines of this function. Default is CDF_th = 10^-3.
            g) [Switanek et al. (2017), step 3)]
                Get the scaling between the model projected period and 
                historical period distributions.
            h) Interpolate observed and modeled CDF of the historical period to
               the length of the projected period.
            i) [Switanek et al. (2017), steps 4 and 5)]
                Get recurrence intervals and its scaled version for the 
                observed, historical modeled and projected period modeled CDFs.
            j) [Switanek et al. (2017), step 6)]
                Get the initial bias corrected values. For precipitation, these
                values are interpolated to the length of the expected raindays.
            k) [Switanek et al. (2017), step 7)]
                Bias corrected values are placed back matching the higher bias
                corrected values with the higher rainday or detrended values.
                For temperature, the trend of the projected period is added
                back.
            l) If the projected period is the historical period (period 0, j=0)
                save the complete bias corrected series.
                If the projected period is not the historical period (j>0),
                save the value of the last year.
    
    Inputs:
        obs:    A column vector of monthly or annual observed data (temperature
                or precipitation). If monthly frequency is specified, 
                the length of this vector is 12 times the number of observed 
                years [12 x y_obs, 1]. If annual frequency is specified, the 
                length of this vector is equal to the number of observed years
                [y_obs, 1].

        mod:    A column vector of monthly or annual modeled data (temperature
                or precipitation). If monthly frequency is specified, 
                the length of this vector is 12 times the number of observed
                years [12 x y_mod, 1]. If annual frequency is specified, the 
                length of this vector is equal to the number of observed years
                [y_mod, 1].

        var:    A flag that identifies if data are temperature or 
                precipitation. This flag tells the getDist function if it has 
                to discard distribution functions that allow negative numbers,
                and if the terms in the correction equations are
                multiplied/divided or added/subtracted.
                    Temperature:   var = 0
                    Precipitation: var = 1

        NOTE: This routine considers that obs and mod series start in january
        of the first year and ends in december of the last year.

    Optional inputs:
        frq:    A string specifying if the input is annual or monthly data. If
                not specified, it is set monthly as default.
                    Monthly:   frq = 'M'
                    Annual:    frq = 'A'
                    
        pp_threshold:    A float indicating the threshold to consider 
                         physically null precipitation values.
                
        pp_factor:       A float indicating the maximum value of the random
                         values that replace physically null precipitation 
                         values.

    Output:
        SDM_series:  A column vector of monthly or annual modeled data 
                    (temperature or precipitation) corrected by the SDM method. 
                    If monthly frequency is specified, the length of this 
                    vector is 12 times the number of observed years 
                    [12 x y_mod, 1]. If annual frequency is specified, the 
                    length of this vector is equal to the number of observed 
                    years [y_mod, 1].
    
    References:
        Switanek, B. M., Troch, P., Castro, L. C., Leuprecht, A., Chang, H. I.,
        Mukherjee, R., and Demaria, M. C. E. (2017) Scaled distribution
        mapping: A bias correction method that preserves raw climate model
        projected changes. Hydrology &amp; Earth System Sciences, 21,
        2649-2666, https://doi.org/10.5194/hess-21-2649-2017.
 
    """
    
    lower_lim = pp_threshold
    CDF_th = 1e-3
    
    # 0) Format inputs and get statistics of the observed and modeled series of
    #    the historical period (formatQM function of the climQMBC package).
    y_obs,obs_series,mod_series = formatQM(obs, mod, var, frq, pp_threshold, pp_factor)[:3]
    
    SDM = np.zeros((mod_series.shape[0], mod_series.shape[1]-y_obs))    
    SDM_h = np.zeros((obs_series.shape[0], y_obs))    
    for m in range(mod_series.shape[0]):
        # 1) Historical period:

        # a) [Switanek et al. (2017), step 1)]
        #     For temperature, get the detrended modeled and observed series in
        #     the historical period.
        #     For precipitation, get the rainday values and its frequency for
        #     the modeled and observed series in the historical period.
         
        if var == 0: # Temperature
            D_obs = detrend(obs_series[m])
            D_mod = detrend(mod_series[m,:y_obs])
            
            mu_obs = np.nanmean(obs_series[m])
            mu_mod = np.nanmean(mod_series[m,:y_obs])
        else:       # Precipitation
            D_obs = np.sort(obs_series[m][obs_series[m]>lower_lim])
            D_mod = np.sort(mod_series[m,:y_obs][mod_series[m,:y_obs]>lower_lim])
            
            freq_obs = D_obs.shape[0]/obs_series.shape[1]
            freq_mod = D_mod.shape[0]/mod_series[m,:y_obs].shape[0]
    
        # b) [Switanek et al. (2017), step 2)]
        #    For temperature, fit normal probability distribution to the
        #    detrended observed and modeled series in the historical period.
        #    For precipitation, fit gamma distribution parameters by maximum
        #    likelihood to the rainday values of the observed and modeled
        #    series in the historical period.
        #    Using these fitted distributions, get the corresponding cumulative
        #    distribution values. Tail probabilities are set to (0 + threshold)
        #    for temperature and to (1 - threshold) for temperature and
        #    precipitation. Threshold is set in the first lines of this
        #    function. Default is CDF_th = 10^-3.
        if var == 0: # Temperature
            mu_obsD = np.nanmean(D_obs)
            mu_modD = np.nanmean(D_mod)
            
            sigma_obsD = np.nanstd(D_obs, 0, ddof=0)
            sigma_modD = np.nanstd(D_mod, 0, ddof=0)
            
            CDF_obs = stat.norm.cdf(np.sort(D_obs), *(mu_obsD,sigma_obsD))
            CDF_mod = stat.norm.cdf(np.sort(D_mod), *(mu_modD,sigma_modD))
            
            CDF_obs[CDF_obs<CDF_th] = CDF_th
            CDF_mod[CDF_mod<CDF_th] = CDF_th
            CDF_obs[CDF_obs>1-CDF_th] = 1 - CDF_th
            CDF_mod[CDF_mod>1-CDF_th] = 1 - CDF_th
            
        else: # Precipitation
            fit_obs = stat.gamma.fit(D_obs, floc=0)
            fit_mod = stat.gamma.fit(D_mod, floc=0)
            
            CDF_obs = stat.gamma.cdf(D_obs, *fit_obs)
            CDF_mod = stat.gamma.cdf(D_mod, *fit_mod)
            CDF_obs[CDF_obs>1-CDF_th] = 1 - CDF_th
            CDF_mod[CDF_mod>1-CDF_th] = 1 - CDF_th
            
        # 2) Projected periods:
        for j in range(-1, mod_series.shape[1]-y_obs):
            # c) Initialize correction array.
            corr_temp = np.zeros(y_obs)
        
            # d) Define projected window.
            win_series  = mod_series[m,j+1:y_obs+j+1]
            
            # e) [Switanek et al. (2017), step 1)]
            #    For temperature, get the detrended series of the projected
            #    period.
            #    For precipitation, get the rainday values, its frequency, and
            #    expected raindays for the projected period.
            #    Get the index of the sorted detrended or rainday values.
            if var==0: # Temperature
                D_win = detrend(win_series)
                exp_D = win_series.shape[0]
                win_argsort = np.argsort(D_win)
            
                mu_win = np.nanmean(win_series)
            else: # Precipitation
                D_win = np.sort(win_series[win_series>lower_lim])
                exp_D = min(round(D_win.shape[0]*freq_obs/freq_mod), win_series.shape[0])
                win_argsort = np.argsort(win_series)
    
            # f) [Switanek et al. (2017), step 2)]
            #    For temperature, fit normal probability distribution to the
            #    detrended values of the projected period.
            #    For precipitation, fit gamma distribution parameters by
            #    maximum likelihood to the rainday values of the projected
            #    period.
            #    Using these fitted distributions, get the corresponding
            #    cumulative distribution values. Tail probabilities are set to
            #    (0 + threshold) for temperature and to (1 - threshold) for
            #    temperature and precipitation. Threshold is set in the first
            #    lines of this function. Default is CDF_th = 10^-3.
            if var==0: # Temperature
                mu_winD = np.nanmean(D_win)
                sigma_winD = np.nanstd(D_win, 0, ddof=0)
                
                diff_win = win_series - D_win
                
                CDF_win = stat.norm.cdf(np.sort(D_win), *(mu_winD,sigma_winD))
                CDF_win[CDF_win<CDF_th] = CDF_th
                CDF_win[CDF_win>1-CDF_th] = 1 - CDF_th
            
            else: #Precipitation
                fit_win = stat.gamma.fit(D_win, floc=0)
                CDF_win = stat.gamma.cdf(D_win, *fit_win)    
                CDF_win[CDF_win>1-CDF_th] = 1 - CDF_th

            # g) [Switanek et al. (2017), step 3)]
            #    Get the scaling between the model projected period and
            #    historical period distributions.
            if var==0: # Temperature
                SF = (sigma_obsD/sigma_modD)*(stat.norm.ppf(CDF_win, *(mu_winD,sigma_winD)) - stat.norm.ppf(CDF_win, *(mu_modD,sigma_modD)))
            else: # Precipitation
                SF = stat.gamma.ppf(CDF_win, *fit_win)/stat.gamma.ppf(CDF_win, *fit_mod)
    
            # h) Interpolate observed and modeled CDF of the historical period
            #    to the length of the projected period.
            obs_cdf_intpol = np.interp(np.linspace(1, len(D_obs), len(D_win)),
                                       np.linspace(1, len(D_obs), len(D_obs)),
                                       CDF_obs)
            mod_cdf_intpol = np.interp(np.linspace(1, len(D_mod), len(D_win)),
                                       np.linspace(1, len(D_mod), len(D_mod)),
                                       CDF_mod)

            # i) [Switanek et al. (2017), steps 4 and 5)]
            #    Get recurrence intervals and its scaled version for the
            #    observed, historical modeled and projected period modeled
            #    CDFs.
            if var == 0: # Temperature
                RI_obs = 1/(0.5 - np.abs(obs_cdf_intpol - 0.5))
                RI_mod = 1/(0.5 - np.abs(mod_cdf_intpol - 0.5))
                RI_win = 1/(0.5 - np.abs(CDF_win - 0.5))
                RI_scaled = RI_obs*RI_win/RI_mod
                
                CDF_scaled = np.sign(obs_cdf_intpol - 0.5)*(1 - 1/RI_scaled)
                CDF_scaled[CDF_scaled<0] = CDF_scaled[CDF_scaled<0] + 1
                
                CDF_scaled[CDF_scaled<CDF_th] = CDF_th
                CDF_scaled[CDF_scaled>1-CDF_th] = 1 - CDF_th
            else: # Precipitation
                RI_obs = 1/(1 - obs_cdf_intpol)
                RI_mod = 1/(1 - mod_cdf_intpol)
                RI_win = 1/(1 - CDF_win)
                RI_scaled = RI_obs*RI_win/RI_mod
                RI_scaled[RI_scaled<1] = 1
                
                CDF_scaled = np.sort(1 - 1/RI_scaled)
                
            # j) [Switanek et al. (2017), step 6)]
            #    Get the initial bias corrected values. For precipitation,
            #    these values are interpolated to the length of the expected
            #    raindays.
            if var == 0: # Temperature
                xvals = stat.norm.ppf(np.sort(CDF_scaled), *(mu_obsD,sigma_obsD)) + SF
                xvals = xvals - np.mean(xvals) + mu_obs + (mu_win - mu_mod)
            else:
                xvals = stat.gamma.ppf(CDF_scaled, *fit_obs)*SF
                if D_win.shape[0] > exp_D:
                    xvals = np.interp(np.linspace(1, len(D_win), exp_D),
                                      np.linspace(1, len(D_win), len(D_win)),
                                      xvals)
                else:
                    xvals = np.hstack([np.zeros(exp_D-D_win.shape[0]), xvals])
    
            # k) [Switanek et al. (2017), step 7)]
            #    Bias corrected values are placed back matching the higher bias
            #    corrected values with the higher rainday or detrended values.
            #    For temperature, the trend of the projected period is added
            #    back.
            corr_temp[win_argsort[-exp_D:]] = xvals
            if var==0:
                corr_temp = corr_temp + diff_win - mu_win
                
            # l) If the projected period is the historical period (period 0,
            #    j=0) save the complete bias corrected series.
            #    If the projected period is not the historical period (j>0),
            #    save the value of the last year.
            if j==-1:
                SDM_h[m] = corr_temp
            else:
                SDM[m,j] = corr_temp[-1]
            
    SDM_h = SDM_h.reshape(-1, order='F')
    SDM = SDM.reshape(-1, order='F')
    SDM_series = np.hstack([SDM_h, SDM])
    if var==1:
        SDM_series[SDM_series<pp_threshold] = 0

    return SDM_series