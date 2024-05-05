# -*- coding: utf-8 -*-
"""
This script contains the bias correction methods implemented in the climQMBC
package at daily, monthly and annual scale, including the Quantile Mapping (QM),
the Detrended Quantile Mapping (DQM), the Quantile Delta Mapping (QDM), the
Unbiased Quantile Mapping (UQM), and the Scaled Distribution Mapping (SDM).

References:
    Cannon, A. J., S. R. Sobie, and T. Q. Murdock. (2015). Bias correction of 
    GCM precipitation by quantile mapping: How well do methods preserve changes
    in quantiles and extremes? J. Climate, 28(17), 6938-6959, 
    https://doi.org/10.1175/JCLI-D-14-00754.1
    
    Switanek, B. M., Troch, P., Castro, L. C., Leuprecht, A., Chang, H. I., 
    Mukherjee, R., and Demaria, M. C. E. (2017). Scaled distribution mapping: A 
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
Revision: 2, updated Apr 2024
"""

from .utils import formatQM, getStats, getDist, getCDF, getCDFinv
from .utils import projected_backward_moving_window, day_centered_moving_window
from .utils import get_pp_threshold_mod, set_norain_to_nan
from scipy.signal import detrend
import scipy.stats as stat
import numpy as np


def QM(obs, mod, allow_negatives=1, frq='A', pp_threshold=1, pp_factor=1/100,
       day_win=1, user_pdf=False, pdf_obs=None, pdf_mod=None):
    """
    This function performs bias correction of modeled series based on observed
    data by the Quantile Mapping (QM) method, as described by Cannon et al. 
    (2015). Correction is performed to daily, monthly or annual data in a single
    location. An independent probability distribution function is assigned to 
    each sub-annual period based on the Kolmogorov-Smirnov test.

    Inputs:
        obs:    A column vector of daily, monthly or annual observed data,
                without considering leap days. If daily frequency is specified,
                the length of the column vector should by a multiple of 365 and
                for monthly frequency, it should be a multiple of 12.
                [ndata_obs, 1]

        mod:    A column vector of daily, monthly or annual modeled or GCM data,
                without considering leap days. If daily frequency is specified,
                the length of the column vector should by a multiple of 365 and
                for monthly frequency, it should be a multiple of 12.
                [ndata_mod, 1]

        NOTE: This routine considers that obs and mod series start in the same
        day/month/year and are continuous until the end day/month/year.

    Optional inputs:
        allow_negatives: A flag that identifies if data allows negative values
                         and also to replace no-rain values with random small 
                         values (Chadwick et al., 2023) to avoid numerical
                         problems with the probability distribution functions.
                             allow_negatives = 1 or True: Allow negatives (default)
                             allow_negatives = 0 or False: Do not allow negative

        frq:             A string specifying if the input frequency is daily,
                         monthly or annual.
                             Daily:     frq = 'D'
                             Monthly:   frq = 'M'
                             Annual:    frq = 'A' (default)

        pp_threshold:    A float indicating the threshold to consider no-rain
                         values. Default is 1.
                
        pp_factor:       A float which multiplied to pp_threshold indicates the
                         maximum value of the random values that replace no-rain
                         values. Default is 1/100.

        day_win:         (only for frq='D') An integer indicating how many days
                         to consider backwards and forward to get the statistics
                         of each calendar day.The length of the window will be 
                         (2*win_day-1). For example, day_win=15 -> window of 29.
                         Default: win = 1

        user_pdf:        A flag indicating if the user will define the
                         probability distribution functions (pdf) for the
                         observed and modeled series. The distributions will be
                         the same for all periods and sub-periods.
                            user_pdf = 1 or True: User defines the pdf
                            user_pdf = 0 or False: pdf defined by the Kolmogorov
                                                   -Smirnov test (default)

        pdf_obs:         An integer indicating the probability distribution 
                         function (pdf) to be used for the observed data. The
                         pdf will be the same for all periods and sub-periods.
                         Default: None

        pdf_mod:         An integer indicating the probability distribution 
                         function (pdf) to be used for the modeled data. The
                         pdf will be the same for all periods and sub-periods.
                         Default: None
                         
        NOTE: The available distributions for pdf_obs and pdf_mod are:
            0) Normal
            1) Log-Normal
            2) Gamma 2 parameters
            3) Gamma 3 parameters
               (Pearson 3 parameters)
            4) Log-Gamma 3 parameters
               (Log-Pearson 3 parameters)
            5) Gumbel
            6) Exponential

    Output:
        QM_series: A column vector of data bias corrected with the QM method.
                   [ndata_mod, 1]

    References:
        Cannon, A. J., S. R. Sobie, and T. Q. Murdock. (2015). Bias correction
        of GCM precipitation by quantile mapping: How well do methods preserve
        changes in quantiles and extremes? J. Climate, 28(17), 6938-6959, 
        https://doi.org/10.1175/JCLI-D-14-00754.1
    """

    # 0) Only for daily frequency data that does not allow negatives (as
    #    precipitacion). Get a pp_threshold for the modeled series that result
    #    in the same amount of rain-day values in the historical period as in
    #    the observed series
    if (frq=='D') & (not allow_negatives):
        pp_threshold_mod = get_pp_threshold_mod(obs, mod, pp_threshold)
    else:
        pp_threshold_mod = pp_threshold

    # 1) Format series to matrix with rows as sub-periods and columns as years
    #    and, if needed, replace no-rain values with random small values 
    y_obs, obs_series = formatQM(obs, allow_negatives, frq, pp_threshold, pp_factor)
    y_mod, mod_series = formatQM(mod, allow_negatives, frq, pp_threshold_mod, pp_factor)
    
    # 2) Get statistics of the observed and modeled series in the historical
    #    period.
    if frq=='D':
        # Get a 3D array with a moving window centered on each day
        obs_series_moving = day_centered_moving_window(obs_series, day_win)
        mod_series_moving = day_centered_moving_window(mod_series, day_win)
        
        # Reshape to 2D in order to have all the days of the window in each row
        obs_series_moving = obs_series_moving.reshape(365,(2*day_win-1)*y_obs)
        modh_series_moving = mod_series_moving[:,:,:y_obs].reshape(365,(2*day_win-1)*y_obs)
        
        # Replace no-rain values with nans but keep at least 30 values to get
        # the statistics. If fewer than 30 values are rain values, the missing
        # ones will be previously random values replaced
        if not allow_negatives:
            obs_series_moving = set_norain_to_nan(obs_series_moving, pp_threshold, pp_factor)
            modh_series_moving = set_norain_to_nan(modh_series_moving, pp_threshold_mod, pp_factor)
            mod_series[mod_series<pp_threshold_mod] = np.nan
        
        # Get the statistics for each row
        mu_obs, std_obs, skew_obs, skewy_obs = getStats(obs_series_moving)
        mu_mod, std_mod, skew_mod, skewy_mod = getStats(modh_series_moving)
        
    else:
        # Get the statistics for each row
        mu_obs, std_obs, skew_obs, skewy_obs = getStats(obs_series)
        mu_mod, std_mod, skew_mod, skewy_mod = getStats(mod_series[:,:y_obs])
    
    # 3) Assign a probability distribution function to each sub-period for the
    #    observed and modeled data in the historical period.
    if user_pdf==False:
        if frq=='D':
            pdf_obs, ks_fail_obs = getDist(obs_series_moving, allow_negatives, mu_obs, std_obs, skew_obs, skewy_obs)
            pdf_mod, ks_fail_mod = getDist(modh_series_moving, allow_negatives, mu_mod, std_mod, skew_mod, skewy_mod)
        else:
            pdf_obs, ks_fail_obs = getDist(obs_series, allow_negatives, mu_obs, std_obs, skew_obs, skewy_obs)
            pdf_mod, ks_fail_mod = getDist(mod_series[:,:y_obs], allow_negatives, mu_mod, std_mod, skew_mod, skewy_mod)
    else:
        pdf_obs = np.zeros(obs_series.shape[0]) + pdf_obs
        pdf_mod = np.zeros(mod_series.shape[0]) + pdf_mod

    # 4) Apply the cumulative distribution function of the modeled data,
    #    evaluated with the statistics of the modeled data in the historical
    #    period and the modeled data. Equation 1 of Cannon et al. (2015).
    prob = getCDF(pdf_mod, mod_series, mu_mod, std_mod, skew_mod, skewy_mod)
    
    # 5) Apply the inverse cumulative distribution function of the observed
    #    data, evaluated with the statistics of the observed data in the
    #    historical period, to the probabilities obtained from 2). Equation 1 of
    #    Cannon et al. (2015).
    QM_series = getCDFinv(pdf_obs, prob, mu_obs, std_obs, skew_obs, skewy_obs)
    
    # 6) Reshape to a column vector
    QM_series = QM_series.reshape(-1, order='F')
    
    # 7) If does not allow negative values, replace nans (should be the result 
    #    of correcting removed no-rain values) and replace no-rain values with 0
    if not allow_negatives:
        QM_series[np.isnan(QM_series)] = 0
        QM_series[QM_series<pp_threshold] = 0
    
    if user_pdf==False:
        # Check if ks-test failures
        ks_fail = ks_fail_obs + ks_fail_mod
        if ks_fail>0:
            print('QM: Some of the probability distribution functions did not pass the KS-Test')

    return QM_series


def DQM(obs, mod, mult_change=1, allow_negatives=1, frq='A', pp_threshold=1,
        pp_factor=1/100, day_win=1, user_pdf=False, pdf_obs=None, pdf_mod=None):
    """
    This function performs bias correction of modeled series based on observed
    data by the Detrended Quantile Mapping (DQM) method, as described by Cannon
    et al. (2015). Correction is performed to daily, monthly or annual data in a
    single location. An independent probability distribution function is 
    assigned to each sub-annual period based on the Kolmogorov-Smirnov test.
    Each projected period considers a window whose length is equal to the number
    of years of the historical period and ends in the analyzed year.
    
    Correction of the historical period is performed by the Quantile Mapping 
    (QM) method, as described in Cannon et al. (2015) and Chadwick et al. (2023).

    Inputs:
        obs:    A column vector of daily, monthly or annual observed data,
                without considering leap days. If daily frequency is specified,
                the length of the column vector should by a multiple of 365 and
                for monthly frequency, it should be a multiple of 12.
                [ndata_obs, 1]

        mod:    A column vector of daily, monthly or annual modeled or GCM data,
                without considering leap days. If daily frequency is specified,
                the length of the column vector should by a multiple of 365 and
                for monthly frequency, it should be a multiple of 12.
                [ndata_mod, 1]

        NOTE: This routine considers that obs and mod series start in the same
        day/month/year and are continuous until the end day/month/year.

    Optional inputs:
        mult_change:     A flag that indicates if projected changes should be
                         computed as multiplicative (fut = hist*delta) or 
                         additive (fut = hist + delta) changes.
                             mult_change = 1 or True: Multiplicative (default)
                             mult_change = 0 or False: Additive

        allow_negatives: A flag that identifies if data allows negative values
                         and also to replace no-rain values with random small 
                         values (Chadwick et al., 2023) to avoid numerical
                         problems with the probability distribution functions.
                             allow_negatives = 1 or True: Allow negatives (default)
                             allow_negatives = 0 or False: Do not allow negative

        frq:             A string specifying if the input frequency is daily,
                         monthly or annual.
                             Daily:     frq = 'D'
                             Monthly:   frq = 'M'
                             Annual:    frq = 'A' (default)

        pp_threshold:    A float indicating the threshold to consider no-rain
                         values. Default is 1.
                
        pp_factor:       A float which multiplied to pp_threshold indicates the
                         maximum value of the random values that replace no-rain
                         values. Default is 1/100.

        day_win:         (only for frq='D') An integer indicating how many days
                         to consider backwards and forward to get the statistics
                         of each calendar day.The length of the window will be 
                         (2*win_day-1). For example, day_win=15 -> window of 29.
                         Default: win = 1

        user_pdf:        A flag indicating if the user will define the
                         probability distribution functions (pdf) for the
                         observed and modeled series. The distributions will be
                         the same for all periods and sub-periods.
                            user_pdf = 1 or True: User defines the pdf
                            user_pdf = 0 or False: pdf defined by the Kolmogorov
                                                   -Smirnov test (default)

        pdf_obs:         An integer indicating the probability distribution 
                         function (pdf) to be used for the observed data. The
                         pdf will be the same for all periods and sub-periods.
                         Default: None

        pdf_mod:         An integer indicating the probability distribution 
                         function (pdf) to be used for the modeled data. The
                         pdf will be the same for all periods and sub-periods.
                         Default: None

        NOTE: The available distributions for pdf_obs and pdf_mod are:
            0) Normal
            1) Log-Normal
            2) Gamma 2 parameters
            3) Gamma 3 parameters
               (Pearson 3 parameters)
            4) Log-Gamma 3 parameters
               (Log-Pearson 3 parameters)
            5) Gumbel
            6) Exponential

    Output:
        DQM_series: A column vector of data bias corrected with the DQM method.
                    [ndata_mod, 1]
    
    References:
        Cannon, A. J., S. R. Sobie, and T. Q. Murdock. (2015). Bias correction
        of GCM precipitation by quantile mapping: How well do methods preserve
        changes in quantiles and extremes? J. Climate, 28(17), 6938-6959, 
        https://doi.org/10.1175/JCLI-D-14-00754.1
    """

    # 0) Only for daily frequency data that does not allow negatives (as
    #    precipitacion). Get a pp_threshold for the modeled series that result
    #    in the same amount of rain-day values in the historical period as in
    #    the observed series
    if (frq=='D') & (not allow_negatives):
        pp_threshold_mod = get_pp_threshold_mod(obs, mod, pp_threshold)
    else:
        pp_threshold_mod = pp_threshold
    
    # 1) Format series to matrix with rows as sub-periods and columns as years
    #    and, if needed, replace no-rain values with random small values 
    y_obs, obs_series = formatQM(obs, allow_negatives, frq, pp_threshold, pp_factor)
    y_mod, mod_series = formatQM(mod, allow_negatives, frq, pp_threshold_mod, pp_factor)
    
    # 2) Get statistics of the observed and modeled series in the historical
    #    period.
    if frq=='D':
        # Get a 3D array with a moving window centered on each day
        obs_series_moving = day_centered_moving_window(obs_series, day_win)
        mod_series_moving = day_centered_moving_window(mod_series, day_win)
        
        # Get a 4D array with a backward moving window for each projected period
        win_series = projected_backward_moving_window(mod_series_moving, y_obs, frq)

        # Reshape to 2D in order to have all the days of the window in each row
        obs_series_moving = obs_series_moving.reshape(365,(2*day_win-1)*y_obs)
        modh_series_moving = mod_series_moving[:,:,:y_obs].reshape(365,(2*day_win-1)*y_obs)
        
        # Reshape to 3D in order to have all the days of the window in each row
        win_series_moving = win_series.reshape(365,(2*day_win-1)*y_obs,y_mod-y_obs)
        
        # Replace no-rain values with nans but keep at least 30 values to get
        # the statistics. If fewer than 30 values are rain values, the missing
        # ones will be previously random values replaced
        if not allow_negatives:
            obs_series_moving = set_norain_to_nan(obs_series_moving, pp_threshold, pp_factor)
            modh_series_moving = set_norain_to_nan(modh_series_moving, pp_threshold_mod, pp_factor)
            mod_series[mod_series<pp_threshold_mod] = np.nan
            win_series_moving = set_norain_to_nan(win_series_moving, pp_threshold_mod, pp_factor)

        # Get the statistics for each row
        mu_obs, std_obs, skew_obs, skewy_obs = getStats(obs_series_moving)
        mu_mod, std_mod, skew_mod, skewy_mod = getStats(modh_series_moving)
        mu_win, std_win, skew_win, skewy_win = getStats(win_series_moving)
        
    else:
        # Get a 3D array with a backward moving window for each projected period
        win_series = projected_backward_moving_window(mod_series, y_obs, frq)
        
        # Get the statistics for each row
        mu_obs, std_obs, skew_obs, skewy_obs = getStats(obs_series)
        mu_mod, std_mod, skew_mod, skewy_mod = getStats(mod_series[:,:y_obs])
        mu_win, std_win, skew_win, skewy_win = getStats(win_series)
    
    # 3) Assign a probability distribution function to each sub-period for the
    #    observed and modeled data in the historical period.
    if user_pdf==False:
        if frq=='D':
            pdf_obs, ks_fail_obs = getDist(obs_series_moving, allow_negatives, mu_obs, std_obs, skew_obs, skewy_obs)
            pdf_mod, ks_fail_mod = getDist(modh_series_moving, allow_negatives, mu_mod, std_mod, skew_mod, skewy_mod)
        else:
            pdf_obs, ks_fail_obs = getDist(obs_series, allow_negatives, mu_obs, std_obs, skew_obs, skewy_obs)
            pdf_mod, ks_fail_mod = getDist(mod_series[:,:y_obs], allow_negatives, mu_mod, std_mod, skew_mod, skewy_mod)
    else:
        pdf_obs = np.zeros(obs_series.shape[0]) + pdf_obs
        pdf_mod = np.zeros(mod_series.shape[0]) + pdf_mod

    # 4) Compute the linear scaled values (value in square brackets in equation
    #    2 of Cannon et al. (2015)).
    mu_mod_repeated = np.tile(mu_mod, (y_mod-y_obs, 1)).T
    if mult_change:
        scaling_factor = mu_win/mu_mod_repeated
        detrended_series  = mod_series[:,y_obs:]/scaling_factor
    else:
        scaling_factor = mu_win - mu_mod_repeated
        detrended_series  = mod_series[:,y_obs:] - scaling_factor

    # 5) Apply the cumulative distribution function of the modeled data,
    #    evaluated with the statistics of the modeled data in the historical
    #    period and the modeled data. Equation 2 of Cannon et al. (2015).
    prob = getCDF(pdf_mod, detrended_series, mu_mod, std_mod, skew_mod, skewy_mod)
    
    # 6) Apply the inverse cumulative distribution function of the observed
    #    data, evaluated with the statistics of the observed data in the
    #    historical period, to the probabilities obtained from 5). Equation 2 of
    #    Cannon et al. (2015).
    DQM_raw = getCDFinv(pdf_obs, prob, mu_obs, std_obs, skew_obs, skewy_obs)
    
    # 7) Reimpose the trend to the values obtained in 6). Equation 2 of Cannon
    #    et al. (2015).
    if mult_change:
        DQM_series = DQM_raw*scaling_factor
    else:
        DQM_series = DQM_raw + scaling_factor
    
    # 8) Perform QM for the historical period.
    mod_h = mod[:obs.shape[0]]
    QM_series = QM(obs, mod_h, allow_negatives, frq,pp_threshold, pp_factor, day_win, user_pdf, pdf_obs, pdf_mod)
    
    # 9) Reshape to a column vector
    DQM_series = DQM_series.reshape(-1, order='F')
    
    # 10) Concat QM (hsitorical) and DQM (future)
    DQM_series = np.hstack([QM_series, DQM_series])
    
    # 11) If does not allow negative values, replace nans (should be the result
    #     of correcting removed no-rain values) and replace no-rain values with 0
    if not allow_negatives:
        DQM_series[np.isnan(DQM_series)] = 0
        DQM_series[DQM_series<pp_threshold] = 0
    
    if user_pdf==False:
        # Check if ks-test failures
        ks_fail = ks_fail_obs + ks_fail_mod
        if ks_fail>0:
            print('DQM: Some of the probability distribution functions did not pass the KS-Test')
    
    return DQM_series


def QDM(obs, mod, mult_change=1, allow_negatives=1, frq='A', pp_threshold=1,
        pp_factor=1/100, rel_change_th=2, inv_mod_th=None, day_win=1,
        user_pdf=False, pdf_obs=None, pdf_mod=None):
    """
    This function performs bias correction of modeled series based on observed
    data by the Quantile Delta Mapping (QDM) method, as described by Cannon
    et al. (2015). Correction is performed to daily, monthly or annual data in a
    single location. An independent probability distribution function is 
    assigned to each sub-annual period based on the Kolmogorov-Smirnov test.
    Each projected period considers a window whose length is equal to the number
    of years of the historical period and ends in the analyzed year.
    
    Correction of the historical period is performed by the Quantile Mapping 
    (QM) method, as described in Cannon et al. (2015) and Chadwick et al. (2023).
    
    Inputs:
        obs:    A column vector of daily, monthly or annual observed data,
                without considering leap days. If daily frequency is specified,
                the length of the column vector should by a multiple of 365 and
                for monthly frequency, it should be a multiple of 12.
                [ndata_obs, 1]

        mod:    A column vector of daily, monthly or annual modeled or GCM data,
                without considering leap days. If daily frequency is specified,
                the length of the column vector should by a multiple of 365 and
                for monthly frequency, it should be a multiple of 12.
                [ndata_mod, 1]

        NOTE: This routine considers that obs and mod series start in the same
        day/month/year and are continuous until the end day/month/year.

    Optional inputs:
        mult_change:     A flag that indicates if projected changes should be
                         computed as multiplicative (fut = hist*delta) or 
                         additive (fut = hist + delta) changes.
                             mult_change = 1 or True: Multiplicative (default)
                             mult_change = 0 or False: Additive

        allow_negatives: A flag that identifies if data allows negative values
                         and also to replace no-rain values with random small 
                         values (Chadwick et al., 2023) to avoid numerical
                         problems with the probability distribution functions.
                             allow_negatives = 1 or True: Allow negatives (default)
                             allow_negatives = 0 or False: Do not allow negative

        frq:             A string specifying if the input frequency is daily,
                         monthly or annual.
                             Daily:     frq = 'D'
                             Monthly:   frq = 'M'
                             Annual:    frq = 'A' (default)

        pp_threshold:    A float indicating the threshold to consider no-rain
                         values. Default is 1.
                
        pp_factor:       A float which multiplied to pp_threshold indicates the
                         maximum value of the random values that replace no-rain
                         values. Default is 1/100.

        rel_change_th:    A float indicating the maximum scaling factor 
                          (Equation 4 of Cannon et al. (2015)) when the
                          denominator is below inv_mod_th.

        inv_mod_th:      A float indicating the upper threshold of the 
                         denominator of the scaling factor (Equation 4 of
                         Cannon et al. (2015)) to truncate the scaling factor.
                         This parameter is defined as default as the
                         pp_threshold parameter, described above.

        day_win:         (only for frq='D') An integer indicating how many days
                         to consider backwards and forward to get the statistics
                         of each calendar day.The length of the window will be 
                         (2*win_day-1). For example, day_win=15 -> window of 29.
                         Default: win = 1

        user_pdf:        A flag indicating if the user will define the
                         probability distribution functions (pdf) for the
                         observed and modeled series. The distributions will be
                         the same for all periods and sub-periods.
                            user_pdf = 1 or True: User defines the pdf
                            user_pdf = 0 or False: pdf defined by the Kolmogorov
                                                   -Smirnov test (default)

        pdf_obs:         An integer indicating the probability distribution 
                         function (pdf) to be used for the observed data. The
                         pdf will be the same for all periods and sub-periods.
                         Default: None

        pdf_mod:         An integer indicating the probability distribution 
                         function (pdf) to be used for the modeled data. The
                         pdf will be the same for all periods and sub-periods.
                         Default: None

        NOTE: The available distributions for pdf_obs and pdf_mod are:
            0) Normal
            1) Log-Normal
            2) Gamma 2 parameters
            3) Gamma 3 parameters
               (Pearson 3 parameters)
            4) Log-Gamma 3 parameters
               (Log-Pearson 3 parameters)
            5) Gumbel
            6) Exponential

    Output:
        QDM_series: A column vector of data bias corrected with the QDM method.
                   [ndata_mod, 1]
    
    References:
        Cannon, A. J., S. R. Sobie, and T. Q. Murdock. (2015). Bias correction
        of GCM precipitation by quantile mapping: How well do methods preserve
        changes in quantiles and extremes? J. Climate, 28(17), 6938-6959, 
        https://doi.org/10.1175/JCLI-D-14-00754.1
    """
    
    # Sets the threshold of the modeled series to compute de delta of the
    # quantile as pp_threshold.
    if inv_mod_th is None:
        inv_mod_th = pp_threshold
    
    # 0) Only for daily frequency data that does not allow negatives (as
    #    precipitacion). Get a pp_threshold for the modeled series that result
    #    in the same amount of rain-day values in the historical period as in
    #    the observed series
    if (frq=='D') & (not allow_negatives):
        pp_threshold_mod = get_pp_threshold_mod(obs, mod, pp_threshold)
    else:
        pp_threshold_mod = pp_threshold
        
    # 1) Format series to matrix with rows as sub-periods and columns as years
    #    and, if needed, replace no-rain values with random small values 
    y_obs, obs_series = formatQM(obs, allow_negatives, frq, pp_threshold, pp_factor)
    y_mod, mod_series = formatQM(mod, allow_negatives, frq, pp_threshold_mod, pp_factor)
    
    # 2) Get statistics of the observed and modeled series in the historical
    #    period.
    if frq=='D':
        # Get a 3D array with a moving window centered on each day
        obs_series_moving = day_centered_moving_window(obs_series, day_win)
        mod_series_moving = day_centered_moving_window(mod_series, day_win)
        
        # Get a 4D array with a backward moving window for each projected period
        win_series = projected_backward_moving_window(mod_series_moving, y_obs, frq)
        
        # Reshape to 2D in order to have all the days of the window in each row
        obs_series_moving = obs_series_moving.reshape(365,(2*day_win-1)*y_obs)
        modh_series_moving = mod_series_moving[:,:,:y_obs].reshape(365,(2*day_win-1)*y_obs)
        
        # Reshape to 3D in order to have all the days of the window in each row
        win_series_moving = win_series.reshape(365,(2*day_win-1)*y_obs,y_mod-y_obs)
        
        # Replace no-rain values with nans but keep at least 30 values to get
        # the statistics. If fewer than 30 values are rain values, the missing
        # ones will be previously random values replaced
        if not allow_negatives:
            obs_series_moving = set_norain_to_nan(obs_series_moving, pp_threshold, pp_factor)
            modh_series_moving = set_norain_to_nan(modh_series_moving, pp_threshold_mod, pp_factor)
            mod_series[mod_series<pp_threshold_mod] = np.nan
            win_series_moving = set_norain_to_nan(win_series_moving, pp_threshold_mod, pp_factor)
        
        # Get the statistics for each row
        mu_obs, std_obs, skew_obs, skewy_obs = getStats(obs_series_moving)
        mu_mod, std_mod, skew_mod, skewy_mod = getStats(modh_series_moving)
        mu_win, std_win, skew_win, skewy_win = getStats(win_series_moving)
    else:
        # Get a 3D array with a backward moving window for each projected period
        win_series = projected_backward_moving_window(mod_series, y_obs, frq)
        
        # Get the statistics for each row
        mu_obs, std_obs, skew_obs, skewy_obs = getStats(obs_series)
        mu_mod, std_mod, skew_mod, skewy_mod = getStats(mod_series[:,:y_obs])
        mu_win, std_win, skew_win, skewy_win = getStats(win_series)
    
    # 3) Assign a probability distribution function to each sub-period for the
    #    observed and modeled data in the historical period.
    if user_pdf==False:
        if frq=='D':
            pdf_obs, ks_fail_obs = getDist(obs_series_moving, allow_negatives, mu_obs, std_obs, skew_obs, skewy_obs)
            pdf_mod, ks_fail_mod = getDist(modh_series_moving, allow_negatives, mu_mod, std_mod, skew_mod, skewy_mod)
        else:
            pdf_obs, ks_fail_obs = getDist(obs_series, allow_negatives, mu_obs, std_obs, skew_obs, skewy_obs)
            pdf_mod, ks_fail_mod = getDist(mod_series[:,:y_obs], allow_negatives, mu_mod, std_mod, skew_mod, skewy_mod)
    else:
        pdf_obs = np.zeros(obs_series.shape[0]) + pdf_obs
        pdf_mod = np.zeros(mod_series.shape[0]) + pdf_mod
    
    # 4) For each projected period assign a probability distribution function to
    #    to each sub-period and apply the cumulative distribution function to
    #    the modeled data.
    ks_fail_win = 0
    pdf_win = np.zeros((mod_series.shape[0], mod_series.shape[1]-y_obs))
    prob = np.zeros((mod_series.shape[0], y_mod-y_obs))
    for j in range(prob.shape[1]):
        # Assign a probability distribution function to each month.
        if user_pdf==False:
            if frq=='D':
                pdf_win[:,j], ks_fail_wtemp = getDist(win_series_moving[:,:,j], allow_negatives, mu_win[:,j], std_win[:,j], skew_win[:,j], skewy_win[:,j])
            else:
                pdf_win[:,j], ks_fail_wtemp = getDist(win_series[:,:,j], allow_negatives, mu_win[:,j], std_win[:,j], skew_win[:,j], skewy_win[:,j])
            ks_fail_win = ks_fail_win + ks_fail_wtemp
        else:
            pdf_win[:,j] = pdf_mod
        
        # Apply the cumulative distribution function of the projected period,
        # evaluated with the statistics of this period, to the last data of the
        # period. Equation 3 in Cannon et al. (2015).
        if frq=='D':
            prob[:,j] = getCDF(pdf_win[:,j],mod_series[:,y_obs+j:y_obs+j+1], mu_win[:,j], std_win[:,j], skew_win[:,j], skewy_win[:,j])[:,0]
        else:
            prob[:,j] = getCDF(pdf_win[:,j], win_series[:,-1:,j], mu_win[:,j], std_win[:,j], skew_win[:,j], skewy_win[:,j])[:,0]

    # 5) Apply the inverse cumulative distribution function:
    #    a) Of the observed data, evaluated with the statistics of the observed
    #       data in the historical period, to the probabilities obtained from 
    #       4b). Equation 5 of Cannon et al. (2015).
    inv_obs = getCDFinv(pdf_obs, prob, mu_obs, std_obs, skew_obs, skewy_obs)
    
    #    b) Of the modeled data, evaluated with the statistics of the observed
    #       data in the historical period, to the probabilities obtained from 
    #       4b). Equations 4 of Cannon et al. (2015).
    inv_mod = getCDFinv(pdf_mod, prob, mu_mod, std_mod, skew_mod, skewy_mod)
    
    # 6) Get the delta factor or relative change and apply it to the value
    #    obtained in 4b). Equation 4 and 6 of Cannon et al. (2015).
    if mult_change:
        delta_quantile = mod_series[:,y_obs:]/inv_mod
        bool_undefined = (delta_quantile>rel_change_th) & (inv_mod<inv_mod_th)
        delta_quantile[bool_undefined] = rel_change_th
        QDM_series = inv_obs*delta_quantile
    else:
        delta_quantile = mod_series[:,y_obs:] - inv_mod
        QDM_series = inv_obs + delta_quantile
    
    # 7) Perform QM for the historical period.
    mod_h = mod[:obs.shape[0]]
    QM_series = QM(obs, mod_h, allow_negatives, frq,pp_threshold, pp_factor, day_win, user_pdf, pdf_obs, pdf_mod)
    
    # 8) Reshape to a column vector
    QDM_series = QDM_series.reshape(-1,order='F')
    
    # 9) Concat QM (hsitorical) and QDM (future)
    QDM_series = np.hstack([QM_series, QDM_series])
    
    # 10) If does not allow negative values, replace nans (should be the result
    #     of correcting removed no-rain values) and replace no-rain values with 0
    if not allow_negatives:
        QDM_series[np.isnan(QDM_series)] = 0
        QDM_series[QDM_series<pp_threshold] = 0
    
    if user_pdf==False:
        # Check if ks-test failures
        ks_fail = ks_fail_obs + ks_fail_mod + ks_fail_win
        if ks_fail>0:
            print('QDM: Some of the probability distribution functions did not pass the KS-Test')
    
    return QDM_series


def UQM(obs, mod, mult_change=1, allow_negatives=1, frq='A', pp_threshold=1,
        pp_factor=1/100, day_win=1, user_pdf=False, pdf_obs=None, pdf_mod=None):
    """
    This function performs bias correction of modeled series based on observed
    data by the Unbiased Quantile Mapping (UQM) method, as described by Chadwick
    et al. (2023). Correction is performed to daily, monthly or annual data in a
    single location. An independent probability distribution function is 
    assigned to each sub-annual period based on the Kolmogorov-Smirnov test.
    Each projected period considers a window whose length is equal to the number
    of years of the historical period and ends in the analyzed year.
    
    Correction of the historical period is performed by the Quantile Mapping 
    (QM) method, as described in Cannon et al. (2015) and Chadwick et al. (2023).

    Inputs:
       obs:    A column vector of daily, monthly or annual observed data,
               without considering leap days. If daily frequency is specified,
               the length of the column vector should by a multiple of 365 and
               for monthly frequency, it should be a multiple of 12.
               [ndata_obs, 1]

       mod:    A column vector of daily, monthly or annual modeled or GCM data,
               without considering leap days. If daily frequency is specified,
               the length of the column vector should by a multiple of 365 and
               for monthly frequency, it should be a multiple of 12.
               [ndata_mod, 1]

       NOTE: This routine considers that obs and mod series start in the same
       day/month/year and are continuous until the end day/month/year.

    Optional inputs:
        mult_change:     A flag that indicates if projected changes should be
                         computed as multiplicative (fut = hist*delta) or 
                         additive (fut = hist + delta) changes.
                             mult_change = 1 or True: Multiplicative (default)
                             mult_change = 0 or False: Additive

        allow_negatives: A flag that identifies if data allows negative values
                         and also to replace no-rain values with random small 
                         values (Chadwick et al., 2023) to avoid numerical
                         problems with the probability distribution functions.
                             allow_negatives = 1 or True: Allow negatives (default)
                             allow_negatives = 0 or False: Do not allow negative

        frq:             A string specifying if the input frequency is daily,
                         monthly or annual.
                             Daily:     frq = 'D'
                             Monthly:   frq = 'M'
                             Annual:    frq = 'A' (default)

        pp_threshold:    A float indicating the threshold to consider no-rain
                         values. Default is 1.
                
        pp_factor:       A float which multiplied to pp_threshold indicates the
                         maximum value of the random values that replace no-rain
                         values. Default is 1/100.

        day_win:         (only for frq='D') An integer indicating how many days
                         to consider backwards and forward to get the statistics
                         of each calendar day.The length of the window will be 
                         (2*win_day-1). For example, day_win=15 -> window of 29.
                         Default: win = 1

        user_pdf:        A flag indicating if the user will define the
                         probability distribution functions (pdf) for the
                         observed and modeled series. The distributions will be
                         the same for all periods and sub-periods.
                            user_pdf = 1 or True: User defines the pdf
                            user_pdf = 0 or False: pdf defined by the Kolmogorov
                                                   -Smirnov test (default)

        pdf_obs:         An integer indicating the probability distribution 
                         function (pdf) to be used for the observed data. The
                         pdf will be the same for all periods and sub-periods.
                         Default: None

        pdf_mod:         An integer indicating the probability distribution 
                         function (pdf) to be used for the modeled data. The
                         pdf will be the same for all periods and sub-periods.
                         Default: None

        NOTE: The available distributions for pdf_obs and pdf_mod are:
            0) Normal
            1) Log-Normal
            2) Gamma 2 parameters
            3) Gamma 3 parameters
               (Pearson 3 parameters)
            4) Log-Gamma 3 parameters
               (Log-Pearson 3 parameters)
            5) Gumbel
            6) Exponential

    Output:
        UQM_series: A column vector of data bias corrected with the UQM method.
                   [ndata_mod, 1]
    
    References:
        Chadwick, C., Gironás, J., González-Leiva, F., and Aedo, S. (2023). Bias
        adjustment to preserve changes in variability: the unbiased mapping of
        GCM changes. Hydrological Sciences Journal,
        https://doi.org/10.1080/02626667.2023.2201450
    """
    
    # 0) Only for daily frequency data that does not allow negatives (as
    #    precipitacion). Get a pp_threshold for the modeled series that result
    #    in the same amount of rain-day values in the historical period as in
    #    the observed series
    if (frq=='D') & (not allow_negatives):
        pp_threshold_mod = get_pp_threshold_mod(obs, mod, pp_threshold)
    else:
        pp_threshold_mod = pp_threshold
    
    # 1) Format series to matrix with rows as sub-periods and columns as years
    #    and, if needed, replace no-rain values with random small values 
    y_obs, obs_series = formatQM(obs, allow_negatives, frq, pp_threshold, pp_factor)
    y_mod, mod_series = formatQM(mod, allow_negatives, frq, pp_threshold_mod, pp_factor)
    
    # 2) Get statistics of the observed and modeled series in the historical
    #    period.
    if frq=='D':
        # Get a 3D array with a moving window centered on each day
        obs_series_moving = day_centered_moving_window(obs_series, day_win)
        mod_series_moving = day_centered_moving_window(mod_series, day_win)
        
        # Get a 4D array with a backward moving window for each projected period
        win_series = projected_backward_moving_window(mod_series_moving, y_obs, frq)
        
        # Reshape to 2D in order to have all the days of the window in each row
        obs_series_moving = obs_series_moving.reshape(365,(2*day_win-1)*y_obs)
        modh_series_moving = mod_series_moving[:,:,:y_obs].reshape(365,(2*day_win-1)*y_obs)
        
        # Reshape to 3D in order to have all the days of the window in each row
        win_series_moving = win_series.reshape(365,(2*day_win-1)*y_obs,y_mod-y_obs)
        
        # Replace no-rain values with nans but keep at least 30 values to get
        # the statistics. If fewer than 30 values are rain values, the missing
        # ones will be previously random values replaced
        if not allow_negatives:
            obs_series_moving = set_norain_to_nan(obs_series_moving, pp_threshold, pp_factor)
            modh_series_moving = set_norain_to_nan(modh_series_moving, pp_threshold_mod, pp_factor)
            mod_series[mod_series<pp_threshold_mod] = np.nan
            win_series_moving = set_norain_to_nan(win_series_moving, pp_threshold_mod, pp_factor)

        # Get the statistics for each row
        mu_obs, std_obs, skew_obs, skewy_obs = getStats(obs_series_moving)
        mu_mod, std_mod, skew_mod, skewy_mod = getStats(modh_series_moving)
        mu_win, std_win, skew_win, skewy_win = getStats(win_series_moving)
    else:
        # Get a 3D array with a backward moving window for each projected period
        win_series = projected_backward_moving_window(mod_series, y_obs, frq)
        
        # Get the statistics for each row
        mu_obs, std_obs, skew_obs, skewy_obs = getStats(obs_series)
        mu_mod, std_mod, skew_mod, skewy_mod = getStats(mod_series[:,:y_obs])
        mu_win, std_win, skew_win, skewy_win = getStats(win_series)
    
    # 3) Assign a probability distribution function to each sub-period for the
    #    observed and modeled data in the historical period.
    if user_pdf==False:
        if frq=='D':
            pdf_obs, ks_fail_obs = getDist(obs_series_moving, allow_negatives, mu_obs, std_obs, skew_obs, skewy_obs)
        else:
            pdf_obs, ks_fail_obs = getDist(obs_series, allow_negatives, mu_obs, std_obs, skew_obs, skewy_obs)
    else:
        pdf_obs = np.zeros(obs_series.shape[0]) + pdf_obs
        pdf_mod = np.zeros(mod_series.shape[0]) + pdf_mod
    
    # 4) For each projected period, get the delta factor (delta) and time
    #    dependent (projected) statistics (mean, standard deviation, skewness, 
    #    and log skewness). Equations 13 to 14 in Chadwick et al. (2023).
    if mult_change:  # Precipitation
        delta_mu = mu_win/np.tile(mu_mod,(y_mod-y_obs,1)).T
        delta_sigma = std_win/np.tile(std_mod,(y_mod-y_obs,1)).T
        delta_skew = skew_win/np.tile(skew_mod,(y_mod-y_obs,1)).T
        delta_skewy = skewy_win/np.tile(skewy_mod,(y_mod-y_obs,1)).T
        
        mu_projected = np.tile(mu_obs,(y_mod-y_obs,1)).T*delta_mu
        sigma_projected = np.tile(std_obs,(y_mod-y_obs,1)).T*delta_sigma
        skew_projected = np.tile(skew_obs,(y_mod-y_obs,1)).T*delta_skew
        skewy_projected = np.tile(skewy_obs,(y_mod-y_obs,1)).T*delta_skewy
        
    else:
        delta_mu = mu_win - np.tile(mu_mod,(y_mod-y_obs,1)).T
        delta_sigma = std_win - np.tile(std_mod,(y_mod-y_obs,1)).T
        delta_skew = skew_win - np.tile(skew_mod,(y_mod-y_obs,1)).T
        delta_skewy = skewy_win- np.tile(skewy_mod,(y_mod-y_obs,1)).T
    
        mu_projected = np.tile(mu_obs,(y_mod-y_obs,1)).T + delta_mu
        sigma_projected = np.tile(std_obs,(y_mod-y_obs,1)).T + delta_sigma
        skew_projected = np.tile(skew_obs,(y_mod-y_obs,1)).T + delta_skew
        skewy_projected = np.tile(skewy_obs,(y_mod-y_obs,1)).T + delta_skewy

    # 5) For each projected period assign a probability distribution function to
    #    to each sub-period and apply the cumulative distribution function to
    #    the modeled data.
    ks_fail_win = 0
    pdf_win = np.zeros((mod_series.shape[0], y_mod-y_obs))
    prob = np.zeros((mod_series.shape[0], y_mod-y_obs))
    UQM_series = np.zeros((mod_series.shape[0], y_mod-y_obs))
    for j in range(prob.shape[1]):
        
        # Dwin_series = detrend(win_series[:,:,j]) + np.tile(mu_win[:,j],(y_obs,1)).T
        # mu_win_, std_win_, skew_win_, skewy_win_ = getStats(Dwin_series)
        
        # Assign a probability distribution function to each month.
        if user_pdf==False:
            if frq=='D':
                pdf_win[:,j], ks_fail_wtemp = getDist(win_series_moving[:,:,j], allow_negatives, mu_win[:,j], std_win[:,j], skew_win[:,j], skewy_win[:,j])
            else:
                pdf_win[:,j], ks_fail_wtemp = getDist(win_series[:,:,j], allow_negatives, mu_win[:,j], std_win[:,j], skew_win[:,j], skewy_win[:,j])
            ks_fail_win = ks_fail_win + ks_fail_wtemp
        else:
            pdf_win[:,j] = pdf_mod
        
        # Apply the cumulative distribution function of the projected period,
        # evaluated with the statistics of this period, to the last data of the
        # period. Equation 3 in Chadwick et al. (2023).
        if frq=='D':
            prob[:,j] = getCDF(pdf_win[:,j],mod_series[:,y_obs+j:y_obs+j+1], mu_win[:,j], std_win[:,j], skew_win[:,j], skewy_win[:,j])[:,0]
        else:
            prob[:,j] = getCDF(pdf_win[:,j], win_series[:,-1:,j], mu_win[:,j], std_win[:,j], skew_win[:,j], skewy_win[:,j])[:,0]
        # Taot[:,j] = getCDF(pdf_win[:,j], Dwin_series[:,-1:], mu_win_, std_win_, skew_win_, skewy_win_)[:,0]
        
        # 6) Apply the inverse cumulative distribution function of the observed
        #    data, evaluated with the time dependent statistics, to the values
        #    obtained in 3) (getCDFinv function of the climQMBC package). Equation
        #    15 in Chadwick et al. (2023).
        UQM_series[:,j] = getCDFinv(pdf_obs, prob[:,j:j+1], mu_projected[:,j], sigma_projected[:,j], skew_projected[:,j], skewy_projected[:,j])[:,0]
    
    # 7) Perform QM for the historical period.
    mod_h = mod[:obs.shape[0]]
    QM_series = QM(obs, mod_h, allow_negatives, frq,pp_threshold, pp_factor, day_win, user_pdf, pdf_obs, pdf_mod)
    
    # 8) Reshape to a column vector
    UQM_series = UQM_series.reshape(-1, order='F')
    
    # 9) Concat QM (hsitorical) and UQM (future)
    UQM_series = np.hstack([QM_series, UQM_series])
    
    # 10) If does not allow negative values, replace nans (should be the result
    #     of correcting removed no-rain values) and replace no-rain values with 0
    if not allow_negatives:
        UQM_series[np.isnan(UQM_series)] = 0
        UQM_series[UQM_series<pp_threshold] = 0
        
    if user_pdf==False:
        # Check if ks-test failures
        ks_fail = ks_fail_obs + ks_fail_win
        if ks_fail>0:
            print('UQM: Some of the probability distribution functions did not pass the KS-Test')
    
    return UQM_series


def SDM(obs, mod, SDM_var, frq='A', pp_threshold=1, pp_factor=1/100, day_win=1):
    """
    This function performs bias correction of modeled series based on observed
    data by the Scaled Distribution Mapping (SDM) method, as described by 
    Switanek et al. (2017). Correction is performed to daily, monthly or annual 
    precipitation or temperature data in a single location. The historical 
    period of the modeled series is bias corrected as a scenario where the 
    complete series is replaced by the bias corrected series. On the other 
    hand, for each year after the historical period the method considers a 
    projected period as a window whose length is equal to the number of years 
    of the historical period and ends in the analyzed year. From the bias 
    corrected series, only the last year is saved and all the others are 
    discarded. This approach aims to build a continuum for the bias corrected     
    series.
    
    Inputs:
        obs:     A column vector of daily, monthly or annual observed data,
                 without considering leap days. If daily frequency is specified,
                 the length of the column vector should by a multiple of 365 and
                 for monthly frequency, it should be a multiple of 12.
                 [ndata_obs, 1]

        mod:     A column vector of daily, monthly or annual modeled or GCM data,
                 without considering leap days. If daily frequency is specified,
                 the length of the column vector should by a multiple of 365 and
                 for monthly frequency, it should be a multiple of 12.
                 [ndata_mod, 1]

        SDM_var: A flag that identifies if data are temperature or precipitation.
                     Temperature:   var = 0
                     Precipitation: var = 1

        NOTE: This routine considers that obs and mod series start in the same
        day/month/year and are continuous until the end day/month/year.

    Optional inputs:
        frq:             A string specifying if the input frequency is daily,
                         monthly or annual.
                             Daily:     frq = 'D'
                             Monthly:   frq = 'M'
                             Annual:    frq = 'A' (default)

        pp_threshold:    A float indicating the threshold to consider no-rain
                         values. Default is 1.
                
        pp_factor:       A float which multiplied to pp_threshold indicates the
                         maximum value of the random values that replace no-rain
                         values. Default is 1/100.
                         
        day_win:         (only for frq='D') An integer indicating how many days
                         to consider backwards and forward to get the statistics
                         of each calendar day.The length of the window will be 
                         (2*win_day-1). For example, day_win=15 -> window of 29.
                         Default: win = 1

    Output:
        SDM_series: A column vector of data bias corrected with the QM method.
                    [ndata_mod, 1]
    
    References:
        Switanek, B. M., Troch, P., Castro, L. C., Leuprecht, A., Chang, H. I.,
        Mukherjee, R., and Demaria, M. C. E. (2017) Scaled distribution
        mapping: A bias correction method that preserves raw climate model
        projected changes. Hydrology &amp; Earth System Sciences, 21,
        2649-2666, https://doi.org/10.5194/hess-21-2649-2017.
    """
    # Define the no-rain value and the threshold to truncate the tail of the
    # probability ditribution functions
    lower_lim = pp_threshold
    CDF_th = 1e-3
    
    # 0) Format series to matrix with rows as sub-periods and columns as years
    #    and, if needed, replace no-rain values with random small values 
    y_obs, obs_series = formatQM(obs, SDM_var, frq, pp_threshold, pp_factor)
    y_mod, mod_series = formatQM(mod, SDM_var, frq, pp_threshold, pp_factor)
    modh_series= mod_series[:,:y_obs]
    
    # If frequency is daily, consider a moving window for each day to fit
    # the statistics
    if frq=='D':
        # Get a 3D array with a moving window centered on each day
        obs_series = day_centered_moving_window(obs_series, win=day_win)
        mod_series = day_centered_moving_window(mod_series, win=day_win)
        
        # Get a 4D array with a backward moving window for each projected period
        win_series = projected_backward_moving_window(mod_series, y_obs, frq)
        
        # Reshape to 2D in order to have all the days of the window in each row
        obs_series = obs_series.reshape(365,(2*day_win-1)*y_obs)
        modh_series= mod_series[:,:,:y_obs].reshape(365,(2*day_win-1)*y_obs)
        
        # Reshape to 3D in order to have all the days of the window in each row
        win_series= win_series.reshape(365,(2*day_win-1)*y_obs,y_mod-y_obs)
        
        # Add the historical period to the moving window
        win_series_= np.dstack([modh_series, win_series])
        
    else:
        # For non-daily frequency, make sure that day_win=1
        day_win=1

        # Get a 3D array with a backward moving window for each projected period
        win_series = projected_backward_moving_window(mod_series, y_obs, frq)
        
        # Add the historical period to the moving window
        win_series_= np.dstack([modh_series, win_series])
    
    SDM_series = np.zeros((mod_series.shape[0], y_mod-y_obs))    
    SDM_h_series = np.zeros((obs_series.shape[0], y_obs))
    for m in range(mod_series.shape[0]):
        # 1) Historical period:

        # a) [Switanek et al. (2017), step 1)]
        #     For temperature, get the detrended modeled and observed series in
        #     the historical period.
        #     For precipitation, get the rainday values and its frequency for
        #     the modeled and observed series in the historical period.
         
        if SDM_var==0: # Temperature
            D_obs = detrend(obs_series[m])
            D_mod = detrend(modh_series[m])
            
            mu_obs = np.nanmean(obs_series[m])
            mu_mod = np.nanmean(modh_series[m])
        else:       # Precipitation
            D_obs = np.sort(obs_series[m][obs_series[m]>lower_lim])
            D_mod = np.sort(modh_series[m][modh_series[m]>lower_lim])
            
            freq_obs = D_obs.shape[0]/obs_series.shape[1]
            freq_mod = D_mod.shape[0]/modh_series[m].shape[0]
            
            if freq_mod==0:
                freq_mod = 1/(365*y_obs)
    
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
        if SDM_var==0: # Temperature
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
        for j in range(-1, y_mod-y_obs):
            # c) Initialize correction array.
            corr_temp = np.zeros(y_obs*(2*day_win-1))
        
            # d) Define projected window.
            win_series = win_series_[m,:,j+1]
            
            # e) [Switanek et al. (2017), step 1)]
            #    For temperature, get the detrended series of the projected
            #    period.
            #    For precipitation, get the rainday values, its frequency, and
            #    expected raindays for the projected period.
            #    Get the index of the sorted detrended or rainday values.
            if SDM_var==0: # Temperature
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
            if SDM_var==0: # Temperature
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
            if SDM_var==0: # Temperature
                SF = (sigma_obsD/sigma_modD)*(stat.norm.ppf(CDF_win, *(mu_winD,sigma_winD)) - stat.norm.ppf(CDF_win, *(mu_modD,sigma_modD)))
            else: # Precipitation
                SF = stat.gamma.ppf(CDF_win, *fit_win)/stat.gamma.ppf(CDF_win, *fit_mod)
    
            # h) Interpolate observed and modeled CDF of the historical period
            #    to the length of the projected period.
            if len(D_obs)>0:
                obs_cdf_intpol = np.interp(np.linspace(1, len(D_obs), len(D_win)),
                                           np.linspace(1, len(D_obs), len(D_obs)),
                                           CDF_obs)
            else:
                obs_cdf_intpol = CDF_win*0
                
            if len(D_mod)>0:
                mod_cdf_intpol = np.interp(np.linspace(1, len(D_mod), len(D_win)),
                                           np.linspace(1, len(D_mod), len(D_mod)),
                                           CDF_mod)
            else:
                mod_cdf_intpol = CDF_win*0

            # i) [Switanek et al. (2017), steps 4 and 5)]
            #    Get recurrence intervals and its scaled version for the
            #    observed, historical modeled and projected period modeled
            #    CDFs.
            if SDM_var==0: # Temperature
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
            if SDM_var==0: # Temperature
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
            if exp_D>0: # if exp_D==0 -> corr_temp stays a a zero vector
                corr_temp[win_argsort[-exp_D:]] = xvals
                
            if SDM_var==0:
                corr_temp = corr_temp + diff_win - mu_win
                
            # l) If the projected period is the historical period (period -1,
            #    j=-1) save the complete bias corrected series.
            #    If the projected period is not the historical period (j>-1),
            #    save the value of the last year.
            if j==-1:
                SDM_h_series[m] = corr_temp[(day_win-1)::(2*day_win-1)]
            else:
                SDM_series[m,j] = corr_temp[-1]
    
    # Reshape to a column vector
    SDM_h_series = SDM_h_series.reshape(-1, order='F')
    SDM_series = SDM_series.reshape(-1, order='F')
    
    # Concat historical and future
    SDM_series = np.hstack([SDM_h_series, SDM_series])
    
    # If precipitation, replace no-rain values with 0
    if SDM_var==1:
        SDM_series[SDM_series<pp_threshold] = 0

    return SDM_series