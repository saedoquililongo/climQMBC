# -*- coding: utf-8 -*-
"""
This script contains the main functions requiered by the methods implemented in
the climQMBC package, including the formatQM, the getStats, the getDist,
the getCDF, the getCDFinv and projected_backward_moving_window functions. For
daily precipitation data, it also contains the get_pp_threshold_mod and
 set_norain_to_nan functions.


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

import scipy.stats as stat
import numpy as np


def get_pp_threshold_mod(obs, mod, pp_threshold):
    """
    This functions gets a no-rain value threshold for the modeled series by 
    mathcing the no-rain days of the modeled series in the historical period
    with the observed series.

    Inputs:
       obs:             A column vector of daily observed data, without 
                        considering leap  days. The length of the column vector 
                        should by a multiple of 365. [ndata_obs, 1]

       mod:             A column vector of daily modeled data, without 
                        considering leap  days. The length of the column vector 
                        should by a multiple of 365. [ndata_mod, 1]

        pp_threshold:   A float indicating the threshold to consider no-rain
                        values in the observed data.
    
    Output:
        pp_threshold_mod:   A float indicating the threshold to consider no-rain
                            values in the modeled data.
    """
    
    # Get the rain days of the observed series
    obs_rainday_hist = np.sum(obs>pp_threshold)
    
    # Define the pp_threshold of the modeled series by matching the rain days
    # of the observed series
    mod_sort_descending = np.sort(mod[:obs.shape[0]],0)[::-1]
    days_kept = min(obs_rainday_hist,len(obs)-1)
    pp_threshold_mod = mod_sort_descending[days_kept,0]
    
    return pp_threshold_mod


def formatQM(series_, allow_negatives, frq, pp_threshold, pp_factor):
    """
    This function formats a time series from a column vector to a matrix. The
    number of rows depends on the frequency of the data (1, 12 and 365 for 
    annual, monthyl and daily data, respectively). If negative values are not
    allowed, determined by allow_negatives=1, no-rain values (values below
    pp_threshold) are replaced by random positive values below the product
    pp_threshold*pp_factor.

    Inputs:
        series_:         A column vector of daily data, without considering leap
                         days. The length of the column vector should by a
                         multiple of 365. [ndata_obs, 1]

        allow_negatives: A flag that identifies if data allows negative values
                         and also to replace no-rain values with random small 
                         values (Chadwick et al., 2023) to avoid numerical
                         problems with the probability distribution functions.
                             allow_negatives = 1 or True: Allow negatives
                             allow_negatives = 0 or False: Do not allow negative

        frq:             A string specifying if the input frequency is daily,
                         monthly or annual.
                             Daily:     frq = 'D'
                             Monthly:   frq = 'M'
                             Annual:    frq = 'A'

        pp_threshold:    A float indicating the threshold to consider no-rain
                         values.
                
        pp_factor:       A float which multiplied to pp_threshold indicates the
                         maximum value of the random values that replace no-rain
                         values.

    Output:
        years:          Number of years of the series.

        series_matrix:  A matrix of daily, monthly or annual observed data,
                        without considering leap days. For daily, monthly and
                        annual frequency, the number of rows will be 365, 12 and
                        1, respectively. The number of columns will be equal to
                        the number of years.
                        [1, 12 or 365, years]
    """

    # This prevents modyfing the original input passed by reference
    ## Must look for a more elegant way to skip this step
    series_matrix = series_.copy()

    # Check frequency and define the number of rows
    if frq=='D':
        I = 365
    elif frq=='M':
        I = 12
    elif frq=='A':
        I = 1
    else:
        I = 1

    # If variable is precipitation, determined by allow_negatives=1, replace low
    # values with random values.
    if not allow_negatives:
        bool_low = series_matrix<pp_threshold
        series_matrix[bool_low] = np.random.rand(bool_low.sum())*pp_factor*pp_threshold

    # Get number of years of the series.
    years = int(series_matrix.shape[0]/I)
    
    # Reshape to get a matrix of shape [I , years]
    series_matrix = series_matrix.reshape((I,years), order='F')

    return years, series_matrix


def getStats(series):
    """
    This function computes the mean, standard deviation, skewness and skewness
    of the logarithmic values, for each sub-period within the year, according
    to the frequency initially defined. Statistics are computed along axis 1. 

    Inputs:
        series:         An array of daily, monthly or annual data, without 
                        considering leap days. Possible input dimensions:
                        2D array:
                            [sub-periods, years]
                            [sub-periods + days window, years]
                        3D array:
                            [sub_periods, projected periods window, years]
                            [sub-periods + days window, projected periods window, years]

    Output:
        mu:         Mean values of each sub-period (and projected period if
                    input is a 3D array).

        sigma:      Standar deviation values of each sub-period (and projected 
                    period if input is a 3D array).

        skew:       Skewness values of each sub-period (and projected period if
                    input is a 3D array).

        skewy:      Skewness values of the logarithm of theseries of each 
                    sub-period (and projected period if input is a 3D array).
    
        NOTE: If input series is a 2D array, outputs will be a column vector. If
        input series is a 3D array, outputs will be a 2D array.
    """

    # Get the mean, standard deviation, skewness and skewness of the
    # logarithmic values of each year sub-period of the series.
    mu  = np.nanmean(series, 1)                     # Mean
    sigma = np.nanstd(series, 1, ddof=1)            # Standard deviation
    skew = stat.skew(series, 1, bias=False, nan_policy='omit')         # Skewness
    series_log  = np.log(series)
    series_log[np.isinf(series_log)] = np.log(np.random.rand(np.isinf(series_log).sum())*0.01)
    skewy = stat.skew(series_log, 1, bias=False, nan_policy='omit')    # Skewness of the log

    return mu, sigma, skew, skewy


def day_centered_moving_window(series, win):
    """
    This function transform a 2D matrix of dimension [days in year, years] to a 
    3D array, where the new dimension is a centered moving window for each day
    of the year.

    Inputs:
        series_matrix:  A matrix of daily data, without considering leap days. 
                        The number of rows is 365 and number of columns is the 
                        number of years of the series.
                        [365, years]
                            
        day_win:        An integer indicating how many days to consider 
                        backwards and forward to get the statistics of each 
                        calendar day.The length of the window will be 
                        (2*win_day-1). For example, day_win=15 -> window of 29.

    Output:
        series_moving:  A 3D array of daily data, without considering leap days,
                        with a dimension that considers a centered moving window
                        for each day of the year.
                        [365, 2*win-1, years]
    """
        
    series_moving = np.vstack([series[-win:],np.tile(series,(win*2,1)),series[:win]])
    series_moving = series_moving.reshape(series.shape[0]+1,win*2,series.shape[1], order='F')[:-1,1:]
    
    return series_moving


def projected_backward_moving_window(series, projected_win, frq):
    """
    This function transform a 2D or 3D array of dimensions [days in year, year]
    or [days in year, centered window for each day, years] to a 3D or 4D array
    adding a new dimension to represent a backward moving window for each 
    projected period.

    Inputs:
        series:         A 2D or 3D array of daily, monthly or annual observed 
                        data, without considering leap days. For daily, monthly
                        and annual frequency, the number of rows will be 365, 12
                        and 1, respectively. If daily data, it is a 3D array
                        that has a  dimension with a centered moving window for
                        each day of the year.
                        [1, 12 or 365, years] or [365, day centered window, years]

        projected_win:  An integer the length of the backward moving window.
                            
        frq:             A string specifying if the input frequency is daily,
                         monthly or annual.
                             Daily:     frq = 'D'
                             Monthly:   frq = 'M'
                             Annual:    frq = 'A'

    Output:
        win_series:     An array of daily, monthly or annual future data,
                        without considering leap days, with a dimension
                        associated to a backward moving window for each 
                        projected period. Possible dimensions:
                        3D array:
                            [sub-periods, projected periods window, years]
                        4D array:
                            [sub-periods, days window, projected periods window, years]
    """
    
    y_mod = series.shape[-1]

    if frq=='D':
        day_win = int((series.shape[1]+1)/2)
        win_series = np.dstack([np.tile(series,(1,1,projected_win)),np.zeros((series.shape[0],2*day_win-1,projected_win))])
        win_series = win_series.reshape(series.shape[0],2*day_win-1,projected_win,y_mod+1)[:,:,:,1:-projected_win]
    else:
        win_series = np.hstack([np.tile(series,(1,projected_win)), np.zeros((series.shape[0],projected_win))])
        win_series = win_series.reshape(series.shape[0],projected_win,y_mod+1)[:,:,1:-projected_win]
        
    return win_series


def set_norain_to_nan(series_moving, pp_threshold, pp_factor, min_rainday=30):
    """
    This function replace no-rain values with nans, leaving a minimum amout
    of values (min_rainday) to fit a distribution. If fewer than min_rainday
    values are above pp_threshold, random small values will be added.

    Inputs:
        series_moving:  A 2D or 3D array of daily, monthly or annual observed 
                        data, without considering leap days. For daily, monthly
                        and annual frequency, the number of rows will be 365, 12
                        and 1, respectively. If daily data, it is a 3D array
                        that has a  dimension with a centered moving window for
                        each day of the year.
                        [1, 12 or 365, years] or [365, day centered window, years]

        pp_factor:       A float which multiplied to pp_threshold indicates the
                         maximum value of the random values that replace no-rain
                         values.

        pp_threshold:    A float indicating the threshold to consider no-rain
                         values.

    Optional inputs:
        min_rainday:    Minimum amount of values to keep for each sub-period and
                        projected period, to ensure a minimum amount to fit a 
                        distribution. Default is 30

    Output:
        series_moving:  The input, but with no-rain values replaced with nans.
    """

    rainday_count = np.sum(series_moving>pp_threshold,1)
    
    if len(rainday_count.shape)==1:
        for per in range(series_moving.shape[0]):
            bool_low = series_moving[per]<pp_threshold
            replace_values_nans = np.random.rand(bool_low.sum())*pp_factor*pp_threshold
            replace_values_nans[max(0,min_rainday-rainday_count[per]):] = np.nan
            series_moving[per,bool_low] = replace_values_nans
    else:
        for per1 in range(rainday_count.shape[0]):
            for per2 in range(rainday_count.shape[1]):
                bool_low = series_moving[per1,:,per2]<pp_threshold
                replace_values_nans = np.random.rand(bool_low.sum())*pp_factor*pp_threshold
                replace_values_nans[max(0,min_rainday-rainday_count[per1,per2]):] = np.nan
                series_moving[per1,bool_low,per2] = replace_values_nans
    
    return series_moving


def getDist(series, allow_negatives, mu, sigma, skew, skewy):
    """
    This function assigns an independent probability distribution function to
    each row of the input series by comparing the empirical probability
    distribution function with seven distributions based on the
    Kolmogorov-Smirnov (KS) test. The available distributions 
    are:
        1) Normal
        2) Log-Normal
        3) Gamma 2 parameters
        4) Gamma 3 parameters
           (Pearson 3 parameters)
        5) Log-Gamma 3 parameters
           (Log-Pearson 3 parameters)
        6) Gumbel
        7) Exponential
        
        For allow_negatives=0, only 2), 3) and 5) are considered (1, 4, 6, and 7
        are discarded). For series with negative values, only 1), 3), 4), 6),
        and 7) are considered (2, 3 and 5 are discarded).

    Input:
        series:         An array of daily, monthly or annual data, without 
                        considering leap days. Possible input dimensions:
                        2D array:
                            [sub-periods, years]
                            [sub-periods + days window, years]
                            
        allow_negatives: A flag that identifies if data allows negative values
                         and also to replace no-rain values with random small 
                         values (Chadwick et al., 2023) to avoid numerical
                         problems with the probability distribution functions.
                             allow_negatives = 1 or True: Allow negatives
                             allow_negatives = 0 or False: Do not allow negative
        
        mu:         A vector with mean values of each sub-period.

        sigma:      A vector with standard deviation values of each sub-period.

        skew:       A vector with skewness values of each sub-period.

        skewy:      A vector with skewness values of the logarithm of the series
                    of each sub-period.

    Output:
        pdf:    A vector with an ID for the resulting distribution from
                the KS test. The ID is related to  the numeration of the
                distribution listed in the description of this function. This ID
                is used in the getCDF and getCDFinv  functions of the climQMBC
                package.
    """
    
    # Get the number of rows of the input series.
    nrows = series.shape[0]
    
    # Initialize column vectors for the statistics needed for the available
    # probability distribution functions.
    pdf = np.zeros(nrows)
    
    # Initialize a counter of failures of the KS test
    ks_fail = 0
    
    # Perform the Kolmogorov-Smirnov test for sub-period or row.
    for sp in range(nrows):
        series_sub = series[sp,:]
        
        # Get the number of years with valid data (for allow_negatives=0, the 
        # series might have nan values to ignore them)
        series_sub = series_sub[~np.isnan(series_sub)]
        y_series = len(series_sub)

        # a) Get empirical distribution.
        sortdata = np.sort(series_sub) 
        probEmp = np.linspace(1/(y_series+1),
                              y_series/(y_series+1),
                              y_series)

        # b) Compare distributions.
        # i) Normal distribution.
        normal = stat.norm.cdf(sortdata, mu[sp], sigma[sp])
        KSnormal = max(abs(normal-probEmp))

        # ii) Log Normal distribution.
        if (series_sub<0).any():
            KSlognormal = 1
        else:
            sigmay = np.sqrt(np.log(1 + (sigma[sp]/mu[sp])**2))
            muy = np.log(mu[sp]) - (sigmay**2)/2
            lognormal = stat.lognorm.cdf(sortdata, scale=np.exp(muy), s=sigmay)
            KSlognormal = max(abs(lognormal-probEmp))
        
        # iii) Gamma 2 parameters distribution.
        if (series_sub<0).any():
            KSgammaII = 1
        else:
            A = (sigma[sp]**2)/mu[sp] 
            B = (mu[sp]/sigma[sp])**2
            GammaII = stat.gamma.cdf(sortdata, a=B, scale=A)
            KSgammaII = max(abs(GammaII-probEmp))
    
        # iv) Gamma 3 parameters distribution.
        #     (Pearson 3 parameters distribution)
        Bet = (2/skew[sp])**2
        Alp = sigma[sp]/np.sqrt(Bet)
        Gam = mu[sp] - (Alp*Bet)
        GammaIII = stat.gamma.cdf(sortdata-Gam, a=Bet, scale=Alp)
        KSgammaIII = max(abs(GammaIII-probEmp))
    
        # v) Log-Gamma 3 parameters distribution.
        #    (Log-Pearson 3 parameters distribution)
        if (series_sub<0).any():
            KSLpIII = 1
        else:
            sigmay = np.sqrt(np.log(1 + (sigma[sp]/mu[sp])**2))
            muy = np.log(mu[sp]) - (sigmay**2)/2
            Bety = (2/skewy[sp])**2
            Alpy = sigmay/np.sqrt(Bety)
            Gamy = muy-(Alpy*Bety)
            sortdata_log = np.log(sortdata)
            sortdata_log[np.isinf(sortdata_log)] = np.log(np.random.rand(np.isinf(sortdata_log).sum())*0.01)
            LpIII = stat.gamma.cdf(sortdata_log-Gamy, a=Bety, scale=Alpy)
            KSLpIII = max(abs(LpIII-probEmp))

        # vi) Gumbel distribution.
        Sn = np.pi/np.sqrt(6)
        yn = 0.5772
        a = Sn/sigma[sp]
        u = mu[sp] - (yn/a)
        gumbel = np.exp(-np.exp(-a*(sortdata - u)))
        KSgumbel = max(abs(gumbel-probEmp))
    
        # vii) Exponential distribution.
        gamexp = mu[sp] - sigma[sp]
        exponential = 1-np.exp(-1/sigma[sp]*(sortdata - gamexp))
        exponential[exponential<0] = 0
        KSexponential= max(abs(exponential-probEmp))
        
        # c) allow_negatives=0, set KS=1 to distributions that allow
        #    negative values (this will discard those distributions).
        if not allow_negatives:
            KSnormal = 1
            KSgammaIII = 1
            KSgumbel = 1
            KSexponential = 1

        # d) The distribution with lower KS value is considered for each month.
        KS_vals = [KSnormal, KSlognormal, KSgammaII, KSgammaIII, KSLpIII, KSgumbel, KSexponential]
        min_KS = min(KS_vals)
        bestPDF = np.argmin(KS_vals)
        
        ks_crit = 1.3581/np.sqrt(y_series)
        ks_fail = ks_fail + (np.sign(min_KS-ks_crit)+1)/2

        pdf[sp] = bestPDF

    return pdf, ks_fail


def getCDF(pdf, series, mu, sigma, skew, skewy):
    """
    This function evaluates each row of the series in the respective cumulative
    distribution function assigned by the Kolmogorov-Smirnov (KS) test in the 
    getDist function of the climQMBC package. The available distributions are:
        1) Normal
        2) Log-Normal
        3) Gamma 2 parameters
        4) Gamma 3 parameters
           (Pearson 3 parameters)
        5) Log-Gamma 3 parameters
           (Log-Pearson 3 parameters)
        6) Gumbel
        7) Exponential

    Note that each row is associated to an independent distribution function.

    Input:
        pdf:     A vector with an ID for the resulting distribution from
                 the KS test. The ID is related to  the numeration of the
                 distribution listed in the description of this function.
        
        series:  An array of daily, monthly or annual data, without 
                 considering leap days. Possible input dimensions:
                 2D array:
                     [sub-periods, years]
                     [sub-periods + days window, years]
        
        mu:      A vector with mean values of each sub-period.

        sigma:   A vector with standard deviation values of each sub-period.

        skew:    A vector with skewness values of each sub-period.

        skewy:   A vector with skewness values of the logarithm of the series
                 of each sub-period.

    Output:
        Taot:   A matrix with the non-exceedance probability for value
                row of the input series.
    """

    # Get the number of rows and columns of the series.
    rows = series.shape[0]
    cols = series.shape[1]
    
    # Compute the cumulative distribution function to the values of each 
    # row, based on the distribution assigned in the getDist function of the
    # climQMBC package.
    prob = np.zeros((rows,cols))
    for sp in range(rows):
        series_sub = series[sp,:]
        
        if pdf[sp]==0: # i) Normal distribution.
            prob[sp,:] = stat.norm.cdf(series_sub, mu[sp], sigma[sp])
            
        elif pdf[sp]==1: # ii) Log-Normal distribution.
            sigmay = np.sqrt(np.log(1 + (sigma[sp]/mu[sp])**2))
            muy = np.log(mu[sp])-(sigmay**2)/2
            prob[sp,:] = stat.lognorm.cdf(series_sub, scale=np.exp(muy), s=sigmay)
            
        elif pdf[sp]==2: # iii) Gamma 2 parameters distribution.
            A = (sigma[sp]**2)/mu[sp] 
            B = (mu[sp]/sigma[sp])**2
            prob[sp,:] = stat.gamma.cdf(series_sub, a=B, scale=A)
            
        elif pdf[sp]==3: # iv) Gamma 3 parameters distribution.
            Bet = (2/skew[sp])**2
            Alp = sigma[sp]/np.sqrt(Bet)
            Gam = mu[sp]-(Alp*Bet)
            prob[sp,:] = stat.gamma.cdf(series_sub-Gam, a=Bet, scale=Alp)
            
        elif pdf[sp]==4: # v) Log-Gamma 3 parameters distribution.
            sigmay = np.sqrt(np.log(1 + (sigma[sp]/mu[sp])**2))
            
            muy = np.log(mu[sp]) - (sigmay**2)/2
            Bety = (2/skewy[sp])**2
            Alpy = sigmay/np.sqrt(Bety)
            Gamy = muy - (Alpy*Bety)
            Lnsortdata = np.log(series_sub)
            Lnsortdata[np.isinf(Lnsortdata)] = np.log(np.random.rand(np.isinf(Lnsortdata).sum())*0.01)
            prob[sp,:]  = stat.gamma.cdf(Lnsortdata-Gamy,a=Bety,scale=Alpy)
        elif pdf[sp]==5: # vi) Gumbel distribution.
            Sn = np.pi/np.sqrt(6)
            yn = 0.5772
            a = Sn/sigma[sp]
            u = mu[sp]-(yn/a)
            prob[sp,:] = np.exp(-np.exp(-a*(series_sub - u)))
    
        elif pdf[sp]==6: # vii) Exponential distribution.
            gamexp = mu[sp]-sigma[sp]
            exponential = 1-np.exp(-1/sigma[sp]*(series_sub-gamexp))
            exponential[exponential<0] = 0
            prob[sp,:] = exponential
    
    # Tail probabilities are set to (0 +  threshold) and (1 - threshold) to
    # avoid numerical errors.
    th = 1e-3
    prob[prob>1-th] = 1-th
    prob[prob<th] = th
    
    return prob


def getCDFinv(PDF, prob, mu, sigma, skew, skewy):
    """
    This function evaluates the probability, prob, in the respective inverse
    cumulative distribution function assigned by the Kolmogorov-Smirnov (KS)
    test in the getDist function of the climQMBC package. The available
    distributions are:
        1) Normal
        2) Log-Normal
        3) Gamma 2 parameters
        4) Gamma 3 parameters
           (Pearson 3 parameters)
        5) Log-Gamma 3 parameters
           (Log-Pearson 3 parameters)
        6) Gumbel
        7) Exponential
        
    Note that each row is associated to an independent distribution function.
        
    Input:
        pdf:     A vector with an ID for the resulting distribution from
                 the KS test. The ID is related to  the numeration of the
                 distribution listed in the description of this function.
        
        prob:    A matrix with the non-exceedance probability for value
                 row of the input series.
        
        mu:      A vector with mean values of each sub-period.

        sigma:   A vector with standard deviation values of each sub-period.

        skew:    A vector with skewness values of each sub-period.

        skewy:   A vector with skewness values of the logarithm of the series
                 of each sub-period.
            
    Output:
        xhat:   A matrix with the values obtained when the inverse
                cumulative distribution function is applied to prob.
    """

    # Get the number of rows and columns of the series.
    rows = prob.shape[0]
    cols = prob.shape[1]
    
    # Compute the cumulative distribution function to the values of each 
    # row, based on the distribution assigned in the getDist function of the
    # climQMBC package.
    xhat = np.zeros((rows,cols))
    for sp in range(rows):
        
        if PDF[sp]==0: # i) Normal distribution.
            xhat[sp,:] = stat.norm.ppf(prob[sp,:],mu[sp],sigma[sp])
            
        elif PDF[sp]==1: # ii) Log-Normal distribution.
            sigmay = np.sqrt(np.log(1 + (sigma[sp]/mu[sp])**2))
            muy = np.log(mu[sp])-(sigmay**2)/2
            xhat[sp,:] = stat.lognorm.ppf(prob[sp,:],scale=np.exp(muy),s=sigmay)
            
        elif PDF[sp]==2: # iii) Gamma 2 parameters distribution.
            A = (sigma[sp]**2)/mu[sp] 
            B = (mu[sp]/sigma[sp])**2    
            xhat[sp,:] = stat.gamma.ppf(prob[sp,:],a=B,scale=A)
            
        elif PDF[sp]==3: # iv) Gamma 3 parameters distribution.
            Bet = (2/skew[sp])**2
            Alp = sigma[sp]/np.sqrt(Bet)
            Gam = mu[sp]-(Alp*Bet)
            xhat[sp,:] = stat.gamma.ppf(prob[sp,:],a=Bet,scale=Alp) + Gam
            
        elif PDF[sp]==4: # v) Log-Gamma 3 parameters distribution.
            sigmay = np.sqrt(np.log(1 + (sigma[sp]/mu[sp])**2))
            muy = np.log(mu[sp])-(sigmay**2)/2
            Bety = (2/skewy[sp])**2
            Alpy = sigmay/np.sqrt(Bety)
            Gamy = muy-(Alpy*Bety)
            xhat[sp,:]  = np.exp(stat.gamma.ppf(prob[sp,:],a=Bety,scale=Alpy) + Gamy)
            
        elif PDF[sp]==5: # vi) Gumbel distribution.
            Sn = np.pi/np.sqrt(6)
            yn = 0.5772
            a = Sn/sigma[sp]
            u = mu[sp]-(yn/a)
            xhat[sp,:] = (u - np.log(-np.log(prob[sp,:]))/a)
    
        elif PDF[sp]==6: # vii) Exponential distribution.
            gamexp = mu[sp]-sigma[sp]
            xhat[sp,:] = (gamexp - (sigma[sp]*np.log(1 - prob[sp,:])))
    
    return xhat