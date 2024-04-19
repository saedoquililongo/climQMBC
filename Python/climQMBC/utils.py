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
                             allow_negatives = 1 or True: Allow negatives (default)
                             allow_negatives = 0 or False: Do not allow negative

        frq:             A string specifying if the input frequency is daily,
                         monthly or annual.
                             Daily:     frq = 'D'
                             Monthly:   frq = 'M'
                             Annual:    frq = 'A' (default)

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
        series:         An array of daily, monthly or annual observed data,
                        without considering leap days. Possible input dimensions:
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
                        Default: win = 1

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
        series_matrix:  A matrix of daily, monthly or annual observed data,
                        without considering leap days. For daily, monthly and
                        annual frequency, the number of rows will be 365, 12 and
                        1, respectively. The number of columns will be equal to
                        the number of years.
                        [1, 12 or 365, years]
                        
                        A 3D array of daily data, without considering leap days,
                        with a dimension that considers a centered moving window
                        for each day of the year.
                        [365, 2*win-1, years]

        projected_win:  An integer indicating how many days to consider 
                        backwards and forward to get the statistics of each 
                        calendar day.The length of the window will be 
                        (2*win_day-1). For example, day_win=15 -> window of 29.
                        Default: win = 1
                            
        frq:             A string specifying if the input frequency is daily,
                         monthly or annual.
                             Daily:     frq = 'D'
                             Monthly:   frq = 'M'
                             Annual:    frq = 'A' (default)

    Output:
        win_series:     An array of daily, monthly or annual observed data,
                        without considering leap days. Possible input dimensions:
                        2D array:
                            [sub-periods, years]
                            [sub-periods + days window, years]
                        3D array:
                            [sub_periods, projected periods window, years]
                            [sub-periods + days window, projected periods window, years]
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


def set_norain_to_nan(series_moving, pp_factor, pp_threshold):
    min_rainday = 30
    # Minimos necesarios  -> min_rainday
    # Dias que tengo      -> rainday_count[per]
    # Valores a modificar -> replace_values_nans
    # Valores no-nans     -> max(Minimos_necesarios - dias que tengo,0)
    #                            >0: Faltan valores  ->  No todos son nans
    #                            <0: Sobran valores  ->  Todos son nans
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
    Kolmogorov-Smirnov (KS) test. If the series consider monthly data, it will
    have 12 rows and each row will represent a month. For annual data the 
    series will have only one row. Only strictly positive distributions are 
    considered for precipitation and strictly positive distributions are 
    discarded if the series has negative values. The available distributions 
    are:
        1) Normal distribution
        2) Log-Normal distribution
        3) Gamma 2 parameters distribution
        4) Gamma 3 parameters distribution
           (Pearson 3 parameters distribution)
        5) Log-Gamma 3 parameters distribution
           (Log-Pearson 3 parameters distribution)
        6) Gumbel distribution
        7) Exponential distribution
        
        For precipitation, only 2), 3) and 5) are considered (1, 4, 6, and 7 
        are discarded). For series with negative values, only 1), 3), 4), 6),
        and 7) are considered (2, 3 and 5 are discarded).

    Input:
        series:     A matrix of monthly or annual data (temperature or 
                    precipitation). If the series consider monthly data, it 
                    will have 12 rows and each row will represent a month. 
                    For annual data the series will have only one row.
        
        mu:         A column vector of mean values of the series. [12,1] if the
                    series consider monthly data and [1,1] if the series 
                    consider annual data.
        
        sigma:      A column vector of standard deviation of the series. [12,1]
                    if the series consider monthly data and [1,1] if the series
                    consider annual data.
        
        skew:       A column vector of skewness of the series. [12,1] if the 
                    series consider monthly data and [1,1] if the series 
                    consider annual data.
        
        skewy:      A column vector of skewness of the logarithm of the series.
                    [12,1] if the series consider monthly data and [1,1] if the
                    series consider annual data.
        
        var:        A flag that identifies if data are temperature or
                    precipitation. This flag tells the getDist function if it
                    has to discard distribution functions that allow negative
                    numbers.
                        Temperature:   var = 0
                        Precipitation: var = 1

    Output:
        PDF:    A column vector with an ID for the resulting distribution from
                the KS test. [12,1] if the series consider monthly data and 
                [1,1] if the series consider annual data. The ID is related to 
                the numeration of the distribution listed in the description of
                this function. This ID is used in the getCDF and getCDFinv 
                functions of the climQMBC package.
                
        mu_obs:         If monthly frequency is specified, a column vector of 
                        monthly mean of observed data [12,1]. If annual 
                        frequency is specified, the mean of the observed data 
                        (float).
        
        mu_mod:         If monthly frequency is specified, a column vector of 
                        monthly mean of modeled data of the historical period 
                        [12,1]. If annual frequency is specified, the mean of 
                        the modeled data of the historical period(float).
        
        sigma_obs:      If monthly frequency is specified, a column vector of 
                        monthly standard deviation of observed data [12,1]. If
                        annual frequency is specified, the standard deviation 
                        of the observed data (float).
        
        sigma_mod:      If monthly frequency is specified, a column vector of 
                        monthly standard deviation of modeled data of the 
                        historical period [12,1]. If annual frequency is 
                        specified, the standard deviation of the modeled data 
                        of the historical period(float).
        
        skew_obs:       If monthly frequency is specified, a column vector of 
                        monthly skewness of observed data [12,1]. If annual 
                        frequency is specified, the skewness of the observed 
                        data (float).
        
        skew_mod:       If monthly frequency is specified, a column vector of 
                        monthly skewness of modeled data of the historical 
                        period [12,1]. If annual frequency is specified, the 
                        skewness of the modeled data of the historical period
                        (float).
        
        skewy_obs:      If monthly frequency is specified, a column vector of 
                        monthly skewness of the logarithm of observed data 
                        [12,1]. If annual frequency is specified, the skewness
                        of the logarithm of the observed data (float).
        
        skewy_mod:      If monthly frequency is specified, a column vector of 
                        monthly skewness of the logarithm of modeled data of 
                        the historical period [12,1]. If annual frequency is 
                        specified, the skewness of the logarithm of the modeled
                        data of the historical period(float).
    """
    
    # 1) Get the number of years to compute the empirical distribution in step 
    #    3) and get the number of rows of the input series.
    y_series = series.shape[1]
    n = series.shape[0]
    
    # 2) Initialize column vectors for the statistics needed for the available
    #    probability distribution functions.
    PDF = np.zeros(n)
    ks_fail = 0
    
    # 3) Perform the Kolmogorov-Smirnov test for each row.
    for m in range(n):
        series_sub = series[m,:]
        
        ######
        series_sub = series_sub[~np.isnan(series_sub)]
        y_series = len(series_sub)
        #####
        
        # a) Get empirical distribution.
        sortdata = np.sort(series_sub) 
        probEmp = np.arange(1/(y_series+1),
                            y_series/(y_series+1) + 1/(y_series+1),
                            1/(y_series+1))
        
        probEmp = np.linspace(1/(y_series+1),
                              y_series/(y_series+1) + 1/(y_series+1),
                              y_series)

        # b) Compare distributions.
        # i) Normal distribution.
        normal = stat.norm.cdf(sortdata, mu[m], sigma[m])
        KSnormal = max(abs(normal-probEmp))

        # ii) Log Normal distribution.
        if (series_sub<0).any():
            KSlognormal = 1
        else:
            sigmay = np.sqrt(np.log(1 + (sigma[m]/mu[m])**2))
            muy = np.log(mu[m]) - (sigmay**2)/2
            lognormal = stat.lognorm.cdf(sortdata, scale=np.exp(muy), s=sigmay)
            KSlognormal = max(abs(lognormal-probEmp))
        
        # iii) Gamma 2 parameters distribution.
        if (series_sub<0).any():
            KSgammaII = 1
        else:
            A = (sigma[m]**2)/mu[m] 
            B = (mu[m]/sigma[m])**2
            GammaII = stat.gamma.cdf(sortdata, a=B, scale=A)
            KSgammaII = max(abs(GammaII-probEmp))
    
        # iv) Gamma 3 parameters distribution.
        #     (Pearson 3 parameters distribution)
        Bet = (2/skew[m])**2
        Alp = sigma[m]/np.sqrt(Bet)
        Gam = mu[m] - (Alp*Bet)
        GammaIII = stat.gamma.cdf(sortdata-Gam, a=Bet, scale=Alp)
        KSgammaIII = max(abs(GammaIII-probEmp))
    
        # v) Log-Gamma 3 parameters distribution.
        #    (Log-Pearson 3 parameters distribution)
        if (series_sub<0).any():
            KSLpIII = 1
        else:
            sigmay = np.sqrt(np.log(1 + (sigma[m]/mu[m])**2))
            Bety = (2/skewy[m])**2
            Alpy = sigmay/np.sqrt(Bety)
            Gamy = muy-(Alpy*Bety)
            sortdata_log = np.log(sortdata)
            sortdata_log[np.isinf(sortdata_log)] = np.log(np.random.rand(np.isinf(sortdata_log).sum())*0.01)
            LpIII = stat.gamma.cdf(sortdata_log-Gamy, a=Bety, scale=Alpy)
            KSLpIII = max(abs(LpIII-probEmp))

        # vi) Gumbel distribution.
        Sn = np.pi/np.sqrt(6)
        yn = 0.5772
        a = Sn/sigma[m]
        u = mu[m] - (yn/a)
        gumbel = np.exp(-np.exp(-a*(sortdata - u)))
        KSgumbel = max(abs(gumbel-probEmp))
    
        # vii) Exponential distribution.
        gamexp = mu[m] - sigma[m]
        exponential = 1-np.exp(-1/sigma[m]*(sortdata - gamexp))
        exponential[exponential<0] = 0
        KSexponential= max(abs(exponential-probEmp))
        
        # c) If variable is precipitation, set KS=1 to distributions that allow
        #    negative values (this will discard those distributions).        
        if not allow_negatives:
            KSnormal = 1
            KSgammaIII = 1
            KSgumbel = 1
            KSexponential = 1

        # d) The distribution with lower KS value is considered for each month.
        KS_vals = [KSnormal, KSlognormal, KSgammaII, KSgammaIII, KSLpIII, KSgumbel, KSexponential]
        bestPDF = np.argmin(KS_vals)
        
        ks_crit = stat.ksone.ppf(1-0.05/2, y_series)
        ks_fail = ks_fail + (np.sign(min(KS_vals)-ks_crit)+1)/2

        PDF[m] = bestPDF
    
    # # Check if ks-test failures
    # if ks_fail>0:
    #     print(f'KS-Test failed: {int(ks_fail)} times out of {int(n)}')
        
    return PDF


def getCDF(PDF, series, mu, sigma, skew, skewy):
    """
    This function evaluates each row of the series in the respective cumulative
    distribution function assigned by the Kolmogorov-Smirnov (KS) test in the 
    getDist function of the climQMBC package. The available distributions are:
        1) Normal distribution
        2) Log-Normal distribution
        3) Gamma 2 parameters distribution
        4) Gamma 3 parameters distribution
           (Pearson 3 parameters distribution)
        5) Log-Gamma 3 parameters distribution
           (Log-Pearson 3 parameters distribution)
        6) Gumbel distribution
        7) Exponential distribution

    Note that each row is associated to an independent distribution function.

    Input:
        PDF:        A column vector with an ID for the resulting distribution 
                    from the KS test. [12,1] if the series consider monthly 
                    data and [1,1] if the series consider annual data. The ID 
                    is related to the numeration of the distribution listed in 
                    the description of this function. 
        
        series:     A matrix of monthly or annual data (temperature or 
                    precipitation). If the series consider monthly data, it 
                    will have 12 rows and each row will represent a month. For
                    annual data the series will have only one row.
        
        mu:         A column vector of mean values of the series. [12,1] if the
                    series consider monthly data and [1,1] if the series
                    consider annual data.
        
        sigma:      A column vector of standard deviation of the series. [12,1]
                    if the series consider monthly data and [1,1] if the series
                    consider annual data.
        
        skew:       A column vector of skewness of the series. [12,1] if the 
                    series consider monthly data and [1,1] if the series 
                    consider annual data.
        
        skewy:      A column vector of skewness of the logarithm of the series.
                    [12,1] if the series consider monthly data and [1,1] if the
                    series consider annual data.

    Output:
        Taot:   A column vector with the non-exceedance probability for each
                row of the input series. [12,1] if the series consider monthly
                data and [1,1] if the series consider annual data.

    """

    # 1) Get the number of rows and years of the series.
    n_m = series.shape[0]
    n_y = series.shape[1]
    
    # 2) Compute the cumulative distribution function to the values of each 
    #    row, based on the distribution assigned in the getDist function of the
    #    climQMBC package.
    Taot = np.zeros((n_m,n_y))
    for m in range (n_m):
        series_sub = series[m,:]
        
        if PDF[m]==0: # i) Normal distribution.
            Taot[m,:] = stat.norm.cdf(series_sub, mu[m], sigma[m])
            
        elif PDF[m]==1: # ii) Log-Normal distribution.
            sigmay = np.sqrt(np.log(1 + (sigma[m]/mu[m])**2))
            muy   = np.log(mu[m])-(sigmay**2)/2
            Taot[m,:] = stat.lognorm.cdf(series_sub, scale=np.exp(muy), s=sigmay)
            
        elif PDF[m]==2: # iii) Gamma 2 parameters distribution.
            A = (sigma[m]**2)/mu[m] 
            B = (mu[m]/sigma[m])**2
            Taot[m,:] = stat.gamma.cdf(series_sub, a=B, scale=A)
            
        elif PDF[m]==3: # iv) Gamma 3 parameters distribution.
            Bet = (2/skew[m])**2
            Alp = sigma[m]/np.sqrt(Bet)
            Gam = mu[m]-(Alp*Bet)
            Taot[m,:] = stat.gamma.cdf(series_sub-Gam, a=Bet, scale=Alp)
            
        elif PDF[m]==4: # v) Log-Gamma 3 parameters distribution.
            sigmay = np.sqrt(np.log(1 + (sigma[m]/mu[m])**2))
            
            muy = np.log(mu[m]) - (sigmay**2)/2
            Bety = (2/skewy[m])**2
            Alpy = sigmay/np.sqrt(Bety)
            Gamy = muy - (Alpy*Bety)
            Lnsortdata = np.log(series_sub)
            Lnsortdata[np.isinf(Lnsortdata)] = np.log(np.random.rand(np.isinf(Lnsortdata).sum())*0.01)
            Taot[m,:]  = stat.gamma.cdf(Lnsortdata-Gamy,a=Bety,scale=Alpy)
            
        elif PDF[m]==5: # vi) Gumbel distribution.
            Sn = np.pi/np.sqrt(6)
            yn = 0.5772
            a = Sn/sigma[m]
            u = mu[m]-(yn/a)
            Taot[m,:] = np.exp(-np.exp(-a*(series_sub - u)))
    
        elif PDF[m]==6: # vii) Exponential distribution.
            gamexp = mu[m]-sigma[m]
            exponential = 1-np.exp(-1/sigma[m]*(series_sub-gamexp))
            exponential[exponential<0] = 0
            Taot[m,:] = exponential
    
    # 3) Tail probabilities are set to (0 +  threshold) and (1 - threshold) to
    #    avoid numerical errors.
    th = 1e-3
    Taot[Taot>1-th] = 1-th
    Taot[Taot<th] = th
    
    return Taot


def getCDFinv(PDF, Taot, mu, sigma, skew, skewy):
    """
    This function evaluates the probability, Taot, in the respective inverse
    cumulative distribution function assigned by the Kolmogorov-Smirnov (KS)
    test in the getDist function of the climQMBC package. The available
    distributions are:
        1) Normal distribution
        2) Log-Normal distribution
        3) Gamma 2 parameters distribution
        4) Gamma 3 parameters distribution
           (Pearson 3 parameters distribution)
        5) Log-Gamma 3 parameters distribution
           (Log-Pearson 3 parameters distribution)
        6) Gumbel distribution
        7) Exponential distribution
        
    Note that each row is associated to an independent distribution function.
        
    Input:
        PDF:    A column vector with an ID for the resulting distribution from
                the KS test. [12,1] if the series consider monthly data and 
                [1,1] if the series consider annual data. The ID is related to
                the numeration of the distribution listed in the description 
                of this function. 
        
        Taot:   A column vector with the non-exceedance probability for each 
                row of the input series. [12,1] if the series consider monthly 
                data and [1,1] if the series consider annual data.
        
        mu:     A column vector of mean values of the series. [12,1] if the 
                series consider monthly data and [1,1] if the series consider 
                annual data.
        
        sigma:  A column vector of standard deviation of the series. [12,1] if 
                the series consider monthly data and [1,1] if the series 
                consider annual data.
        
        skew:   A column vector of skewness of the series. [12,1] if the serie
                s consider monthly data and [1,1] if the series consider annual
                data.
        
        skewy:  A column vector of skewness of the logarithm of the series. 
                [12,1] if the series consider monthly data and [1,1] if the 
                series consider annual data.
            
    Output:
        xhat:   A column vector with the values obtained when the inverse
                cumulative distribution function is applied. [12,1] if the
                series consider monthly data and [1,1] if the series consider
                annual data.
    """

    # 1) Get the number of rows and years of the series.
    n_m = Taot.shape[0]
    n_y = Taot.shape[1]
    
    # 2) Compute the cumulative distribution function to the values of each 
    #    row, based on the distribution assigned in the getDist function of the
    #    climQMBC package.
    xhat = np.zeros((n_m,n_y))
    for m in range (n_m):
        
        if PDF[m]==0: # i) Normal distribution.
            xhat[m,:] = stat.norm.ppf(Taot[m,:],mu[m],sigma[m])
            
        elif PDF[m]==1: # ii) Log-Normal distribution.
            sigmay = np.sqrt(np.log(1 + (sigma[m]/mu[m])**2))
            muy = np.log(mu[m])-(sigmay**2)/2
            xhat[m,:] = stat.lognorm.ppf(Taot[m,:],scale=np.exp(muy),s=sigmay)
            
        elif PDF[m]==2: # iii) Gamma 2 parameters distribution.
            A = (sigma[m]**2)/mu[m] 
            B = (mu[m]/sigma[m])**2    
            xhat[m,:] = stat.gamma.ppf(Taot[m,:],a=B,scale=A)
            
        elif PDF[m]==3: # iv) Gamma 3 parameters distribution.
            Bet = (2/skew[m])**2
            Alp = sigma[m]/np.sqrt(Bet)
            Gam = mu[m]-(Alp*Bet)
            xhat[m,:] = stat.gamma.ppf(Taot[m,:],a=Bet,scale=Alp) + Gam
            
        elif PDF[m]==4: # v) Log-Gamma 3 parameters distribution.
            sigmay = np.sqrt(np.log(1 + (sigma[m]/mu[m])**2))
            muy = np.log(mu[m])-(sigmay**2)/2
            Bety = (2/skewy[m])**2
            Alpy = sigmay/np.sqrt(Bety)
            Gamy = muy-(Alpy*Bety)
            xhat[m,:]  = np.exp(stat.gamma.ppf(Taot[m,:],a=Bety,scale=Alpy) + Gamy)
            
        elif PDF[m]==5: # vi) Gumbel distribution.
            Sn = np.pi/np.sqrt(6)
            yn = 0.5772
            a = Sn/sigma[m]
            u = mu[m]-(yn/a)
            xhat[m,:] = (u - np.log(-np.log(Taot[m,:]))/a)
    
        elif PDF[m]==6: # vii) Exponential distribution.
            gamexp = mu[m]-sigma[m]
            xhat[m,:] = (gamexp - (sigma[m]*np.log(1 - Taot[m,:])))
    
    return xhat