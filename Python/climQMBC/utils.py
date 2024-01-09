# -*- coding: utf-8 -*-
"""
This script contains the main functions requiered by the methods implemented in
the climQMBC package, including the formatQM, the getStats, the getDist,
the getCDF, and the getCDFinv functions.

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
Revision: 2, updated Jan 2024
"""


import scipy.stats as stat
import numpy as np


def formatQM(series_, allow_negatives, frq, pp_threshold, pp_factor):
    """
    This function formats a time series and gets basic statistics for the 
    different Quantile Mapping methods (QM, DQM, QDM, UQM and SDM) available in
    the climQMBC package. If monthly data is specified, the input series will 
    be reshaped to a matrix of 12 rows and several columns equal to the number 
    of years of the series. If annual data is specified, the input is reshaped
    to a row vector with same entries as the input series. For precipitation,
    determined by allow_negatives=1, physically null values (values below
    pp_threshold) are replaced by random positive values below pp_factor.

    Description:
        0) Check frequency.

        1) If variable is precipitation, determined by allow_negatives=1, 
           replace low values with random values.

        2) Get number of years of the series.

        3) If monthly data is specified, reshape the input series to a matrix
           of 12 rows and several columns equal to the number of years of each
           series. If annual data is specified, reshape the input to a row
           vector with same entries as the input series.

    Inputs:
        series_:         A column vector of monthly or annual data. If monthly
                         frequency is specified, the length of this vector
                         should be 12 times the number of years in the series
                         [12 x years, 1]. If annual frequency is specified, the
                         length of this vector should be equal to the number of
                         years of the series [years, 1].

        allow_negatives: A flag that identifies if the series allow for negative
                         values or not. If negative values are not allowed, the
                         series will be treated as a precipitation series,
                         replacing the values considered as physicaly null 
                         precipitation by near-zero random values.
                             allow_negatives = 1   : Allow negative values
                             allow_negatives = 0   : Does not allow negative
                                                     values

        frq:             A string specifying if the input has annual or monthly
                         data. If not specified, it is set annual as default.
                             Monthly:   frq = 'M'
                             Annual:    frq = 'A'

        pp_threshold:    A float indicating the threshold to consider 
                         physically null precipitation values.

        pp_factor:       A float indicating the maximum value of the random
                         values that replace physically null precipitation 
                         values.


    Output:
        years:          Number of years in the series.

        series    :     A matrix of monthly or annual data. If monthly frequency
                        is specified, the number of rows of the matrix is 12
                        and the number of columns is the number of years
                        [12, years]. If annual frequency is specified, the
                        number of rows of the matrix is 1 and the number of
                        columns is the number of years [1, years].
    """

    # This prevents modyfing the original input passed by reference
    ## Must look for a more elegant way to skip this step
    series = series_.copy()

    # 0) Check frequency.
    if frq=='D':
        I = 365
    elif frq=='M':
        I = 12
    elif frq=='A':
        I = 1
    else:
        I = 1

    # 1) If variable is precipitation, determined by allow_negatives=1, replace
    #    low values with random values.
    if not allow_negatives:
        bool_low = series<pp_threshold
        series[bool_low] = np.random.rand(bool_low.sum())*pp_factor

    # 2) Get number of years of the series.
    years = int(series.shape[0]/I)
    
    # 3) If monthly data is specified, reshape the input series to a matrix of
    #    12 rows and several columns equal to the number of years of each
    #    series. If annual data is specified, reshape the input to a row vector
    #    with same entries as the input series.
    series = series.reshape((I,int(series.shape[0]/I)), order='F')

    return years, series


def getStats(series):
    """

    """
    # 4) If monthly data is specified, get monthly mean, standard deviation and 
    #    skewness for the historical period of the observed
    #    and modeled series. If annually data is specified, get monthly mean,
    #    standard deviation, skewness, and log-skewness for the historical
    #    period of the observed and modeled series.
    mu  = np.nanmean(series, 1)     # Mean
    sigma = np.nanstd(series, 1, ddof=1)    # Standard deviation
    skew = stat.skew(series, 1, bias=False)     # Skewness
    series_log  = np.log(series)
    series_log[np.isinf(series_log)] = np.log(0.01)
    skewy = stat.skew(series_log, 1, bias=False)

    return mu, sigma, skew, skewy


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

    Description:
        1) Get the number of years to compute the empirical distribution in 
           step 3) and get the number of rows of the input series.
                                                              
        4) If monthly data is specified, get monthly mean, standard deviation,
           skewness, and log-skewness for the historical period of the observed
           and modeled series. If annual data is specified, get monthly mean, 
           standard deviation, skewness, and log-skewness for the historical 
           period of the observed and modeled series.

        2) Initialize column vectors for the statistics needed for the 
           available probability distribution functions.

        3) Perform the Kolmogorov-Smirnov test for each row.
            a) Get empirical distribution.
            b) Compare distributions.
            c) If variable is precipitation, discard distribution that allows
               negative values.
            d) Compare KS values.


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
        
        # a) Get empirical distribution.
        sortdata = np.sort(series_sub) 
        probEmp = np.arange(1/(y_series+1),
                            y_series/(y_series+1) + 1/(y_series+1),
                            1/(y_series+1))

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
            # Lnsortdata[np.isinf(Lnsortdata)] = np.log(0.01)
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
    
    # Check if ks-test failures
    if ks_fail>0:
        print(f'KS-Test failed: {int(ks_fail)} times out of {int(n)}')
        
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

    Description:
        1) Get the number of rows and years of the series.
        
        2) Compute the cumulative distribution function to the values of each 
           row, based on the distribution assigned in the getDist function of 
           the climQMBC package.
        
        3) Tail probabilities are set to (0 + threshold) and (1 - threshold) 
           to avoid numerical errors. Threshold is set in the last lines of 
           this function. Default is th = 10^-3.

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
            Lnsortdata[np.isinf(Lnsortdata)] = np.log(0.01)
            Lnsortdata[np.imag(Lnsortdata)!=0] = np.log(0.01)
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

    Description:
        1) Get the number of rows and years of the series.
    
        2) Compute the inverse cumulative distribution function to the values
           of each row, based on the distribution assigned in the getDist 
           function of the climQMBC package.
        
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

    
