# -*- coding: utf-8 -*-
"""
This script contains the main functions requiered by the methods implemented in
the climQMBC package, including the formatQM, the getDist, the getCDF, and the
getCDFinv functions.

Written by Sebastian Aedo Quililongo (1)
           Cristian Chadwick         (2)
           Fernando Gonzalez-Leiva   (1)
           Jorge Gironas             (1)
           
  (1) Pontificia Universidad Catolica de Chile, Santiago, Chile
      Department of Environmental and Hydraulic Engineering
  (2) Universidad Adolfo Ibanez, Santiago, Chile
      Faculty of Engineering and Sciences
Maintainer contact: slaedo@uc.cl
Revision: 0, updated Dec 2021
"""

import scipy.stats as stat
import numpy as np

def formatQM(obs,mod,frq):
    """
    This function formats the inputs and gets basic statistics for the 
    different Quantile Mapping (QM, DQM, QDM, UQM and SDM) methods available in
    the climQMBC package. If monthly data is specified, the input series will 
    be reshaped to a matrix of 12 rows and several columns equal to the number 
    of years of each series. If annual data is specified, the input is reshaped
    to a row vector with same entries as the input series.

    Description:
        0) Check if annual or monthly data is specified.
        
        1) Get number of years of the observed period.
        
        2) If monthly data is specified, reshape the input series to a matrix
           of 12 rows and several columns equal to the number of years of each
           series. If annual data is specified, reshape the input to a row
           vector with same entries as the input series.
        
        3) If monthly data is specified, get monthly mean, standard deviation,
           skewness, and log-skewness for the historical period of the observed
           and modeled series. If annual data is specified, get monthly mean, 
           standard deviation, skewness, and log-skewness for the historical 
           period of the observed and modeled series.

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

        frq:    A string specifying if the input is annual or monthly data. If
                not specified, it is set monthly as default.
                    Monthly:   frq = 'M'
                    Annual:    frq = 'A'

        NOTE: This routine considers that obs and mod series start in january
        of the first year and ends in december of the last year.


    Output:
        y_obs:          Number of observed years.
        
        obs_series:     A column vector of monthly or annual observed data
                        (temperature or precipitation). If monthly frequency is
                        specified, the length of this vector is 12 times the 
                        number of observed years [12, y_obs]. If annual 
                        frequency is specified, the length of this vector is 
                        equal to the number of observed years [1, y_obs].
        
        mod_series:     A column vector of monthly or annual modeled data 
                        (temperature or precipitation). If monthly frequency is
                        specified, the length of this vector is 12 times the 
                        number of observed years [12, y_mod]. If annual 
                        frequency is specified, the length of this vector is 
                        equal to the number of observed years [1, y_mod].
        
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
    
    # 0) Check if annually or monthly data is specified.
    if frq == 'A':
        I = 1
    else:
        I = 12
    
    # 1) Get number of years of the observed period.
    y_obs = int(obs.shape[0]/I)
    
    # 2) If monthly data is specified, reshape the input series to a matrix of
    #    12 rows and several columns equal to the number of years of each 
    #    series. If annually data is specified, reshape the input to a row 
    #    vector with same entries as the input series.
    obs_series = obs.reshape((I,int(obs.shape[0]/I)),order='F')
    mod_series = mod.reshape((I,int(mod.shape[0]/I)),order='F')
        
    # 3) If monthly data is specified, get monthly mean, standard deviation, 
    #    skewness, and log-skewness for the historical period of the observed
    #    and modeled series. If annually data is specified, get monthly mean,
    #    standard deviation, skewness, and log-skewness for the historical
    #    period of the observed and modeled series.
    mu_obs  = np.nanmean(obs_series,1)     # Mean
    sigma_obs = np.nanstd(obs_series,1,ddof=1)    # Standard deviation
    skew_obs = stat.skew(obs_series,1,bias=False)     # Skewness
    with np.errstate(divide = 'ignore'):
        Ln_obs  = np.log(obs_series)
    Ln_obs[np.imag(Ln_obs)!=0] = 0
    Ln_obs[np.isinf(Ln_obs)] = np.log(0.01)
    skewy_obs = stat.skew(Ln_obs,1,bias=False)        # Log-skewness    
    
    mu_mod  = np.nanmean(mod_series[:,:y_obs],1)     # Mean
    sigma_mod = np.nanstd(mod_series[:,:y_obs],1,ddof=1)    # Standard deviation
    skew_mod = stat.skew(mod_series[:,:y_obs],1,bias=False)     # Skewness
    with np.errstate(divide = 'ignore'):
        Ln_mod  = np.log(mod_series[:,:y_obs])
    Ln_mod[np.imag(Ln_mod)!=0] = 0
    Ln_mod[np.isinf(Ln_mod)] = np.log(0.01)
    skewy_mod = stat.skew(Ln_mod,1,bias=False)        # Log-skewness    

    return y_obs,obs_series,mod_series,mu_obs,sigma_obs,skew_obs,skewy_obs,mu_mod,sigma_mod,skew_mod,skewy_mod


def getDist(series,mu,sigma,skew,skewy,var):
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

    """
    
    # 1) Get the number of years to compute the empirical distribution in step 
    #    3) and get the number of rows of the input series.
    y_series = series.shape[1]
    n = series.shape[0]
    
    # 2) Initialize column vectors for the statistics needed for the available
    #    probability distribution functions.
    PDF    = np.zeros(n)
    sigmay = np.zeros(n)
    muy    = np.zeros(n)
    A      = np.zeros(n)
    B      = np.zeros(n)
    Alp    = np.zeros(n)
    Bet    = np.zeros(n)
    Gam    = np.zeros(n)
    Alpy   = np.zeros(n)
    Bety   = np.zeros(n)
    Gamy   = np.zeros(n)
    a      = np.zeros(n)
    u      = np.zeros(n)
    
    # 3) Perform the Kolmogorov-Smirnov test for each row.
    for m in range(n):
        
        # a) Get empirical distribution.
        sortdata = np.sort(series[m,:]) 
        probEmp = np.arange(1/(y_series+1),y_series/(y_series+1) + 1/(y_series+1),1/(y_series+1))
    
        # b) Compare distributions.
        # i) Normal distribution.
        normal = stat.norm.cdf(sortdata,mu[m],sigma[m])
        KSnormal = max(abs(normal-probEmp))
    
        # ii) Log Normal distribution.
        if (series[m,:] < 0).any():
            KSlognormal = 1
        else:
            sigmay[m]= np.sqrt(np.log(1+(sigma[m]/mu[m])**2))
            muy[m]   = np.log(mu[m])-(sigmay[m]**2)/2                            
            lognormal  = stat.lognorm.cdf(sortdata,scale=np.exp(muy[m]),s=sigmay[m])
            KSlognormal = max(abs(lognormal-probEmp))
        
        # iii) Gamma 2 parameters distribution.
        if (series[m,:] < 0).any():
            KSgammaII = 1
        else:
            A[m]   = (sigma[m]**2)/mu[m] 
            B[m]   = (mu[m]/sigma[m])**2
            GammaII  = stat.gamma.cdf(sortdata,a=B[m],scale=A[m])
            KSgammaII = max(abs(GammaII-probEmp))
    
        # iv) Gamma 3 parameters distribution.
        #     (Pearson 3 parameters distribution)
        Bet[m]  = (2/skew[m])**2
        Alp[m]  = sigma[m]/np.sqrt(Bet[m])
        Gam[m]  = mu[m]-(Alp[m]*Bet[m])
        GammaIII  = stat.gamma.cdf(sortdata-Gam[m],a=Bet[m],scale=Alp[m])
        KSgammaIII = max(abs(GammaIII-probEmp))
    
        # v) Log-Gamma 3 parameters distribution.
        #    (Log-Pearson 3 parameters distribution)
        if (series[m,:] < 0).any():
            KSLpIII = 1
        else:
            Bety[m] = (2/skewy[m])**2
            Alpy[m] = sigmay[m]/np.sqrt(Bety[m])
            Gamy[m] = muy[m]-(Alpy[m]*Bety[m])
            Lnsortdata = np.log(sortdata)
            Lnsortdata[np.isinf(Lnsortdata)] = np.log(0.01)
            LpIII  = stat.gamma.cdf(Lnsortdata-Gamy[m],a=Bety[m],scale=Alpy[m])            
            KSLpIII = max(abs(LpIII-probEmp))

        # vi) Gumbel distribution.
        Sn    = np.pi/np.sqrt(6)
        yn    = 0.5772
        a[m] = Sn/sigma[m]
        u[m] = mu[m]-(yn/a[m])
        gumbel = np.exp(-np.exp(-a[m]*(sortdata-u[m])))
        KSgumbel = max(abs(gumbel-probEmp))
    
        # vii) Exponential distribution.
        gamexp = mu[m]-sigma[m]
        exponential = 1-np.exp(-1/sigma[m]*(sortdata-gamexp))
        exponential[exponential<0] = 0
        KSexponential= max(abs(exponential-probEmp))
        
        # c) If variable is precipitation, set KS=1 to distributions that allow
        #    negative values (this will discard those distributions).
        if var==1:
            KSnormal=1
            KSgammaIII=1
            KSgumbel=1
            KSexponential=1
        
        # d) The distribution with lower KS value is considered for each month.
        bestPDF= np.argmin([KSnormal,KSlognormal,KSgammaII,KSgammaIII,KSLpIII,KSgumbel,KSexponential])
    
        PDF[m]= bestPDF
    return PDF


def getCDF(PDF,series,mu,sigma,skew,skewy):
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
        if PDF[m] == 0: # i) Normal distribution.
            Taot[m,:] = stat.norm.cdf(series[m,:],mu[m],sigma[m])
            
        elif PDF[m] == 1: # ii) Log-Normal distribution.
            sigmay= np.sqrt(np.log(1+(sigma[m]/mu[m])**2))
            muy   = np.log(mu[m])-(sigmay**2)/2
            Taot[m,:] = stat.lognorm.cdf(series[m,:],scale=np.exp(muy),s=sigmay)
            
        elif PDF[m] == 2: # iii) Gamma 2 parameters distribution.
            A   = (sigma[m]**2)/mu[m] 
            B   = (mu[m]/sigma[m])**2    
            Taot[m,:] = stat.gamma.cdf(series[m,:],a=B,scale=A)
            
        elif PDF[m] == 3: # iv) Gamma 3 parameters distribution.
            Bet  = (2/skew[m])**2
            Alp  = sigma[m]/np.sqrt(Bet)
            Gam  = mu[m]-(Alp*Bet)
            Taot[m,:]  = stat.gamma.cdf(series[m,:]-Gam,a=Bet,scale=Alp)
            
        elif PDF[m] == 4: # v) Log-Gamma 3 parameters distribution.
            sigmay= np.sqrt(np.log(1+(sigma[m]/mu[m])**2))
            muy   = np.log(mu[m])-(sigmay**2)/2
            Bety = (2/skewy[m])**2
            Alpy = sigmay/np.sqrt(Bety)
            Gamy = muy-(Alpy*Bety)
            Lnsortdata = np.log(series[m,:])
            Lnsortdata[np.isinf(Lnsortdata)] = np.log(0.01)
            Lnsortdata[np.imag(Lnsortdata)!=0] = np.log(0.01)
            Taot[m,:]  = stat.gamma.cdf(Lnsortdata-Gamy,a=Bety,scale=Alpy)
            
        elif PDF[m] == 5: # vi) Gumbel distribution.
            Sn    = np.pi/np.sqrt(6)
            yn    = 0.5772
            a = Sn/sigma[m]
            u = mu[m]-(yn/a)
            Taot[m,:] = np.exp(-np.exp(-a*(series[m,:]-u)))
    
        elif PDF[m] == 6: # vii) Exponential distribution.
            gamexp = mu[m]-sigma[m]
            exponential = 1-np.exp(-1/sigma[m]*(series[m,:]-gamexp))
            exponential[exponential<0] = 0
            Taot[m,:] = exponential
    
    # 3) Tail probabilities are set to (0 +  threshold) and (1 - threshold) to
    #    avoid numerical errors.
    th = 1e-3
    Taot[Taot>1-th] = 1-th
    Taot[Taot<th] = th
    
    return Taot


def getCDFinv(PDF,Taot,mu,sigma,skew,skewy):
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
        if PDF[m] == 0: # i) Normal distribution.
            xhat[m,:] = stat.norm.ppf(Taot[m,:],mu[m],sigma[m])
            
        elif PDF[m] == 1: # ii) Log-Normal distribution.
            sigmay= np.sqrt(np.log(1+(sigma[m]/mu[m])**2))
            muy   = np.log(mu[m])-(sigmay**2)/2
            xhat[m,:] = stat.lognorm.ppf(Taot[m,:],scale=np.exp(muy),s=sigmay)
            
        elif PDF[m] == 2: # iii) Gamma 2 parameters distribution.
            A   = (sigma[m]**2)/mu[m] 
            B   = (mu[m]/sigma[m])**2    
            xhat[m,:] = stat.gamma.ppf(Taot[m,:],a=B,scale=A)
            
        elif PDF[m] == 3: # iv) Gamma 3 parameters distribution.
            Bet  = (2/skew[m])**2
            Alp  = sigma[m]/np.sqrt(Bet)
            Gam  = mu[m]-(Alp*Bet)
            xhat[m,:]  = stat.gamma.ppf(Taot[m,:],a=Bet,scale=Alp) + Gam
            
        elif PDF[m] == 4: # v) Log-Gamma 3 parameters distribution.
            sigmay= np.sqrt(np.log(1+(sigma[m]/mu[m])**2))
            muy   = np.log(mu[m])-(sigmay**2)/2
            Bety = (2/skewy[m])**2
            Alpy = sigmay/np.sqrt(Bety)
            Gamy = muy-(Alpy*Bety)
            xhat[m,:]  = np.exp(stat.gamma.ppf(Taot[m,:],a=Bety,scale=Alpy) + Gamy)
            
        elif PDF[m] == 5: # vi) Gumbel distribution.
            Sn    = np.pi/np.sqrt(6)
            yn    = 0.5772
            a = Sn/sigma[m]
            u = mu[m]-(yn/a)
            xhat[m,:] = (u-np.log(-np.log(Taot[m,:]))/a)
    
        elif PDF[m] == 6: # vii) Exponential distribution.
            gamexp = mu[m]-sigma[m]
            xhat[m,:] = (gamexp-(sigma[m]*np.log(1-Taot[m,:])))
    
    return xhat

    
