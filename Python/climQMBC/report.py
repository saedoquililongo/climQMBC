# -*- coding: utf-8 -*-
import matplotlib.pylab as plt
import pandas as pd
import numpy as np

def report(obs, mod, var, fun=['QM','DQM','QDM','UQM','SDM'], y_init=0, y_wind=0, user_pdf=False, pdf_obs=None, pdf_mod=None):
    """
    This function generates two report of the performance of the different
    methods (QM, DQM, QDM, UQM and SDM) available in the climQMBC package.
    These reports are based on the mean and standard deviation of the series
    in the historical and future projected periods.

    The first report is a summary table with the overall performance of the QM
    method in the historical period and of the different methods in future
    projected periods. The methods performance is addressed by comparing its
    variations in the mean and standard deviation with the ones of the modeled
    data. 

    The second report consist of three figures. Figure 1 shows the cumulative 
    distribution of the observed, modeled, and corrected series in the 
    historical and future period. This figure also shows the complete time 
    series. Figure 2 and 3 shows the monthly mean and standard deviation, 
    respectively, of each series in the historical and future projected periods.
    In these two figures, in the future projected periods, the observed series 
    is replaced by an objective series which is computed as the observed series
    (monthly mean or standard deviation), scaled by the variation between the 
    projected and the historical period of the modeled series. Each projected 
    period is centered in a moving window whose length is equal to the length 
    of the historical period.

    NOTE: 
        1) Even if a set of methods are specified in the optional inputs, this
           function returns all five methods implemented in the climQMBC
           package.
    
        2) This report function was built for monthly data.
    
    Description:
        0) Get the number of observed and modeled years.
        
        1) Set non-declared arguments.
        
        2) Apply bias correction methods.
        
        3) Get observed, modeled and bias corrected statistics of the 
           historical and complete future period.
        
        4) Get observed, modeled and bias corrected statistics of the projected
           periods.
        
        5) Get delta mean and delta standard deviation.
        
        6) Display report:
            a) Report description.
            b) Table.
            c) Figures.
            
    Input:
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
    
    
    Optional inputs:
        fun:    A list of strings with the desired bias correction methods to
        be reported. If this input is not recieved by the function, all bias
        correction methods available in the climQMBC package will be reported.
        The methods supported are:
            a) 'QM' : Quantile Mapping
            b) 'DQM': Detrended Quantile Mapping
            c) 'QDM': Quantile Delta Mapping
            d) 'UQM': Unbiased Quantile Mapping
            Default: fun = ['QM','DQM','QDM','UQM','SDM']
    
        y_init:     First year of the observed and modeled series (integer).
                    (Default: y_init = 0)
    
        y_wind:     A list of integers with the year of the center of the
                    projected periods to be reported.
                    Default: y_wind = [int(y_obs+y_obs/2), int(y_mod-y_obs/2)]
                    This value sets a first projected period just after the end 
                    of historical period, and a second projected period just 
                    before the end of the modeled series.
    
    Output:
        QM_series:  A column vector of monthly or annual modeled data 
                    (temperature or precipitation) corrected by the QM method. 
                    If monthly frequency is specified, the length of this 
                    vector is 12 times the number of observed years 
                    [12 x y_mod, 1]. If annual frequency is specified, the 
                    length of this vector is equal to the number of observed 
                    years [y_mod, 1].
                    
        DQM_series: A column vector of monthly or annual modeled data 
                    (temperature or precipitation) corrected by the DQM method. 
                    If monthly frequency is specified, the length of this 
                    vector is 12 times the number of observed years 
                    [12 x y_mod, 1]. If annual frequency is specified, the 
                    length of this vector is equal to the number of observed 
                    years [y_mod, 1].
                    
        QDM_series: A column vector of monthly or annual modeled data 
                    (temperature or precipitation) corrected by the QDM method. 
                    If monthly frequency is specified, the length of this 
                    vector is 12 times the number of observed years 
                    [12 x y_mod, 1]. If annual frequency is specified, the 
                    length of this vector is equal to the number of observed 
                    years [y_mod, 1].
                    
        UQM_series: A column vector of monthly or annual modeled data 
                    (temperature or precipitation) corrected by the UQM method. 
                    If monthly frequency is specified, the length of this 
                    vector is 12 times the number of observed years 
                    [12 x y_mod, 1]. If annual frequency is specified, the 
                    length of this vector is equal to the number of observed 
                    years [y_mod, 1].
                    
        SDM_series: A column vector of monthly or annual modeled data 
                    (temperature or precipitation) corrected by the SDM method. 
                    If monthly frequency is specified, the length of this 
                    vector is 12 times the number of observed years 
                    [12 x y_mod, 1]. If annual frequency is specified, the 
                    length of this vector is equal to the number of observed 
                    years [y_mod, 1].
    
        NOTE: This function returns all five bias correction methods,
              independently of which methods are specified for this report.
    
    
    Written by Sebastian Aedo Quililongo (1*)
               Cristian Chadwick         (2)
               Fernando Gonzalez-Leiva   (3)
               Jorge Gironas             (3)
               
      (1) Centro de Cambio Global UC, Pontificia Universidad Catolica de Chile,
          Santiago, Chile
      (2) Faculty of Engineering and Sciences, Universidad Adolfo Ibanez,
          Santiago, Chile
      (3) Department of Hydraulics and Environmental Engineering, Pontificia
          Universidad Catolica de Chile, Santiago, Chile
          
    *Maintainer contact: slaedo@uc.cl
    Revision: 0, updated Dec 2021
    """
    
    from .methods import QM, DQM, QDM, UQM, SDM
    
    
    def var_op(a,b,var):
        if var==0:
            return a - b
        else:
            return a/b
    
    def ecdf(x):
        return [(1 + i)/len(x) for i in range(len(x))]


    # 0) Get the number of observed and modeled years
    n_obs = obs.shape[0]
    y_obs = int(n_obs/12)
    
    n_mod = mod.shape[0]
    y_mod = int(n_mod/12)
    
    # 1) Set headers and indexes
    df_ind = ['Modeled','QM','DQM','QDM','UQM','SDM']
    df_m = pd.DataFrame(index=df_ind)
    df_s = pd.DataFrame(index=df_ind)
    
    if y_wind==0:
        y_wind = np.array([int(y_obs+y_obs/2), int(y_mod-y_obs/2)])
        w_label = ['First period', 'Last period']
        rep_head = ['First', 'Last']
    else:
        w_label = np.array(y_wind)
        rep_head = np.array(y_wind)
        y_wind = np.array(y_wind) - y_init
        
    fun = ['Modeled']+fun
    
    # 2) Apply QM methods
    QM_series = QM(obs, mod, var, user_pdf=user_pdf, pdf_obs=pdf_obs, pdf_mod=pdf_mod)
    DQM_series = DQM(obs, mod, var, user_pdf=user_pdf, pdf_obs=pdf_obs, pdf_mod=pdf_mod)
    QDM_series = QDM(obs, mod, var, user_pdf=user_pdf, pdf_obs=pdf_obs, pdf_mod=pdf_mod)
    UQM_series = UQM(obs, mod, var, user_pdf=user_pdf, pdf_obs=pdf_obs, pdf_mod=pdf_mod)
    SDM_series = SDM(obs, mod, var)
    
    # 3) Get observed, modeled and bias corrected statistics
    # a) Get obs, mod, and QMs as [12 x n] matrix
    obs_M = obs.reshape((12, y_obs), order='F')
    mod_M = mod.reshape((12, y_mod), order='F')
    QM_M = QM_series.reshape((12, y_mod), order='F')
    DQM_M = DQM_series.reshape((12, y_mod), order='F')
    QDM_M = QDM_series.reshape((12, y_mod), order='F')
    UQM_M = UQM_series.reshape((12, y_mod), order='F')
    SDM_M = SDM_series.reshape((12, y_mod), order='F')
    
    ## s: List of historical and future periods
    ## v: List of modeled and bias corrected series
    s = [obs_M, mod_M[:,:y_obs], mod_M[:,y_obs:], QM_M[:,:y_obs], QM_M[:,y_obs:],
         DQM_M[:,y_obs:], QDM_M[:,y_obs:], UQM_M[:,y_obs:], SDM_M[:,y_obs:]]
    v = [mod_M, QM_M, DQM_M, QDM_M, UQM_M, SDM_M]
    mu = []
    std = []
    mu_M = []
    std_M = []
    for i in s:
        # b) Mean
        mu.append(np.nanmean(i))
        
        # c) Monthly mean
        mu_M.append(np.nanmean(i, 1))
        
        # d) Standard deviation
        std.append(np.nanstd(i, ddof=1))
        
        # e) Monthly standard deviation
        std_M.append(np.nanstd(i, 1, ddof=1))
        
    # 4) Get observed, modeled and bias corrected statistics of the projected
    #    periods
    mu_w = []
    std_w = []
    mu_Mw = []
    std_Mw = []
    for w in range(len(y_wind)):
        mu_w.append([])
        mu_Mw.append([])
        std_w.append([])
        std_Mw.append([])
        for i in v:
            # a) Projected period mean
            mu_w[w].append(np.nanmean(i[:,int(y_wind[w] - y_obs/2)-1:int(y_wind[w] + y_obs/2)]))
            
            # b) Projected period monthly mean
            mu_Mw[w].append(np.nanmean(i[:,int(y_wind[w] - y_obs/2)-1:int(y_wind[w] + y_obs/2)],1))
            
            # c) Projected period standard deviation
            std_w[w].append(np.nanstd(i[:,int(y_wind[w] - y_obs/2)-1:int(y_wind[w] + y_obs/2)],ddof=1))
            
            # d) Projected period monthly standard deviation
            std_Mw[w].append(np.nanstd(i[:,int(y_wind[w] - y_obs/2)-1:int(y_wind[w] + y_obs/2)],1,ddof=1))
        
    # 5) Get delta mean and delta standard deviation
    SF_m = []
    SF_s = []
    SF_m.append(var_op(mu_M[2],mu_M[1],var))
    SF_s.append(var_op(std_M[2],std_M[1],var))
    dm = []
    ds = []
    dm.append(var_op(mu[3],mu[0],var))
    dm.append(var_op(mu[4],mu[0],var))
    dm.append(var_op(mu[2],mu[1],var))
    ds.append(var_op(std[3],std[0],var))
    ds.append(var_op(std[4],std[0],var))
    ds.append(var_op(std[2],std[1],var))
    for i in range(4,len(s)):
        dm.append(var_op(mu[i],mu[3],var))
        ds.append(var_op(std[i],std[3],var))
    
    dm_w = []
    ds_w = []
    for w in range(len(y_wind)):
        SF_m.append(var_op(mu_Mw[w][0],mu_M[1],var))
        SF_s.append(var_op(std_Mw[w][0],std_M[1],var))
        dm_w.append([])
        ds_w.append([])
        dm_w[w].append(var_op(mu_w[w][0],mu[1],var))
        ds_w[w].append(var_op(std_w[w][0],std[1],var))
        for i in range(1,len(v)):
            dm_w[w].append(var_op(mu_w[w][i],mu[3],var))
            ds_w[w].append(var_op(std_w[w][i],std[3],var))

    # 6) Display report
    print('Description of this report:\n')
    print('Description of the table')
    print('The first line shows the difference (relative or absolute, according'\
          'to the variable analyzed) between the QM method and the observed data'\
          'in the historical period. For precipitation, a value of 1 means that'\
          'the mean precipitation of the observed data and QM series are the same'\
          'for the historical period. For temperature, this is achieved with a'\
          'value of 0. The second line shows the headers of the periods reported.'\
          'Remember that the moving window is centered in the period and its'\
          'length is equal to the length of the historical period. The first '\
          'column (Hist.) shows the difference between the future and the'\
          'historical period. The next columns show the difference between the'\
          'corresponding period and the historical period. The third and the'\
          'following lines show the difference in each period for the'\
          'corresponding series. As a general rule, it is desirable that the'\
          'difference between the future periods and the historical period of the'\
          'bias corrected series match the differences of the modeled series.\n')

    print('Description of the figures')
    print('Figure 1: (left) shows the cumulative distribution functions of each'\
          'series split between historical and future period. (right) shows the'\
          'observed, modeled, and corrected series.')
    print('Figure 2: shows the monthly mean of the historical and future period.'\
          'Additionally, monthly mean of the projected periods is displayed.')
    print('Figure 3: shows the monthly standard deviation of the historical and'\
          'future period. Additionally, monthly mean of the projected periods is'\
          'displayed.')
    print('* In Fig. 2 and 3, the obj (red) line represents objective values.'\
          'This objective series is computed as the observed series (monthly mean'\
          'or standard deviation), scaled by the variation between the projected'\
          'and the historical period of the modeled series. Each projected period'\
          'is centered in  a moving window whose length is equal to the length of'\
          'the historical period.\n')
    df_m['Historical'] = dm[2:]
    df_s['Historical'] = ds[2:]
    for i in range(len(rep_head)):
        df_m[rep_head[i]] = dm_w[i]
        df_s[rep_head[i]] = ds_w[i]
    print('Mean')
    print(f'QM performance in historical period: {round(dm[0],3)}')
    print(df_m.round(3).loc[fun])
    print('')
    print('Standard deviation')
    print(f'QM performance in historical period: {round(ds[0],3)}')
    print(df_s.round(3).loc[fun])
    
    # =============================================================================
    # Figura 1: CDF and series
    # =============================================================================
    plt.figure(figsize=(13,6))
    plt.subplot(1,2,1)
    plt.title('Empirical distribution functions')
    if var == 0:
        plt.xlabel('Temperature (°C)')
    else:
        plt.xlabel('Precipitation (mm)')
    plt.ylabel('Probability')
    plt.grid()
    plt.plot(np.sort(obs[:,0]),ecdf(obs),'r',lw=1,label=r'Obs$_h$')
    plt.plot(np.sort(mod[:y_obs*12,0]),ecdf(mod[:y_obs*12,0]),'b',lw=1,label=r'Mod$_h$')
    plt.plot(np.sort(mod[y_obs*12:,0]),ecdf(mod[y_obs*12:,0]),'b--',lw=1,label=r'Mod$_f$')
    plt.plot(np.sort(QM_series[:y_obs*12]),ecdf(QM_series[:y_obs*12]),'k',lw=1,label=r'QM$_h$')
    if 'QM' in fun:
        plt.plot(np.sort(QM_series[y_obs*12:]),ecdf(QM_series[y_obs*12:]),'k--',lw=1,label=r'QM$_f$')
    if 'DQM' in fun:
        plt.plot(np.sort(DQM_series[y_obs*12:]),ecdf(DQM_series[y_obs*12:]),'g--',lw=1,label=r'DQM$_f$')
    if 'QDM' in fun:
        plt.plot(np.sort(QDM_series[y_obs*12:]),ecdf(QDM_series[y_obs*12:]),'y--',lw=1,label=r'QDM$_f$')
    if 'UQM' in fun:
        plt.plot(np.sort(UQM_series[y_obs*12:]),ecdf(UQM_series[y_obs*12:]),'c--',lw=1,label=r'UQM$_f$')
    if 'SDM' in fun:
        plt.plot(np.sort(SDM_series[y_obs*12:]),ecdf(SDM_series[y_obs*12:]),'m--',lw=1,label=r'SDM$_f$')
    plt.legend(loc='lower right')
    
    plt.subplot(1,2,2)
    plt.title('Time series')
    if var == 0:
        plt.ylabel('Temperature (°C)')
    else:
        plt.ylabel('Precipitation (mm)')
    plt.xlabel('Month since starting date')
    plt.grid()
    plt.plot(mod[:,0],'b',lw=1,label='Mod')
    if 'QM' in fun:
        plt.plot(QM_series,'k',lw=1,label='QM')
    if 'DQM' in fun:
        plt.plot(DQM_series,'g',lw=1,label='DQM')
    if 'QDM' in fun:
        plt.plot(QDM_series,'y',lw=1,label='QDM')
    if 'UQM' in fun:
        plt.plot(UQM_series,'c',lw=1,label='UQM')
    if 'SDM' in fun:
        plt.plot(SDM_series,'m',lw=1,label='SDM')
    plt.plot(obs[:,0],'r',lw=1,label='Obs')
    plt.legend(loc='upper center',ncol=7)
    plt.tight_layout()
    plt.show()
    
    # =============================================================================
    # Figura 2: Monthly mean
    # =============================================================================
    plt.figure(figsize=(13,6))
    plt.suptitle('Monthly mean',fontweight='bold')    
    
    plt.subplot(1,2+len(y_wind),1)
    plt.grid()
    plt.xlim(0,11)
    plt.title('Historical period')
    plt.xlabel('Month')
    if var == 0:
        plt.ylabel('Temperature (°C)')
    else:
        plt.ylabel('Precipitation (mm)')
    plt.plot(mu_M[1],'b',lw=1,label=r'Mod$_h$')
    plt.plot(mu_M[3],'k',lw=1,label=r'QM$_h$')
    plt.plot(mu_M[0],'r',lw=1,label=r'Obs')
    plt.legend(loc='best')
    plt.xticks(range(12),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],rotation=70)
    
    plt.subplot(1,2+len(y_wind),2)
    plt.grid()
    plt.xlim(0,11)
    plt.title('Future period')
    plt.xlabel('Month')
    if var == 0:
        plt.ylabel('Temperature (°C)')
    else:
        plt.ylabel('Precipitation (mm)')
    plt.plot(mu_M[2],'b',lw=1,label='Mod')
    if 'QM' in fun:
        plt.plot(mu_M[4],'k',lw=1,label='QM')
    if 'DQM' in fun:
        plt.plot(mu_M[5],'g',lw=1,label='DQM')
    if 'QDM' in fun:
        plt.plot(mu_M[6],'y',lw=1,label='QDM')
    if 'UQM' in fun:
        plt.plot(mu_M[7],'c',lw=1,label='UQM')
    if 'SDM' in fun:
        plt.plot(mu_M[8],'m',lw=1,label='SDM')
    if var == 0:
        plt.plot(mu_M[0]+SF_m[0],'r',lw=1,label=r'Obj')
    else:
        plt.plot(mu_M[0]*SF_m[0],'r',lw=1,label=r'Obj')
    plt.legend(loc='best')
    plt.xticks(range(12),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],rotation=70)
        
    for w in range(len(y_wind)):
        plt.subplot(1,2+len(y_wind),3+w)
        plt.grid()
        plt.xlim(0,11)
        plt.title(w_label[w])
        if var == 0:
            plt.ylabel('Temperature (°C)')
        else:
            plt.ylabel('Precipitation (mm)')
        plt.plot(mu_Mw[w][0],'b',lw=1,label='Mod')
        if 'QM' in fun:
            plt.plot(mu_Mw[w][1],'k',lw=1,label='QM')
        if 'DQM' in fun:
            plt.plot(mu_Mw[w][2],'g',lw=1,label='DQM')
        if 'QDM' in fun:
            plt.plot(mu_Mw[w][3],'y',lw=1,label='QDM')
        if 'UQM' in fun:
            plt.plot(mu_Mw[w][4],'c',lw=1,label='UQM')
        if 'SDM' in fun:
            plt.plot(mu_Mw[w][5],'m',lw=1,label='SDM')
        if var == 0:
            plt.plot(mu_M[0]+SF_m[w+1],'r',lw=1,label='Obj')
        else:
            plt.plot(mu_M[0]*SF_m[w+1],'r',lw=1,label='Obj')
        plt.legend(loc='best')
        plt.xlabel('Month')
        plt.xticks(range(12),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],rotation=70)
        
    plt.tight_layout()
    plt.show()
    
    # =============================================================================
    # Figura 3: Monthly std
    # =============================================================================
    plt.figure(figsize=(13,6))
    plt.suptitle('Monthly standard deviation',fontweight='bold')    
    plt.subplot(1,2+len(y_wind),1)
    plt.grid()
    plt.xlim(0,11)
    plt.title('Historical period')
    plt.xlabel('Month')
    if var == 0:
        plt.ylabel('Temperature (°C)')
    else:
        plt.ylabel('Precipitation (mm)')
    plt.plot(std_M[1],'b',lw=1,label=r'Mod$_h$')
    plt.plot(std_M[3],'k',lw=1,label=r'QM$_h$')
    plt.plot(std_M[0],'r',lw=1,label=r'Obs')
    plt.legend(loc='best')
    plt.xticks(range(12),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],rotation=70)
    
    plt.subplot(1,2+len(y_wind),2)
    plt.grid()
    plt.xlim(0,11)
    plt.title('Future period')
    plt.xlabel('Month')
    if var == 0:
        plt.ylabel('Temperature (°C)')
    else:
        plt.ylabel('Precipitation (mm)')
    plt.plot(std_M[2],'b',lw=1,label='Mod')
    if 'QM' in fun:
        plt.plot(std_M[4],'k',lw=1,label='QM')
    if 'DQM' in fun:
        plt.plot(std_M[5],'g',lw=1,label='DQM')
    if 'QDM' in fun:
        plt.plot(std_M[6],'y',lw=1,label='QDM')
    if 'UQM' in fun:
        plt.plot(std_M[7],'c',lw=1,label='UQM')
    if 'SDM' in fun:
        plt.plot(std_M[8],'m',lw=1,label='SDM')
    if var == 0:
        plt.plot(std_M[0]+SF_s[0],'r',lw=1,label=r'Obj')
    else:
        plt.plot(std_M[0]*SF_s[0],'r',lw=1,label=r'Obj')
    plt.legend(loc='best')
    plt.xticks(range(12),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],rotation=70)
    
    for w in range(len(y_wind)):
        plt.subplot(1,2+len(y_wind),3+w)
        plt.grid()
        plt.xlim(0,11)
        plt.title(w_label[w])
        if var == 0:
            plt.ylabel('Temperature (°C)')
        else:
            plt.ylabel('Precipitation (mm)')
        plt.plot(std_Mw[w][0],'b',lw=1,label='Mod')
        if 'QM' in fun:
            plt.plot(std_Mw[w][1],'k',lw=1,label='QM')
        if 'DQM' in fun:
            plt.plot(std_Mw[w][2],'g',lw=1,label='DQM')
        if 'QDM' in fun:
            plt.plot(std_Mw[w][3],'y',lw=1,label='QDM')
        if 'UQM' in fun:
            plt.plot(std_Mw[w][4],'c',lw=1,label='UQM')
        if 'SDM' in fun:
            plt.plot(std_Mw[w][5],'m',lw=1,label='SDM')
        if var == 0:
            plt.plot(std_M[0]+SF_s[w+1],'r',lw=1,label=r'Obj')
        else:
            plt.plot(std_M[0]*SF_s[w+1],'r',lw=1,label=r'Obj')
        plt.legend(loc='best')
        plt.xlabel('Month')
        plt.xticks(range(12),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],rotation=70)
        
    plt.tight_layout()
    plt.show()
    
    return QM_series,DQM_series,QDM_series,UQM_series,SDM_series