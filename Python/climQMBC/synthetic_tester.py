from .methods import QM, DQM, QDM, UQM, SDM
from scipy.stats import skew as sk
import numpy as np
import pandas as pd

"""
This script contains the synthetic tester function to compare the performance 
of different methods based on a synthetic example rather than the original
data, to compare the theoretical performance.

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
Revision: 0, updated Aug 2026
"""
    

def synthetic_tester(obs, mod, allow_negatives, mult_change, SDM_var, yn=10000):
    # Internal functions to generate the synthetic gamma and normal values
    def get_synthetic_gamma(mu,sigma,yn):
        k = (mu/sigma)**2
        th = (sigma**2)/mu
        return np.random.gamma(k,th,yn)


    def get_synthetic_normal(mu,sigma,yn):
        return np.random.normal(mu,sigma,yn)

    # Add bias considering an equiprobability (normal) transform
    def add_bias(series, mu_bias, sigma_bias, mult_change):
        mu = series.mean()
        sigma = series.std()
        
        if mult_change == 1:
            mu_scaled = mu*(1 + mu_bias)
            sigma_scaled = sigma*(1 + sigma_bias)
        else:
            mu_scaled = mu + mu_bias
            sigma_scaled = sigma + sigma_bias
            
        biased_series = mu_scaled + (series - mu)*sigma_scaled/sigma
        return biased_series

    # Compute NSE, KGE and ratios in statistical moments and percentiles
    def get_tester_stats(bc_series,obj_series,mult_change):
        bc_mean = np.mean(bc_series)
        obj_mean = np.mean(obj_series)
        
        a =  np.sum((bc_series-bc_mean)*(obj_series-obj_mean))
        b = np.sqrt(np.sum((bc_series-bc_mean)**2))
        c = np.sqrt(np.sum((obj_series-obj_mean)**2))
        
        r2 = (a/(b*c))**2
        
        if mult_change == 1:
            mu_ratio = bc_series.mean()/obj_series.mean()
            sigma_ratio = bc_series.std()/obj_series.std()
            sk_ratio = sk(bc_series)/sk(obj_series)
            p05_ratio = np.percentile(bc_series,5)/np.percentile(obj_series,5)
            p10_ratio = np.percentile(bc_series,10)/np.percentile(obj_series,10)
            p25_ratio = np.percentile(bc_series,25)/np.percentile(obj_series,25)
            p50_ratio = np.percentile(bc_series,50)/np.percentile(obj_series,50)
            p75_ratio = np.percentile(bc_series,75)/np.percentile(obj_series,75)
            p90_ratio = np.percentile(bc_series,90)/np.percentile(obj_series,90)
            p95_ratio = np.percentile(bc_series,95)/np.percentile(obj_series,95)
            
            kge = 1-np.sqrt((1-r2)**2 + (1-mu_ratio)**2 + (1-sigma_ratio)**2)
            nse = 1-np.sum((bc_series-obj_series)**2)/np.sum((obj_series.mean()-obj_series)**2)
            
        else:
            mu_ratio = bc_series.mean()-obj_series.mean()
            sigma_ratio = bc_series.std()-obj_series.std()
            sk_ratio = sk(bc_series)-sk(obj_series)
            p05_ratio = np.percentile(bc_series,5)-np.percentile(obj_series,5)
            p10_ratio = np.percentile(bc_series,10)-np.percentile(obj_series,10)
            p25_ratio = np.percentile(bc_series,25)-np.percentile(obj_series,25)
            p50_ratio = np.percentile(bc_series,50)-np.percentile(obj_series,50)
            p75_ratio = np.percentile(bc_series,75)-np.percentile(obj_series,75)
            p90_ratio = np.percentile(bc_series,90)-np.percentile(obj_series,90)
            p95_ratio = np.percentile(bc_series,95)-np.percentile(obj_series,95)      
            
            kge = np.sqrt((1-r2)**2 + (mu_ratio)**2 + (sigma_ratio)**2)
            nse = np.sum((bc_series-obj_series)**2)/np.sum((obj_series.mean()-obj_series)**2)

        li_stats = [nse,kge,mu_ratio,sigma_ratio,sk_ratio,
                    p05_ratio,p10_ratio,p25_ratio,p50_ratio,
                    p75_ratio,p90_ratio,p95_ratio]
        
        return np.round(li_stats,3)   
    
    # Get series statistics
    mu_obs_h = obs.mean()
    sigma_obs_h = obs.std()

    mu_mod_h = mod[:obs.shape[0]].mean()
    sigma_mod_h = mod[:obs.shape[0]].std()

    mu_mod_f = mod[obs.shape[0]:].mean()
    sigma_mod_f = mod[obs.shape[0]:].std()

    # Get bias and future values
    if mult_change == 1:
        mu_bias = mu_mod_h/mu_obs_h - 1
        sigma_bias = sigma_mod_h/sigma_obs_h - 1
        
        mu_change = mu_mod_f/mu_mod_h - 1
        sigma_change = sigma_mod_f/sigma_mod_h - 1
        
        mu_obj_f = mu_obs_h*(1+mu_change)
        sigma_obj_f = sigma_obs_h*(1+sigma_change)

    else:
        mu_bias = mu_mod_h-mu_obs_h
        sigma_bias = sigma_mod_h-sigma_obs_h
        
        mu_change = mu_mod_f-mu_mod_h
        sigma_change = sigma_mod_f-sigma_mod_h
        
        mu_obj_f = mu_obs_h + mu_change
        sigma_obj_f = sigma_obs_h + sigma_change
        
    # Get synthetic series
    if allow_negatives == 0:
        obs_h = get_synthetic_gamma(mu_obs_h,sigma_obs_h,yn)
        obs_f = get_synthetic_gamma(mu_obj_f,sigma_obj_f,yn)
        
    else:
        obs_h = get_synthetic_normal(mu_obs_h,sigma_obs_h,yn)
        obs_f = get_synthetic_normal(mu_obj_f,sigma_obj_f,yn)
        
    mod_h = add_bias(obs_h,mu_bias,sigma_bias,mult_change)  
    mod_f = add_bias(obs_f,mu_bias,sigma_bias,mult_change)


    # Concat historical and future period
    ## Future period is concated twice so that only the second future is
    ## analyzed, without moving windows
    obs = np.hstack([obs_h])
    mod = np.hstack([mod_h,mod_f,mod_f])


    # Apply bias correction methods
    QM_h_series = QM(obs_h,mod_h,allow_negatives=allow_negatives,frq='A')
    SDM_h_series = SDM(obs_h,mod_h,SDM_var=SDM_var,frq='A')
    QM_f_series = QM(obs,mod,allow_negatives=allow_negatives,frq='A')[2*yn:]
    DQM_f_series = DQM(obs,mod,allow_negatives=allow_negatives,mult_change=mult_change,frq='A')[2*yn:]
    QDM_f_series = QDM(obs,mod,allow_negatives=allow_negatives,mult_change=mult_change,frq='A')[2*yn:]
    UQM_f_series = UQM(obs,mod,allow_negatives=allow_negatives,mult_change=mult_change,frq='A')[2*yn:]
    SDM_f_series = SDM(obs,mod,SDM_var=SDM_var,frq='A')[2*yn:]

    li_stats_qm_h = get_tester_stats(QM_h_series,obs_h,mult_change)
    li_stats_sdm_h = get_tester_stats(SDM_h_series,obs_h,mult_change)
    li_stats_qm_f = get_tester_stats(QM_f_series,obs_f,mult_change)
    li_stats_dqm_f = get_tester_stats(DQM_f_series,obs_f,mult_change)
    li_stats_qdm_f = get_tester_stats(QDM_f_series,obs_f,mult_change)
    li_stats_uqm_f = get_tester_stats(UQM_f_series,obs_f,mult_change)
    li_stats_sdm_f = get_tester_stats(SDM_f_series,obs_f,mult_change)

    df_tester_stats = pd.DataFrame([li_stats_qm_h,li_stats_sdm_h,
                                    li_stats_qm_f,li_stats_dqm_f,li_stats_qdm_f,li_stats_uqm_f,li_stats_sdm_f],
                                   index=['QM_h','SDM_h','QM_f','DQM_f','QDM_f','UQM_f','SDM_f'],
                                   columns=['NSE','KGE','Mean','Std. Dev','Skew','P05','P10','P25','P50','P75','P90','P95'])
    
    return df_tester_stats