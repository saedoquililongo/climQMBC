#' Quantile Mapping
#'
#' This function performs bias correction of modeled series based on observed data by the Quantile Mapping (QM) method, as described by Cannon et al. (2015). Correction is performed to daily, monthly or annual data in a single location. An independent probability distribution function is assigned to  each sub-annual period based on the Kolmogorov-Smirnov test.
#'
#' The available distributions are: 1) Normal ; 2) Log-Normal ; 3) Gamma 2 parameters ; 4) Gamma 3 parameters ; (Pearson 3 parameters) ; 5) Log-Gamma 3 parameters ; 6) Gumbel ; 7) Exponential
#'
#' NOTE: This routine considers that obs and mod series start in the same day/month/year and are continuous until the end day/month/year.
#'
#' Cannon, A. J., S. R. Sobie, and T. Q. Murdock, 2015: Bias correction of GCM precipitation by quantile mapping: How well do methods preserve changes in quantiles and extremes? J. Climate, 28(17), 6938-6959, https://doi.org/10.1175/JCLI-D-14-00754.1
#'
#' @param obs A column vector of daily, monthly or annual observed data, without considering leap days. If daily frequency is specified, the length of the column vector should by a multiple of 365 and for monthly frequency, it should be a multiple of 12. [ndata_obs, 1]
#' @param mod A column vector of daily, monthly or annual modeled or GCM data, without considering leap days. If daily frequency is specified, the length of the column vector should by a multiple of 365 and for monthly frequency, it should be a multiple of 12. [ndata_mod, 1]
#' @param allow_negatives (Optional) A flag that identifies if data allows negative values and also to replace no-rain values with random small values (Chadwick et al., 2023) to avoid numerical problems with the probability distribution functions. allow_negatives = 1 or True: Allow negatives (default) ; allow_negatives = 0 or False: Do not allow negative
#' @param frq (Optional) A string specifying if the input frequency is daily, monthly or annual. Daily:     frq = 'D' ; Monthly:   frq = 'M' ; Annual:    frq = 'A' (default)
#' @param pp_threshold (Optional) A float indicating the threshold to consider no-rain values. Default is 1.
#' @param pp_factor (Optional) A float which multiplied to pp_threshold indicates the maximum value of the random values that replace no-rain values. Default is 1/100.
#' @param day_win (Optional) (only for frq='D') An integer indicating how many days to consider backwards and forward to get the statistics of each calendar day.The length of the window will be  (2*win_day-1). For example, day_win=15 -> window of 29. Default: win = 1
#' @param user_pdf (Optional) A flag indicating if the user will define the probability distribution functions (pdf) for the observed and modeled series. The distributions will be the same for all periods and sub-periods. user_pdf = 1 or True: User defines the pdf ; user_pdf = 0 or False: pdf defined by the Kolmogorov -Smirnov test (default)
#' @param pdf_obs (Optional) An integer indicating the probability distribution function (pdf) to be used for the observed data. The pdf will be the same for all periods and sub-periods. Default: None
#' @param pdf_mod (Optional) An integer indicating the probability distribution function (pdf) to be used for the modeled data. The pdf will be the same for all periods and sub-periods. Default: None
#'
#' @return QM_series: A column vector of data bias corrected with the QM method. [ndata_mod, 1]
#' @export
#'
#' @examples QM(obs,mod)
#' @examples QM(obs,mod,frq='A')
#' @examples QM(obs,mod,frq='M')
#' @examples QM(obs,mod,frq='D')
#' @examples QM(obs,mod,allow_negatives=1)
#' @examples QM(obs,mod,allow_negatives=0)
#' @examples QM(obs,mod,allow_negatives=0,frq='D',pp_factor=1/10000,day_win=15)
QM <- function(obs, mod, allow_negatives, frq, pp_threshold, pp_factor, day_win,
               user_pdf, pdf_obs, pdf_mod){
  
  # Check for optional arguments
  if(missing(allow_negatives)) {
    allow_negatives <- 1
  }
  
  if(missing(frq)) {
    frq <- 'A'
  }
  
  if(missing(pp_threshold)) {
    pp_threshold <- 1
  }
  
  if(missing(pp_factor)) {
    pp_factor <- 1/100
  }
  
  if(missing(day_win)) {
    day_win <- 1
  }
  
  if(missing(user_pdf)) {
    user_pdf <- FALSE
  }
  
  if(missing(pdf_obs)) {
    pdf_obs <- FALSE
  }
  
  if(missing(pdf_mod)) {
    pdf_mod <- FALSE
  }
  
  # 0) Only for daily frequency data that does not allow negatives (as
  #    precipitacion). Get a pp_threshold for the modeled series that result
  #    in the same amount of rain-day values in the historical period as in
  #    the observed series
  if (frq=='D' & allow_negatives==FALSE) {
    pp_threshold_mod <- get_pp_threshold_mod(obs, mod, pp_threshold)
  } else {
    pp_threshold_mod <- pp_threshold
  }
  
  
  # 1) Format series to matrix with rows as sub-periods and columns as years
  #    and, if needed, replace no-rain values with random small values 
  format_list_obs <- formatQM(obs, allow_negatives, frq, pp_threshold, pp_factor)
  y_obs <- format_list_obs[[1]]
  obs_series <- format_list_obs[[2]]
  
  format_list_mod <- formatQM(mod, allow_negatives, frq, pp_threshold_mod, pp_factor)
  y_mod <- format_list_mod[[1]]
  mod_series <- format_list_mod[[2]]
  
  # 2) Get statistics of the observed and modeled series in the historical
  #    period.
  if (frq=='D') {
    # Get a 3D array with a moving window centered on each day
    obs_series_moving <- day_centered_moving_window(obs_series, day_win)
    mod_series_moving <- day_centered_moving_window(mod_series, day_win)

    # Reshape to 2D in order to have all the days of the window in each row
    obs_series_moving <- array(aperm(obs_series_moving, c(1,3,2)), c(365,(2*day_win-1)*y_obs))
    modh_series_moving <- array(aperm(mod_series_moving,c(1,3,2))[,1:y_obs,], c(365,(2*day_win-1)*y_obs))
    
    # Replace no-rain values with nans but keep at least 30 values to get
    # the statistics. If fewer than 30 values are rain values, the missing
    # ones will be previously random values replaced
    if (allow_negatives==FALSE) {
      obs_series_moving <- set_norain_to_nan(obs_series_moving, pp_threshold, pp_factor)
      modh_series_moving <- set_norain_to_nan(modh_series_moving, pp_threshold_mod, pp_factor)
      mod_series[mod_series<pp_threshold_mod] <- NaN
    }
    
    # Get the statistics for each row
    stats_list_obs <- getStats(obs_series_moving, frq)
    mu_obs <- stats_list_obs[[1]]
    std_obs <- stats_list_obs[[2]]
    skew_obs <- stats_list_obs[[3]]
    skewy_obs <- stats_list_obs[[4]]
    
    stats_list_mod <- getStats(modh_series_moving, frq)
    mu_mod <- stats_list_mod[[1]]
    std_mod <- stats_list_mod[[2]]
    skew_mod <- stats_list_mod[[3]]
    skewy_mod <- stats_list_mod[[4]]
    
  } else {
    # Get the statistics for each row
    stats_list_obs <- getStats(obs_series, frq)
    mu_obs <- stats_list_obs[[1]]
    std_obs <- stats_list_obs[[2]]
    skew_obs <- stats_list_obs[[3]]
    skewy_obs <- stats_list_obs[[4]]
    
    stats_list_mod <- getStats(matrix(mod_series[,1:y_obs],nrow=dim(mod_series)[1]), frq)
    mu_mod <- stats_list_mod[[1]]
    std_mod <- stats_list_mod[[2]]
    skew_mod <- stats_list_mod[[3]]
    skewy_mod <- stats_list_mod[[4]]
  }

  # 3) Assign a probability distribution function to each sub-period for the
  #    observed and modeled data in the historical period.
  if (user_pdf==FALSE){
    if (frq=='D') {
      list_getdist <- getDist(obs_series_moving, allow_negatives, mu_obs,std_obs,skew_obs,skewy_obs)
      pdf_obs <- list_getdist[[1]] 
      ks_fail_obs <- list_getdist[[2]]
      
      list_getdist <- getDist(modh_series_moving, allow_negatives, mu_mod,std_mod,skew_mod,skewy_mod)
      pdf_mod <- list_getdist[[1]] 
      ks_fail_mod <- list_getdist[[2]] 
    } else {
      list_getdist <- getDist(obs_series, allow_negatives, mu_obs,std_obs,skew_obs,skewy_obs)
      pdf_obs <- list_getdist[[1]] 
      ks_fail_obs <- list_getdist[[2]]
      
      list_getdist <- getDist(matrix(mod_series[,1:y_obs],nrow=dim(mod_series)[1]),allow_negatives, mu_mod,std_mod,skew_mod,skewy_mod)
      pdf_mod <- list_getdist[[1]] 
      ks_fail_mod <- list_getdist[[2]] 
    }
    
  } else {
    pdf_obs <- array(0, dim(obs_series)[1]) + pdf_obs
    pdf_mod <- array(0, dim(mod_series)[1]) + pdf_mod
  }
  
  # 4) Apply the cumulative distribution function of the modeled data,
  #    evaluated with the statistics of the modeled data in the historical
  #    period and the modeled data. Equation 1 of Cannon et al. (2015).
  prob <- getCDF(pdf_mod,mod_series,mu_mod,std_mod,skew_mod,skewy_mod)
  prob[is.nan(mod_series)] <- NaN

  # 5) Apply the inverse cumulative distribution function of the observed
  #    data, evaluated with the statistics of the observed data in the
  #    historical period, to the probabilities obtained from 2). Equation 1 of
  #    Cannon et al. (2015).
  QM_series <- getCDFinv(pdf_obs,prob,mu_obs,std_obs,skew_obs,skewy_obs)
  
  # 6) Reshape to a column vector
  QM_series <- matrix(QM_series)
  
  # 7) If does not allow negative values, replace nans (should be the result 
  #    of correcting removed no-rain values) and replace no-rain values with 0
  if (allow_negatives==0){
    QM_series[is.nan(QM_series)] <- 0
    QM_series[QM_series<pp_threshold] <- 0
  }
  
  if (user_pdf==FALSE){
    ks_fail <- ks_fail_obs + ks_fail_mod
    if (ks_fail>0) {
      print('QM: Some of the probability distribution functions did not pass the KS-Test')
    }
  }

  return(QM_series)
}


#' Detrended Quantile Mapping
#'
#' This function performs bias correction of modeled series based on observed data by the Detrended Quantile Mapping (DQM) method, as described by Cannon et al. (2015). Correction is performed to daily, monthly or annual data in a single location. An independent probability distribution function is assigned to each sub-annual period based on the Kolmogorov-Smirnov test. Each projected period considers a window whose length is equal to the number of years of the historical period and ends in the analyzed year.
#' 
#' Correction of the historical period is performed by the Quantile Mapping (QM) method, as described in Cannon et al. (2015) and Chadwick et al. (2023).
#'
#' The available distributions are: 1) Normal ; 2) Log-Normal ; 3) Gamma 2 parameters ; 4) Gamma 3 parameters ; (Pearson 3 parameters) ; 5) Log-Gamma 3 parameters ; 6) Gumbel ; 7) Exponential
#'
#' NOTE: This routine considers that obs and mod series start in the same day/month/year and are continuous until the end day/month/year.
#' 
#' Cannon, A. J., S. R. Sobie, and T. Q. Murdock, 2015: Bias correction of GCM precipitation by quantile mapping: How well do methods preserve changes in quantiles and extremes? J. Climate, 28(17), 6938-6959, https://doi.org/10.1175/JCLI-D-14-00754.1
#'
#' @param obs A column vector of daily, monthly or annual observed data, without considering leap days. If daily frequency is specified, the length of the column vector should by a multiple of 365 and for monthly frequency, it should be a multiple of 12. [ndata_obs, 1]
#' @param mod A column vector of daily, monthly or annual modeled or GCM data, without considering leap days. If daily frequency is specified, the length of the column vector should by a multiple of 365 and for monthly frequency, it should be a multiple of 12. [ndata_mod, 1]
#' @param mult_change (Optional) A flag that indicates if projected changes should be computed as multiplicative (fut = hist*delta) or  additive (fut = hist + delta) changes. mult_change = 1 or True: Multiplicative (default) ; mult_change = 0 or False: Additive
#' @param allow_negatives (Optional) A flag that identifies if data allows negative values and also to replace no-rain values with random small values (Chadwick et al., 2023) to avoid numerical problems with the probability distribution functions. allow_negatives = 1 or True: Allow negatives (default) ; allow_negatives = 0 or False: Do not allow negative
#' @param frq (Optional) A string specifying if the input frequency is daily, monthly or annual. Daily:     frq = 'D' ; Monthly:   frq = 'M' ; Annual:    frq = 'A' (default)
#' @param pp_threshold (Optional) A float indicating the threshold to consider no-rain values. Default is 1.
#' @param pp_factor (Optional) A float which multiplied to pp_threshold indicates the maximum value of the random values that replace no-rain values. Default is 1/100.
#' @param day_win (Optional) (only for frq='D') An integer indicating how many days to consider backwards and forward to get the statistics of each calendar day.The length of the window will be  (2*win_day-1). For example, day_win=15 -> window of 29. Default: win = 1
#' @param user_pdf (Optional) A flag indicating if the user will define the probability distribution functions (pdf) for the observed and modeled series. The distributions will be the same for all periods and sub-periods. user_pdf = 1 or True: User defines the pdf ; user_pdf = 0 or False: pdf defined by the Kolmogorov -Smirnov test (default)
#' @param pdf_obs (Optional) An integer indicating the probability distribution function (pdf) to be used for the observed data. The pdf will be the same for all periods and sub-periods. Default: None
#' @param pdf_mod (Optional) An integer indicating the probability distribution function (pdf) to be used for the modeled data. The pdf will be the same for all periods and sub-periods. Default: None
#'
#' @return DQM_series: A column vector of data bias corrected with the DQM method. [ndata_mod, 1]
#' @export
#'
#' @examples DQM(obs,mod)
#' @examples DQM(obs,mod,frq='A')
#' @examples DQM(obs,mod,frq='M')
#' @examples DQM(obs,mod,frq='D')
#' @examples DQM(obs,mod,mult_changes=0,allow_negatives=1)
#' @examples DQM(obs,mod,mult_changes=1,allow_negatives=0)
#' @examples DQM(obs,mod,mult_changes=1,allow_negatives=0,frq='D',pp_factor=1/10000,day_win=15)
DQM <- function(obs, mod, mult_change, allow_negatives, frq, pp_threshold,
                pp_factor, day_win, user_pdf, pdf_obs, pdf_mod){
  
  # Check for optional arguments
  if(missing(mult_change)) {
    mult_change <- 1
  }
  
  if(missing(allow_negatives)) {
    allow_negatives <- 1
  }

  if(missing(frq)) {
    frq <- 'A'
  }
  
  if(missing(pp_threshold)) {
    pp_threshold <- 1
  }

  if(missing(pp_factor)) {
    pp_factor <- 1/100
  }
  
  if(missing(day_win)) {
    day_win <- 1
  }
  
  if(missing(user_pdf)) {
    user_pdf <- FALSE
  }
  
  if(missing(pdf_obs)) {
    pdf_obs <- FALSE
  }
  
  if(missing(pdf_mod)) {
    pdf_mod <- FALSE
  }
  
  # 0) Only for daily frequency data that does not allow negatives (as
  #    precipitacion). Get a pp_threshold for the modeled series that result
  #    in the same amount of rain-day values in the historical period as in
  #    the observed series
  if (frq=='D' & allow_negatives==FALSE) {
    pp_threshold_mod <- get_pp_threshold_mod(obs, mod, pp_threshold)
  } else {
    pp_threshold_mod <- pp_threshold
  }
  
  # 1) Format series to matrix with rows as sub-periods and columns as years
  #    and, if needed, replace no-rain values with random small values 
  format_list_obs <- formatQM(obs, allow_negatives, frq, pp_threshold, pp_factor)
  y_obs <- format_list_obs[[1]]
  obs_series <- format_list_obs[[2]]
  
  format_list_mod <- formatQM(mod, allow_negatives, frq, pp_threshold_mod, pp_factor)
  y_mod <- format_list_mod[[1]]
  mod_series <- format_list_mod[[2]]
  
  # 2) Get statistics of the observed and modeled series in the historical
  #    period.
  if (frq=='D') {
    # Get a 3D array with a moving window centered on each day
    obs_series_moving <- day_centered_moving_window(obs_series, day_win)
    mod_series_moving <- day_centered_moving_window(mod_series, day_win)
    
    # Get a 4D array with a backward moving window for each projected period
    win_series <- projected_backward_moving_window(mod_series_moving, y_obs, frq)
    
    # Reshape to 2D in order to have all the days of the window in each row
    obs_series_moving <- array(aperm(obs_series_moving, c(1,3,2)), c(365,(2*day_win-1)*y_obs))
    modh_series_moving <- array(aperm(mod_series_moving,c(1,3,2))[,1:y_obs,], c(365,(2*day_win-1)*y_obs))
    
    # Reshape to 3D in order to have all the days of the window in each row
    win_series_moving <- array(aperm(win_series,c(1,4,2,3)), c(365, (2*day_win-1)*y_obs, y_mod-y_obs))

    # Replace no-rain values with nans but keep at least 30 values to get
    # the statistics. If fewer than 30 values are rain values, the missing
    # ones will be previously random values replaced
    if (allow_negatives==FALSE) {
      obs_series_moving <- set_norain_to_nan(obs_series_moving, pp_threshold, pp_factor)
      modh_series_moving <- set_norain_to_nan(modh_series_moving, pp_threshold_mod, pp_factor)
      mod_series[mod_series<pp_threshold_mod] <- NaN
      win_series_moving <- set_norain_to_nan(win_series_moving, pp_threshold_mod, pp_factor)
    }
    
    # Get the statistics for each row
    stats_list_obs <- getStats(obs_series_moving, frq)
    mu_obs <- stats_list_obs[[1]]
    std_obs <- stats_list_obs[[2]]
    skew_obs <- stats_list_obs[[3]]
    skewy_obs <- stats_list_obs[[4]]
    
    stats_list_mod <- getStats(modh_series_moving, frq)
    mu_mod <- stats_list_mod[[1]]
    std_mod <- stats_list_mod[[2]]
    skew_mod <- stats_list_mod[[3]]
    skewy_mod <- stats_list_mod[[4]]
    
    mu_win <- apply(win_series_moving,c(1,3),mean, na.rm=TRUE)
    
  } else {
    # Get a 3D array with a backward moving window for each projected period
    win_series <- projected_backward_moving_window(mod_series, y_obs, frq)
    
    # Get the statistics for each row
    stats_list_obs <- getStats(obs_series, frq)
    mu_obs <- stats_list_obs[[1]]
    std_obs <- stats_list_obs[[2]]
    skew_obs <- stats_list_obs[[3]]
    skewy_obs <- stats_list_obs[[4]]
    
    stats_list_mod <- getStats(matrix(mod_series[,1:y_obs],nrow=dim(mod_series)[1]), frq)
    mu_mod <- stats_list_mod[[1]]
    std_mod <- stats_list_mod[[2]]
    skew_mod <- stats_list_mod[[3]]
    skewy_mod <- stats_list_mod[[4]]
    
    mu_win <- apply(win_series,c(1,2),mean, na.rm=TRUE)
  }
  
  # 3) Assign a probability distribution function to each sub-period for the
  #    observed and modeled data in the historical period.
  if (user_pdf==FALSE){
    if (frq=='D') {
      list_getdist <- getDist(obs_series_moving, allow_negatives, mu_obs,std_obs,skew_obs,skewy_obs)
      pdf_obs <- list_getdist[[1]] 
      ks_fail_obs <- list_getdist[[2]]
      
      list_getdist <- getDist(modh_series_moving, allow_negatives, mu_mod,std_mod,skew_mod,skewy_mod)
      pdf_mod <- list_getdist[[1]] 
      ks_fail_mod <- list_getdist[[2]] 
    } else {
      list_getdist <- getDist(obs_series, allow_negatives, mu_obs,std_obs,skew_obs,skewy_obs)
      pdf_obs <- list_getdist[[1]] 
      ks_fail_obs <- list_getdist[[2]]
      
      list_getdist <- getDist(matrix(mod_series[,1:y_obs],nrow=dim(mod_series)[1]),allow_negatives, mu_mod,std_mod,skew_mod,skewy_mod)
      pdf_mod <- list_getdist[[1]] 
      ks_fail_mod <- list_getdist[[2]] 
    }
    
  } else {
    pdf_obs <- array(0, dim(obs_series)[1]) + pdf_obs
    pdf_mod <- array(0, dim(mod_series)[1]) + pdf_mod
  }
  
  # 4) Compute the linear scaled values (value in square brackets in equation
  #    2 of Cannon et al. (2015)).
  mu_mod_repeated <- array(rep(mu_mod,(y_mod-y_obs)),c(dim(mod_series)[1],(y_mod-y_obs)))
  if (mult_change == 1){ # Precipitation
    scaling_factor <- mu_win/mu_mod_repeated
    detrended_series <- mod_series[,(y_obs+1):dim(mod_series)[2]]/scaling_factor
  } else{        # Temperature
    scaling_factor <- mu_win - mu_mod_repeated
    detrended_series <- mod_series[,(y_obs+1):dim(mod_series)[2]] - scaling_factor
  }
  

  # 5) Apply the cumulative distribution function of the modeled data,
  #    evaluated with the statistics of the modeled data in the historical
  #    period and the modeled data. Equation 2 of Cannon et al. (2015).
  prob <- getCDF(pdf_mod,detrended_series,mu_mod,std_mod,skew_mod,skewy_mod)
  prob[is.nan(detrended_series)] <- NaN

  # 6) Apply the inverse cumulative distribution function of the observed
  #    data, evaluated with the statistics of the observed data in the
  #    historical period, to the probabilities obtained from 5). Equation 2 of
  #    Cannon et al. (2015).
  DQM_raw <- getCDFinv(pdf_obs,prob,mu_obs,std_obs,skew_obs,skewy_obs)

  # 7) Reimpose the trend to the values obtained in 6). Equation 2 of Cannon
  #    et al. (2015).
  if (mult_change == 1){ # Precipitation
    DQM <- DQM_raw*scaling_factor
  } else{        # Temperature
    DQM <- DQM_raw + scaling_factor
  }
  
  # 8) Perform QM for the historical period.
  mod_h <- mod[1:length(obs)]
  mod_h <- matrix(mod_h)
  QM_series <- QM(obs, mod_h, allow_negatives,frq,pp_threshold, pp_factor, day_win, user_pdf, pdf_obs, pdf_mod)
  
  # 9) Reshape to a column vector
  DQM <- matrix(DQM)

  # 10) Concat QM (hsitorical) and DQM (future)
  DQM_series <- matrix(c(QM_series,DQM))
  
  # 11) If does not allow negative values, replace nans (should be the result
  #     of correcting removed no-rain values) and replace no-rain values with 0
  if (allow_negatives == 0){
    DQM_series[is.nan(DQM_series)] <- 0
    DQM_series[DQM_series<pp_threshold] <- 0
  }
  
  if (user_pdf==FALSE){
    ks_fail <- ks_fail_obs + ks_fail_mod
    if (ks_fail>0) {
      print('DQM: Some of the probability distribution functions did not pass the KS-Test')
    }
  }
  

  return(DQM_series)
}


#' Quantile Delta Mapping
#'
#' This function performs bias correction of modeled series based on observed data by the Quantile Delta Mapping (QDM) method, as described by Cannon et al. (2015). Correction is performed to daily, monthly or annual data in a single location. An independent probability distribution function is  assigned to each sub-annual period based on the Kolmogorov-Smirnov test. Each projected period considers a window whose length is equal to the number of years of the historical period and ends in the analyzed year.
#' 
#' Correction of the historical period is performed by the Quantile Mapping (QM) method, as described in Cannon et al. (2015) and Chadwick et al. (2023).
#'
#' The available distributions are: 1) Normal ; 2) Log-Normal ; 3) Gamma 2 parameters ; 4) Gamma 3 parameters ; (Pearson 3 parameters) ; 5) Log-Gamma 3 parameters ; 6) Gumbel ; 7) Exponential
#'
#' NOTE: This routine considers that obs and mod series start in the same day/month/year and are continuous until the end day/month/year.
#' 
#' Cannon, A. J., S. R. Sobie, and T. Q. Murdock, 2015: Bias correction of GCM precipitation by quantile mapping: How well do methods preserve changes in quantiles and extremes? J. Climate, 28(17), 6938-6959, https://doi.org/10.1175/JCLI-D-14-00754.1
#'
#' @param obs A column vector of daily, monthly or annual observed data, without considering leap days. If daily frequency is specified, the length of the column vector should by a multiple of 365 and for monthly frequency, it should be a multiple of 12. [ndata_obs, 1]
#' @param mod A column vector of daily, monthly or annual modeled or GCM data, without considering leap days. If daily frequency is specified, the length of the column vector should by a multiple of 365 and for monthly frequency, it should be a multiple of 12. [ndata_mod, 1]
#' @param mult_change (Optional) A flag that indicates if projected changes should be computed as multiplicative (fut = hist*delta) or  additive (fut = hist + delta) changes. mult_change = 1 or True: Multiplicative (default) ; mult_change = 0 or False: Additive
#' @param allow_negatives (Optional) A flag that identifies if data allows negative values and also to replace no-rain values with random small values (Chadwick et al., 2023) to avoid numerical problems with the probability distribution functions. allow_negatives = 1 or True: Allow negatives (default) ; allow_negatives = 0 or False: Do not allow negative
#' @param frq (Optional) A string specifying if the input frequency is daily, monthly or annual. Daily:     frq = 'D' ; Monthly:   frq = 'M' ; Annual:    frq = 'A' (default)
#' @param pp_threshold (Optional) A float indicating the threshold to consider no-rain values. Default is 1.
#' @param pp_factor (Optional) A float which multiplied to pp_threshold indicates the maximum value of the random values that replace no-rain values. Default is 1/100.
#' @param rel_change_th (Optional)  A float indicating the maximum scaling factor  (Equation 4 of Cannon et al. (2015)) when the denominator is below inv_mod_th.
#' @param inv_mod_th (Optional) A float indicating the upper threshold of the  denominator of the scaling factor (Equation 4 of Cannon et al. (2015)) to truncate the scaling factor. This parameter is defined as default as the pp_threshold parameter, described above.
#' @param day_win (Optional) (only for frq='D') An integer indicating how many days to consider backwards and forward to get the statistics of each calendar day.The length of the window will be  (2*win_day-1). For example, day_win=15 -> window of 29. Default: win = 1
#' @param user_pdf (Optional) A flag indicating if the user will define the probability distribution functions (pdf) for the observed and modeled series. The distributions will be the same for all periods and sub-periods. user_pdf = 1 or True: User defines the pdf ; user_pdf = 0 or False: pdf defined by the Kolmogorov -Smirnov test (default)
#' @param pdf_obs (Optional) An integer indicating the probability distribution function (pdf) to be used for the observed data. The pdf will be the same for all periods and sub-periods. Default: None
#' @param pdf_mod (Optional) An integer indicating the probability distribution function (pdf) to be used for the modeled data. The pdf will be the same for all periods and sub-periods. Default: None
#'
#' @return QDM_series: A column vector of data bias corrected with the QDM method. [ndata_mod, 1]
#' @export
#'
#' @examples QDM(obs,mod)
#' @examples QDM(obs,mod,frq='A')
#' @examples QDM(obs,mod,frq='M')
#' @examples QDM(obs,mod,frq='D')
#' @examples QDM(obs,mod,mult_changes=0,allow_negatives=1)
#' @examples QDM(obs,mod,mult_changes=1,allow_negatives=0)
#' @examples QDM(obs,mod,mult_changes=1,allow_negatives=0,frq='D',pp_factor=1/10000,day_win=15)
QDM <- function(obs, mod, mult_change, allow_negatives, frq, pp_threshold,
                pp_factor ,rel_change_th, inv_mod_th, day_win, user_pdf,
                pdf_obs, pdf_mod){

  # Check for optional arguments
  if(missing(mult_change)) {
    mult_change <- 1
  }
  
  if(missing(allow_negatives)) {
    allow_negatives <- 1
  }
  
  if(missing(frq)) {
    frq <- 'A'
  }
  
  if(missing(pp_threshold)) {
    pp_threshold <- 1
  }
  
  if(missing(pp_factor)) {
    pp_factor <- 1/100
  }

  if(missing(rel_change_th)) {
    rel_change_th <- 2
  }
  
  if(missing(inv_mod_th)) {
    # Sets the threshold of the modeled series to compute de delta of the
    # quantile as pp_threshold.
    inv_mod_th <- pp_threshold
  }
  
  if(missing(day_win)) {
    day_win <- 1
  }
  
  if(missing(user_pdf)) {
    user_pdf <- FALSE
  }
  
  if(missing(pdf_obs)) {
    pdf_obs <- FALSE
  }
  
  if(missing(pdf_mod)) {
    pdf_mod <- FALSE
  }
  
  # 0) Only for daily frequency data that does not allow negatives (as
  #    precipitacion). Get a pp_threshold for the modeled series that result
  #    in the same amount of rain-day values in the historical period as in
  #    the observed series
  if (frq=='D' & allow_negatives==FALSE) {
    pp_threshold_mod <- get_pp_threshold_mod(obs, mod, pp_threshold)
  } else {
    pp_threshold_mod <- pp_threshold
  }

  # 1) Format series to matrix with rows as sub-periods and columns as years
  #    and, if needed, replace no-rain values with random small values 
  format_list_obs <- formatQM(obs, allow_negatives, frq, pp_threshold, pp_factor)
  y_obs <- format_list_obs[[1]]
  obs_series <- format_list_obs[[2]]
  
  format_list_mod <- formatQM(mod, allow_negatives, frq, pp_threshold_mod, pp_factor)
  y_mod <- format_list_mod[[1]]
  mod_series <- format_list_mod[[2]]
  
  # 2) Get statistics of the observed and modeled series in the historical
  #    period.
  if (frq=='D') {
    # Get a 3D array with a moving window centered on each day
    obs_series_moving <- day_centered_moving_window(obs_series, day_win)
    mod_series_moving <- day_centered_moving_window(mod_series, day_win)
    
    # Get a 4D array with a backward moving window for each projected period
    win_series <- projected_backward_moving_window(mod_series_moving, y_obs, frq)
    
    # Reshape to 2D in order to have all the days of the window in each row
    obs_series_moving <- array(aperm(obs_series_moving, c(1,3,2)), c(365,(2*day_win-1)*y_obs))
    modh_series_moving <- array(aperm(mod_series_moving,c(1,3,2))[,1:y_obs,], c(365,(2*day_win-1)*y_obs))
    
    # Reshape to 3D in order to have all the days of the window in each row
    win_series_moving <- array(aperm(win_series,c(1,4,2,3)), c(365, (2*day_win-1)*y_obs, y_mod-y_obs))
    
    # Replace no-rain values with nans but keep at least 30 values to get
    # the statistics. If fewer than 30 values are rain values, the missing
    # ones will be previously random values replaced
    if (allow_negatives==FALSE) {
      obs_series_moving <- set_norain_to_nan(obs_series_moving, pp_threshold, pp_factor)
      modh_series_moving <- set_norain_to_nan(modh_series_moving, pp_threshold_mod, pp_factor)
      mod_series[mod_series<pp_threshold_mod] <- NaN
      win_series_moving <- set_norain_to_nan(win_series_moving, pp_threshold_mod, pp_factor)
    }
    
    # Get the statistics for each row
    stats_list_obs <- getStats(obs_series_moving, frq)
    mu_obs <- stats_list_obs[[1]]
    std_obs <- stats_list_obs[[2]]
    skew_obs <- stats_list_obs[[3]]
    skewy_obs <- stats_list_obs[[4]]
    
    stats_list_mod <- getStats(modh_series_moving, frq)
    mu_mod <- stats_list_mod[[1]]
    std_mod <- stats_list_mod[[2]]
    skew_mod <- stats_list_mod[[3]]
    skewy_mod <- stats_list_mod[[4]]
    
    stats_list_win <- getStats(win_series_moving, frq)
    mu_win <- stats_list_win[[1]]
    std_win <- stats_list_win[[2]]
    skew_win <- stats_list_win[[3]]
    skewy_win <- stats_list_win[[4]]
    
  } else {
    # Get a 3D array with a backward moving window for each projected period
    win_series <- projected_backward_moving_window(mod_series, y_obs, frq)
    
    # Get the statistics for each row
    stats_list_obs <- getStats(obs_series, frq)
    mu_obs <- stats_list_obs[[1]]
    std_obs <- stats_list_obs[[2]]
    skew_obs <- stats_list_obs[[3]]
    skewy_obs <- stats_list_obs[[4]]
    
    stats_list_mod <- getStats(matrix(mod_series[,1:y_obs],nrow=dim(mod_series)[1]), frq)
    mu_mod <- stats_list_mod[[1]]
    std_mod <- stats_list_mod[[2]]
    skew_mod <- stats_list_mod[[3]]
    skewy_mod <- stats_list_mod[[4]]
    
    stats_list_win <- getStats(win_series, frq)
    mu_win <- stats_list_win[[1]]
    std_win <- stats_list_win[[2]]
    skew_win <- stats_list_win[[3]]
    skewy_win <- stats_list_win[[4]]
  }

  # 3) Assign a probability distribution function to each sub-period for the
  #    observed and modeled data in the historical period.
  if (user_pdf==FALSE){
    if (frq=='D') {
      list_getdist <- getDist(obs_series_moving, allow_negatives, mu_obs,std_obs,skew_obs,skewy_obs)
      pdf_obs <- list_getdist[[1]] 
      ks_fail_obs <- list_getdist[[2]]
      
      list_getdist <- getDist(modh_series_moving, allow_negatives, mu_mod,std_mod,skew_mod,skewy_mod)
      pdf_mod <- list_getdist[[1]] 
      ks_fail_mod <- list_getdist[[2]] 
    } else {
      list_getdist <- getDist(obs_series, allow_negatives, mu_obs,std_obs,skew_obs,skewy_obs)
      pdf_obs <- list_getdist[[1]] 
      ks_fail_obs <- list_getdist[[2]]
      
      list_getdist <- getDist(matrix(mod_series[,1:y_obs],nrow=dim(mod_series)[1]),allow_negatives, mu_mod,std_mod,skew_mod,skewy_mod)
      pdf_mod <- list_getdist[[1]] 
      ks_fail_mod <- list_getdist[[2]] 
    }
    
  } else {
    pdf_obs <- array(0, dim(obs_series)[1]) + pdf_obs
    pdf_mod <- array(0, dim(mod_series)[1]) + pdf_mod
  }
  
  # 4) For each projected period assign a probability distribution function to
  #    to each sub-period and apply the cumulative distribution function to
  #    the modeled data.
  ks_fail_win <- 0
  pdf_win <- matrix(0,dim(mod_series)[1],dim(mod_series)[2]-y_obs)
  prob <- matrix(0,dim(mod_series)[1],dim(mod_series)[2]-y_obs)
  for (j in 1:dim(prob)[2]){
    # Assign a probability distribution function to each month.
    if (user_pdf==FALSE){
      if (frq=='D') {
        list_getdist <- getDist(array(win_series[,,j,], c(365, (2*day_win-1)*y_obs)),
                               allow_negatives,mu_win[,j],std_win[,j],skew_win[,j],skewy_win[,j])
        pdf_win[,j] <-list_getdist[[1]] 
        ks_fail_wtemp <- list_getdist[[2]]
      } else {
        list_getdist <- getDist(matrix(win_series[,j,], nrow=dim(win_series)[1], ncol=dim(win_series)[3]),allow_negatives,mu_win[,j],std_win[,j],skew_win[,j],skewy_win[,j])
        pdf_win[,j] <-list_getdist[[1]] 
        ks_fail_wtemp <- list_getdist[[2]]
      }
      ks_fail_win <- ks_fail_win + ks_fail_wtemp
      
    } else {
      pdf_win[,j] <- pdf_mod
    }
    
    # Apply the cumulative distribution function of the projected period,
    # evaluated with the statistics of this period, to the last data of the
    # period. Equation 3 in Cannon et al. (2015).
    if (frq=='D') {
      prob[,j] <- getCDF(pdf_win[,j],matrix(win_series[,day_win,j,y_obs]),mu_win[,j],std_win[,j],skew_win[,j],skewy_win[,j])
      prob[is.nan(matrix(win_series[,day_win,j,y_obs])),j] <- NaN
    } else {
      prob[,j] <- getCDF(pdf_win[,j],matrix(mod_series[,y_obs+j]),mu_win[,j],std_win[,j],skew_win[,j],skewy_win[,j])
    }
  }

  # 5) Apply the inverse cumulative distribution function:
  #    a) Of the observed data, evaluated with the statistics of the observed
  #       data in the historical period, to the probabilities obtained from 
  #       4b). Equation 5 of Cannon et al. (2015).
  inv_obs <- getCDFinv(pdf_obs,prob,mu_obs,std_obs,skew_obs,skewy_obs)

  #    b) Of the modeled data, evaluated with the statistics of the observed
  #       data in the historical period, to the probabilities obtained from 
  #       4b). Equations 4 of Cannon et al. (2015).
  inv_mod <- getCDFinv(pdf_mod,prob,mu_mod,std_mod,skew_mod,skewy_mod)

  # 6) Get the delta factor or relative change and apply it to the value
  #    obtained in 4b). Equation 4 and 6 of Cannon et al. (2015).
  if mult_change:
  if (mult_change == 1){
    delta_quantile <- mod_series[,(y_obs+1):dim(mod_series)[2]]/inv_mod
    bool_undefined <- (delta_quantile>rel_change_th) & (inv_mod<inv_mod_th)
    delta_quantile[bool_undefined] <- rel_change_th
    QDM <- inv_obs*delta_quantile
  } else{
    delta_quantile <- mod_series[,(y_obs+1):dim(mod_series)[2]] - inv_mod
    QDM <- inv_obs + delta_quantile
  }

  # 7) Perform QM for the historical period.
  mod_h <- mod[1:length(obs)]
  mod_h <- matrix(mod_h)
  QM_series <- QM(obs, mod_h, allow_negatives,frq,pp_threshold, pp_factor, day_win, user_pdf, pdf_obs, pdf_mod)
  
  # 8) Reshape to a column vector
  QDM <- matrix(QDM)
  
  # 9) Concat QM (hsitorical) and QDM (future)
  QDM_series <- matrix(c(QM_series,QDM))
  
  # 10) If does not allow negative values, replace nans (should be the result
  #     of correcting removed no-rain values) and replace no-rain values with 0
  if (allow_negatives == 0){
    QDM_series[is.nan(QDM_series)] <- 0
    QDM_series[QDM_series<pp_threshold] <- 0
  }
  
  if (user_pdf==FALSE){
    ks_fail <- ks_fail_obs + ks_fail_mod + ks_fail_win
    if (ks_fail>0) {
      print('QDM: Some of the probability distribution functions did not pass the KS-Test')
    }
  }
  
  return(QDM_series)
}

#' Unbiased Quantile Mapping
#'
#' This function performs bias correction of modeled series based on observed data by the Unbiased Quantile Mapping (UQM) method, as described by Chadwick et al. (2023). Correction is performed to daily, monthly or annual data in a single location. An independent probability distribution function is assigned to each sub-annual period based on the Kolmogorov-Smirnov test. Each projected period considers a window whose length is equal to the number of years of the historical period and ends in the analyzed year.
#' 
#' Correction of the historical period is performed by the Quantile Mapping (QM) method, as described in Cannon et al. (2015) and Chadwick et al. (2023).
#'
#' The available distributions are: 1) Normal ; 2) Log-Normal ; 3) Gamma 2 parameters ; 4) Gamma 3 parameters ; (Pearson 3 parameters) ; 5) Log-Gamma 3 parameters ; 6) Gumbel ; 7) Exponential
#'
#' NOTE: This routine considers that obs and mod series start in the same day/month/year and are continuous until the end day/month/year.
#' 
#' Chadwick, C., Gironás, J., González-Leiva, F., and Aedo, S. (2023). Bias adjustment to preserve changes in variability: the unbiased mapping of GCM changes. Hydrological Sciences Journal, https://doi.org/10.1080/02626667.2023.2201450
#'
#' @param obs A column vector of daily, monthly or annual observed data, without considering leap days. If daily frequency is specified, the length of the column vector should by a multiple of 365 and for monthly frequency, it should be a multiple of 12. [ndata_obs, 1]
#' @param mod A column vector of daily, monthly or annual modeled or GCM data, without considering leap days. If daily frequency is specified, the length of the column vector should by a multiple of 365 and for monthly frequency, it should be a multiple of 12. [ndata_mod, 1]
#' @param mult_change (Optional) A flag that indicates if projected changes should be computed as multiplicative (fut = hist*delta) or  additive (fut = hist + delta) changes. mult_change = 1 or True: Multiplicative (default) ; mult_change = 0 or False: Additive
#' @param allow_negatives (Optional) A flag that identifies if data allows negative values and also to replace no-rain values with random small values (Chadwick et al., 2023) to avoid numerical problems with the probability distribution functions. allow_negatives = 1 or True: Allow negatives (default) ; allow_negatives = 0 or False: Do not allow negative
#' @param frq (Optional) A string specifying if the input frequency is daily, monthly or annual. Daily:     frq = 'D' ; Monthly:   frq = 'M' ; Annual:    frq = 'A' (default)
#' @param pp_threshold (Optional) A float indicating the threshold to consider no-rain values. Default is 1.
#' @param pp_factor (Optional) A float which multiplied to pp_threshold indicates the maximum value of the random values that replace no-rain values. Default is 1/100.
#' @param day_win (Optional) (only for frq='D') An integer indicating how many days to consider backwards and forward to get the statistics of each calendar day.The length of the window will be  (2*win_day-1). For example, day_win=15 -> window of 29. Default: win = 1
#' @param user_pdf (Optional) A flag indicating if the user will define the probability distribution functions (pdf) for the observed and modeled series. The distributions will be the same for all periods and sub-periods. user_pdf = 1 or True: User defines the pdf ; user_pdf = 0 or False: pdf defined by the Kolmogorov -Smirnov test (default)
#' @param pdf_obs (Optional) An integer indicating the probability distribution function (pdf) to be used for the observed data. The pdf will be the same for all periods and sub-periods. Default: None
#' @param pdf_mod (Optional) An integer indicating the probability distribution function (pdf) to be used for the modeled data. The pdf will be the same for all periods and sub-periods. Default: None
#'
#' @return UQM_series: A column vector of data bias corrected with the UQM method. [ndata_mod, 1]
#' @export
#'
#' @examples UQM(obs,mod)
#' @examples UQM(obs,mod,frq='A')
#' @examples UQM(obs,mod,frq='M')
#' @examples UQM(obs,mod,frq='D')
#' @examples UQM(obs,mod,mult_changes=0,allow_negatives=1)
#' @examples UQM(obs,mod,mult_changes=1,allow_negatives=0)
#' @examples UQM(obs,mod,mult_changes=1,allow_negatives=0,frq='D',pp_factor=1/10000,day_win=15)
UQM <- function(obs, mod, mult_change, allow_negatives, frq, pp_threshold,
                pp_factor, day_win, user_pdf, pdf_obs, pdf_mod){

  # Check for optional arguments
  if(missing(mult_change)) {
    mult_change <- 1
  }
  
  if(missing(allow_negatives)) {
    allow_negatives <- 1
  }
  
  if(missing(frq)) {
    frq <- 'A'
  }
  
  if(missing(pp_threshold)) {
    pp_threshold <- 1
  }
  
  if(missing(pp_factor)) {
    pp_factor <- 1/100
  }
  
  if(missing(day_win)) {
    day_win <- 1
  }
  
  if(missing(user_pdf)) {
    user_pdf <- FALSE
  }
  
  if(missing(pdf_obs)) {
    pdf_obs <- FALSE
  }
  
  if(missing(pdf_mod)) {
    pdf_mod <- FALSE
  }
  
  # 0) Only for daily frequency data that does not allow negatives (as
  #    precipitacion). Get a pp_threshold for the modeled series that result
  #    in the same amount of rain-day values in the historical period as in
  #    the observed series
  if (frq=='D' & allow_negatives==FALSE) {
    pp_threshold_mod <- get_pp_threshold_mod(obs, mod, pp_threshold)
  } else {
    pp_threshold_mod <- pp_threshold
  }
  
  # 1) Format series to matrix with rows as sub-periods and columns as years
  #    and, if needed, replace no-rain values with random small values 
  format_list_obs <- formatQM(obs, allow_negatives, frq, pp_threshold, pp_factor)
  y_obs <- format_list_obs[[1]]
  obs_series <- format_list_obs[[2]]
  
  format_list_mod <- formatQM(mod, allow_negatives, frq, pp_threshold_mod, pp_factor)
  y_mod <- format_list_mod[[1]]
  mod_series <- format_list_mod[[2]]
  
  # 2) Get statistics of the observed and modeled series in the historical
  #    period.
  if (frq=='D') {
    # Get a 3D array with a moving window centered on each day
    obs_series_moving <- day_centered_moving_window(obs_series, day_win)
    mod_series_moving <- day_centered_moving_window(mod_series, day_win)
    
    # Get a 4D array with a backward moving window for each projected period
    win_series <- projected_backward_moving_window(mod_series_moving, y_obs, frq)
    
    # Reshape to 2D in order to have all the days of the window in each row
    obs_series_moving <- array(aperm(obs_series_moving, c(1,3,2)), c(365,(2*day_win-1)*y_obs))
    modh_series_moving <- array(aperm(mod_series_moving,c(1,3,2))[,1:y_obs,], c(365,(2*day_win-1)*y_obs))
    
    # Reshape to 3D in order to have all the days of the window in each row
    win_series_moving <- array(aperm(win_series,c(1,4,2,3)), c(365, (2*day_win-1)*y_obs, y_mod-y_obs))
    
    # Replace no-rain values with nans but keep at least 30 values to get
    # the statistics. If fewer than 30 values are rain values, the missing
    # ones will be previously random values replaced
    if (allow_negatives==FALSE) {
      obs_series_moving <- set_norain_to_nan(obs_series_moving, pp_threshold, pp_factor)
      modh_series_moving <- set_norain_to_nan(modh_series_moving, pp_threshold_mod, pp_factor)
      mod_series[mod_series<pp_threshold_mod] <- NaN
      win_series_moving <- set_norain_to_nan(win_series_moving, pp_threshold_mod, pp_factor)
    }
    
    # Get the statistics for each row
    stats_list_obs <- getStats(obs_series_moving, frq)
    mu_obs <- stats_list_obs[[1]]
    std_obs <- stats_list_obs[[2]]
    skew_obs <- stats_list_obs[[3]]
    skewy_obs <- stats_list_obs[[4]]
    
    stats_list_mod <- getStats(modh_series_moving, frq)
    mu_mod <- stats_list_mod[[1]]
    std_mod <- stats_list_mod[[2]]
    skew_mod <- stats_list_mod[[3]]
    skewy_mod <- stats_list_mod[[4]]
    
    stats_list_win <- getStats(win_series_moving, frq)
    mu_win <- stats_list_win[[1]]
    std_win <- stats_list_win[[2]]
    skew_win <- stats_list_win[[3]]
    skewy_win <- stats_list_win[[4]]
    
  } else {
    # Get a 3D array with a backward moving window for each projected period
    win_series <- projected_backward_moving_window(mod_series, y_obs, frq)
    
    # Get the statistics for each row
    stats_list_obs <- getStats(obs_series, frq)
    mu_obs <- stats_list_obs[[1]]
    std_obs <- stats_list_obs[[2]]
    skew_obs <- stats_list_obs[[3]]
    skewy_obs <- stats_list_obs[[4]]
    
    stats_list_mod <- getStats(matrix(mod_series[,1:y_obs],nrow=dim(mod_series)[1]), frq)
    mu_mod <- stats_list_mod[[1]]
    std_mod <- stats_list_mod[[2]]
    skew_mod <- stats_list_mod[[3]]
    skewy_mod <- stats_list_mod[[4]]
    
    stats_list_win <- getStats(win_series, frq)
    mu_win <- stats_list_win[[1]]
    std_win <- stats_list_win[[2]]
    skew_win <- stats_list_win[[3]]
    skewy_win <- stats_list_win[[4]]
  }

  # 3) Assign a probability distribution function to each sub-period for the
  #    observed and modeled data in the historical period.
  if (user_pdf==FALSE){
    if (frq=='D') {
      list_getdist <- getDist(obs_series_moving, allow_negatives, mu_obs,std_obs,skew_obs,skewy_obs)
      pdf_obs <- list_getdist[[1]] 
      ks_fail_obs <- list_getdist[[2]]
      
    } else {
      list_getdist <- getDist(obs_series, allow_negatives, mu_obs,std_obs,skew_obs,skewy_obs)
      pdf_obs <- list_getdist[[1]] 
      ks_fail_obs <- list_getdist[[2]]

    }
    
  } else {
    pdf_obs <- array(0, dim(obs_series)[1]) + pdf_obs
    pdf_mod <- array(0, dim(mod_series)[1]) + pdf_mod
  }

  mu_mod_repeated <- array(rep(mu_mod,(y_mod-y_obs)),c(dim(mod_series)[1],(y_mod-y_obs)))
  std_mod_repeated <- array(rep(std_mod,(y_mod-y_obs)),c(dim(mod_series)[1],(y_mod-y_obs)))
  skew_mod_repeated <- array(rep(skew_mod,(y_mod-y_obs)),c(dim(mod_series)[1],(y_mod-y_obs)))
  skewy_mod_repeated <- array(rep(skewy_mod,(y_mod-y_obs)),c(dim(mod_series)[1],(y_mod-y_obs)))
  
  mu_obs_repeated <- array(rep(mu_obs,(y_mod-y_obs)),c(dim(mod_series)[1],(y_mod-y_obs)))
  std_obs_repeated <- array(rep(std_obs,(y_mod-y_obs)),c(dim(mod_series)[1],(y_mod-y_obs)))
  skew_obs_repeated <- array(rep(skew_obs,(y_mod-y_obs)),c(dim(mod_series)[1],(y_mod-y_obs)))
  skewy_obs_repeated <- array(rep(skewy_obs,(y_mod-y_obs)),c(dim(mod_series)[1],(y_mod-y_obs)))
  
  # 4) For each projected period, get the delta factor (delta) and time
  #    dependent (projected) statistics (mean, standard deviation, skewness, 
  #    and log skewness). Equations 13 to 14 in Chadwick et al. (2023).
  if (mult_change==1){ # Precipitation
    delta_mu <- mu_win/mu_mod_repeated
    delta_sigma <- std_win/std_mod_repeated
    delta_skew <- skew_win/skew_mod_repeated
    delta_skewy <- skewy_win/skewy_mod_repeated
    
    mu_projected <- mu_obs_repeated*delta_mu
    sigma_projected <- std_obs_repeated*delta_sigma
    skew_projected <- skew_obs_repeated*delta_skew
    skewy_projected <- skewy_obs_repeated*delta_skewy
    
  } else { # Temperature
    delta_mu <- mu_win - mu_mod_repeated
    delta_sigma <- std_win - std_mod_repeated
    delta_skew <- skew_win - skew_mod_repeated
    delta_skewy <- skewy_win - skewy_mod_repeated
    
    mu_projected <- mu_obs_repeated + delta_mu
    sigma_projected <- std_obs_repeated + delta_sigma
    skew_projected <- skew_obs_repeated + delta_skew
    skewy_projected <- skewy_obs_repeated + delta_skewy
  }

  # 5) For each projected period assign a probability distribution function to
  #    to each sub-period and apply the cumulative distribution function to
  #    the modeled data.
  ks_fail_win <- 0
  pdf_win <- matrix(0,dim(mod_series)[1],dim(mod_series)[2]-y_obs)
  prob <- matrix(0,dim(mod_series)[1],dim(mod_series)[2]-y_obs)
  UQM <- matrix(0,dim(mod_series)[1],dim(mod_series)[2]-y_obs)
  for (j in 1:dim(prob)[2]){
    # Assign a probability distribution function to each month.
    if (user_pdf==FALSE){
      if (frq=='D') {
        list_getdist <- getDist(array(win_series[,,j,], c(365, (2*day_win-1)*y_obs)),
                                allow_negatives,mu_win[,j],std_win[,j],skew_win[,j],skewy_win[,j])
        pdf_win[,j] <-list_getdist[[1]] 
        ks_fail_wtemp <- list_getdist[[2]]
      } else {
        list_getdist <- getDist(matrix(win_series[,j,], nrow=dim(win_series)[1], ncol=dim(win_series)[3]),allow_negatives,mu_win[,j],std_win[,j],skew_win[,j],skewy_win[,j])
        pdf_win[,j] <-list_getdist[[1]] 
        ks_fail_wtemp <- list_getdist[[2]]
      }
      ks_fail_win <- ks_fail_win + ks_fail_wtemp
      
    } else {
      pdf_win[,j] <- pdf_mod
    }
    
    # Apply the cumulative distribution function of the projected period,
    # evaluated with the statistics of this period, to the last data of the
    # period. Equation 3 in Chadwick et al. (2023).
    if (frq=='D') {
      prob[,j] <- getCDF(pdf_win[,j],matrix(win_series[,day_win,j,y_obs]),mu_win[,j],std_win[,j],skew_win[,j],skewy_win[,j])
    } else {
      prob[,j] <- getCDF(pdf_win[,j],matrix(mod_series[,y_obs+j]),mu_win[,j],std_win[,j],skew_win[,j],skewy_win[,j])
    }
    
    # 6) Apply the inverse cumulative distribution function of the observed
    #    data, evaluated with the time dependent statistics, to the values
    #    obtained in 3) (getCDFinv function of the climQMBC package). Equation
    #    15 in Chadwick et al. (2023).
    UQM[,j] <- getCDFinv(pdf_obs,matrix(prob[,j]),matrix(mu_projected[,j]),matrix(sigma_projected[,j]),matrix(skew_projected[,j]),matrix(skewy_projected[,j]))
  }

  # 7) Perform QM for the historical period.
  mod_h <- mod[1:length(obs)]
  mod_h <- matrix(mod_h)
  QM_series <- QM(obs, mod_h, allow_negatives,frq,pp_threshold, pp_factor, day_win, user_pdf, pdf_obs, pdf_mod)
  
  # 8) Reshape to a column vector
  UQM <- matrix(UQM)
  
  # 9) Concat QM (hsitorical) and UQM (future)
  UQM_series <- matrix(c(QM_series,UQM))
  
  # 10) If does not allow negative values, replace nans (should be the result
  #     of correcting removed no-rain values) and replace no-rain values with 0
  if (allow_negatives == 0){
    UQM_series[is.nan(UQM_series)] <- 0
    UQM_series[UQM_series<pp_threshold] <- 0
  }
  
  if (user_pdf==FALSE){
    ks_fail <- ks_fail_obs + ks_fail_win
    if (ks_fail>0) {
      print('UQM: Some of the probability distribution functions did not pass the KS-Test')
    }
  }
  
  return(UQM_series)
}


#' Scaled Distribution Mapping
#'
#' This function performs bias correction of modeled series based on observed data by the Scaled Distribution Mapping (SDM) method, as described by Switanek et al. (2017). Correction is performed to daily, monthly or annual precipitation or temperature data in a single location. The historical period of the modeled series is bias corrected as a scenario where the complete series is replaced by the bias corrected series. On the other hand, for each year after the historical period the method considers a projected period as a window whose length is equal to the number of years of the historical period and ends in the analyzed year. From the bias corrected series, only the last year is saved and all the others are discarded. This approach aims to build a continuum for the bias corrected series.
#'
#' NOTE: This routine considers that obs and mod series start in the same day/month/year and are continuous until the end day/month/year.
#'
#' Switanek, B. M., Troch, P., Castro, L. C., Leuprecht, A., Chang, H. I., Mukherjee, R., and Demaria, M. C. E. (2017) Scaled distribution mapping: A bias correction method that preserves raw climate model projected changes. Hydrology &amp; Earth System Sciences, 21, 2649-2666, https://doi.org/10.5194/hess-21-2649-2017.
#'
#' @param obs A column vector of daily, monthly or annual observed data, without considering leap days. If daily frequency is specified, the length of the column vector should by a multiple of 365 and for monthly frequency, it should be a multiple of 12. [ndata_obs, 1]
#' @param mod A column vector of daily, monthly or annual modeled or GCM data, without considering leap days. If daily frequency is specified, the length of the column vector should by a multiple of 365 and for monthly frequency, it should be a multiple of 12. [ndata_mod, 1]
#' @param SDM_var A flag that identifies if data are temperature or precipitation. Temperature:   var = 0 ; Precipitation: var = 1
#' @param frq (Optional) A string specifying if the input frequency is daily, monthly or annual. ; Daily:     frq = 'D' ; Monthly:   frq = 'M' ; Annual:    frq = 'A' (default)
#' @param pp_threshold (Optional) A float indicating the threshold to consider no-rain values. Default is 1.
#' @param pp_factor (Optional) A float which multiplied to pp_threshold indicates the maximum value of the random values that replace no-rain values. Default is 1/100.
#' @param day_win (Optional) (only for frq='D') An integer indicating how many days to consider backwards and forward to get the statistics of each calendar day.The length of the window will be  (2*win_day-1). For example, day_win=15 -> window of 29. Default: win = 1
#'
#' @return SDM_series: A column vector of data bias corrected with the QM method. [ndata_mod, 1]
#' @export
#'
#' @examples SDM(obs,mod,SDM_var)
#' @examples SDM(obs,mod,SDM_var,frq='A')
#' @examples SDM(obs,mod,SDM_var,frq='M')
#' @examples SDM(obs,mod,SDM_var,frq='M',pp_threshold=0.1)
#' @examples SDM(obs,mod,SDM_var,frq='M',pp_factor=1/1000)
#' @examples SDM(obs,mod,SDM_var,frq='D',pp_threshold=0.1,pp_factor=1/1000,day_win=15)
SDM <- function(obs ,mod, SDM_var, frq, pp_threshold, pp_factor, day_win){
  
  # Check for optional arguments
  if(missing(frq)) {
    frq <- 'A'
  }
  
  if(missing(pp_threshold)) {
    pp_threshold <- 1
  }

  if(missing(pp_factor)) {
    pp_factor <- 1/100
  }
  
  if(missing(day_win)) {
    day_win <- 1
  }

  # Define the no-rain value and the threshold to truncate the tail of the
  # probability ditribution functions
  lower_lim <- pp_threshold
  CDF_th <- 10^-3

  # 0) Format series to matrix with rows as sub-periods and columns as years
  #    and, if needed, replace no-rain values with random small values 
  format_list_obs <- formatQM(obs, SDM_var, frq, pp_threshold, pp_factor)
  y_obs <- format_list_obs[[1]]
  obs_series <- format_list_obs[[2]]
  
  format_list_mod <- formatQM(mod, SDM_var, frq, pp_threshold, pp_factor)
  y_mod <- format_list_mod[[1]]
  mod_series <- format_list_mod[[2]]
  
  modh_series <- array(mod_series[,1:y_obs], c(dim(mod_series)[1],y_obs))
  
  # If frequency is daily, consider a moving window for each day to fit
  # the statistics
  if (frq=='D') {
    # Get a 3D array with a moving window centered on each day
    obs_series <- day_centered_moving_window(obs_series, day_win)
    mod_series <- day_centered_moving_window(mod_series, day_win)
    
    # Get a 4D array with a backward moving window for each projected period
    win_series <- projected_backward_moving_window(mod_series, y_obs, frq)
    
    # Reshape to 2D in order to have all the days of the window in each row
    obs_series <- array(aperm(obs_series, c(1,3,2)), c(365,(2*day_win-1)*y_obs))
    modh_series <- array(aperm(mod_series,c(1,3,2))[,1:y_obs,], c(365,(2*day_win-1)*y_obs))
    
    # Reshape to 3D in order to have all the days of the window in each row
    win_series <- array(aperm(win_series,c(1,4,2,3)), c(365, (2*day_win-1)*y_obs, y_mod-y_obs))
    
    # Add the historical period to the moving window
    win_series_ <- abind::abind(modh_series, win_series, along=3)
    
  } else {
    # For non-daily frequency, make sure that day_win=1
    day_win <- 1
    
    # Get a 3D array with a backward moving window for each projected period
    win_series <- projected_backward_moving_window(mod_series, y_obs, frq)
    
    # Add the historical period to the moving window
    win_series_ <- abind::abind(modh_series, aperm(win_series,c(1,3,2)), along=3)
  }

  SDM  <- matrix(0,dim(mod_series)[1],y_mod-y_obs)
  SDM_h  <- matrix(0,dim(obs_series)[1],y_obs)
  for (m in 1:dim(mod_series)[1]){
    # 1) Historical period:

    # a) [Switanek et al. (2017), step 1)]
    #     For temperature, get the detrended modeled and observed series in
    #     the historical period.
    #     For precipitation, get the rainday values and its frequency for
    #     the modeled and observed series in the historical period.
    if (SDM_var == 0){
      D_obs <- pracma::detrend(obs_series[m,])
      D_mod <- pracma::detrend(modh_series[m,])

      mu_obs <- mean(obs_series[m,])
      mu_mod <- mean(modh_series[m,])
    } else{
      D_obs <- sort(obs_series[m,obs_series[m,]>lower_lim])
      D_mod <- sort(modh_series[m,][modh_series[m,]>lower_lim])

      freq_obs <- length(D_obs)/dim(obs_series)[2]
      freq_mod <- length(D_mod)/length(modh_series[m,])
      
      if (freq_mod==0) {
        freq_mod <- 1/(365*y_obs)
      }

    }

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
    if (SDM_var == 0){
      mu_obsD <- mean(D_obs)
      mu_modD <- mean(D_mod)
      sigma_obsD <- sd(D_obs)
      sigma_modD <- sd(D_mod)

      CDF_obs <- pnorm(sort(D_obs),mu_obsD,sigma_obsD)
      CDF_mod <- pnorm(sort(D_mod),mu_modD,sigma_modD)
      CDF_obs[CDF_obs<CDF_th] <- CDF_th
      CDF_mod[CDF_mod<CDF_th] <- CDF_th
      CDF_obs[CDF_obs>1-CDF_th] <- 1-CDF_th
      CDF_mod[CDF_mod>1-CDF_th] <- 1-CDF_th
    } else {
      if (length(D_mod)>0){
        fit_mod <- fitdistrplus::fitdist(D_mod, distr = "gamma", method = "mle")
      } else {
        fit_mod <- fitdistrplus::fitdist(runif(2)*pp_threshold*pp_factor, distr = "gamma", method = "mle")
      }
      
      if (length(D_obs)>0){
        fit_obs <- fitdistrplus::fitdist(D_obs, distr = "gamma", method = "mle")
      } else {
        fit_obs <- fitdistrplus::fitdist(runif(2)*pp_threshold*pp_factor, distr = "gamma", method = "mle")
      }

      CDF_obs <- pgamma(D_obs,shape=fit_obs[[1]][1],rate=fit_obs[[1]][2])
      CDF_mod <- pgamma(D_mod,shape=fit_mod[[1]][1],rate=fit_mod[[1]][2])
      CDF_obs[CDF_obs>1-CDF_th] <- 1-CDF_th
      CDF_mod[CDF_mod>1-CDF_th] <- 1-CDF_th
    }

    # 2) Projected periods:
    for (j in 0:(y_mod-y_obs)){
      # c) Initialize correction array.
      corr_temp <- matrix(0,y_obs*(2*day_win-1),1)

      # d) Define projected window.
      win_series <- win_series_[m,,j+1]

      # e) [Switanek et al. (2017), step 1)]
      #    For temperature, get the detrended series of the projected
      #    period.
      #    For precipitation, get the rainday values, its frequency, and
      #    expected raindays for the projected period.
      #    Get the index of the sorted detrended or rainday values.
      if (SDM_var == 0){
        D_win <- pracma::detrend(win_series)
        exp_D <- length(win_series)
        win_argsort <- sort(D_win,index.return=TRUE)[[2]]
        mu_win <- mean(win_series)
      } else {
        D_win <- sort(win_series[win_series>lower_lim])
        exp_D <- min(round(length(D_win)*freq_obs/freq_mod),length(win_series))

        win_argsort <- sort(win_series,index.return=TRUE)[[2]]
      }

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
      if (SDM_var == 0){
        mu_winD <- mean(D_win,na.rm=TRUE)
        sigma_winD <- sd(D_win,na.rm=TRUE)

        diff_win <- win_series - D_win

        CDF_win <- pnorm(sort(D_win),mu_winD,sigma_winD)
        CDF_win[CDF_win<CDF_th] <- CDF_th
        CDF_win[CDF_win>1-CDF_th] <- 1-CDF_th
      } else {
        if (length(D_win)>0){
          fit_win <- fitdistrplus::fitdist(D_win, distr = "gamma", method = "mle")
        } else {
          fit_win <- fitdistrplus::fitdist(runif(2)*pp_threshold*pp_factor, distr = "gamma", method = "mle")
        }
        
        CDF_win <- pgamma(D_win,shape=fit_win[[1]][1],rate=fit_win[[1]][2])
        CDF_win[CDF_win>1-CDF_th] <- 1-CDF_th
      }

      # g) [Switanek et al. (2017), step 3)]
      #    Get the scaling between the model projected period and
      #    historical period distributions.
      if (SDM_var == 0){
        SF <- (sigma_obsD/sigma_modD)*(qnorm(CDF_win,mu_winD,sigma_winD) - qnorm(CDF_win,mu_modD,sigma_modD))
      } else {
        SF <- qgamma(CDF_win,shape=fit_win[[1]][1],rate=fit_win[[1]][2])/qgamma(CDF_win,shape=fit_mod[[1]][1],rate=fit_mod[[1]][2])
      }

      # h) Interpolate observed and modeled CDF of the historical period
      #    to the length of the projected period.
      if (length(D_obs)>0) {
        obs_cdf_intpol <- pracma::interp1(pracma::linspace(1,length(D_obs),length(D_obs)),CDF_obs,pracma::linspace(1,length(D_obs),length(D_win)))
      } else {
        obs_cdf_intpol = CDF_win*0
      }
      
      if (length(D_mod)>0) {
        mod_cdf_intpol <- pracma::interp1(pracma::linspace(1,length(D_mod),length(D_mod)),CDF_mod,pracma::linspace(1,length(D_mod),length(D_win)))
      } else {
        mod_cdf_intpol = CDF_win*0
      }


      # i) [Switanek et al. (2017), steps 4 and 5)]
      #    Get recurrence intervals and its scaled version for the
      #    observed, historical modeled and projected period modeled
      #    CDFs.
      if (SDM_var == 0){
        RI_obs <- 1/(0.5 - abs(obs_cdf_intpol - 0.5))
        RI_mod <- 1/(0.5 - abs(mod_cdf_intpol - 0.5))
        RI_win <- 1/(0.5 - abs(CDF_win - 0.5))
        RI_scaled <- RI_obs*RI_win/RI_mod

        CDF_scaled <- sign(obs_cdf_intpol - 0.5)*(1 - 1/RI_scaled)
        CDF_scaled[CDF_scaled<0] = CDF_scaled[CDF_scaled<0] + 1

        CDF_scaled[CDF_scaled<CDF_th] <- CDF_th
        CDF_scaled[CDF_scaled>1-CDF_th] <- 1-CDF_th
      } else {
        RI_obs <- 1/(1-obs_cdf_intpol)
        RI_mod <- 1/(1-mod_cdf_intpol)
        RI_win <- 1/(1-CDF_win)
        RI_scaled <- RI_obs*RI_win/RI_mod;
        RI_scaled[RI_scaled<1] <- 1

        CDF_scaled <- sort(1 - 1/(RI_scaled))
      }

      # j) [Switanek et al. (2017), step 6)]
      #    Get the initial bias corrected values. For precipitation,
      #    these values are interpolated to the length of the expected
      #    raindays.
      if (SDM_var == 0){
        xvals <- qnorm(sort(CDF_scaled),mu_obsD,sigma_obsD) + SF
        xvals <- xvals - mean(xvals) + mu_obs + (mu_win - mu_mod)
      } else {
        xvals <- qgamma(CDF_scaled,shape=fit_obs[[1]][1],rate=fit_obs[[1]][2])*SF
        if (length(D_win) > exp_D){
          xvals <- pracma::interp1(pracma::linspace(1,length(D_win),length(D_win)),xvals,pracma::linspace(1,length(D_win),exp_D))
        } else {
          xvals <- c(matrix(0,1,exp_D - length(D_win)),xvals)
        }
      }

      # k) [Switanek et al. (2017), step 7)]
      #    Bias corrected values are placed back matching the higher bias
      #    corrected values with the higher rainday or detrended values.
      #    For temperature, the trend of the projected period is added
      #    back.
      if (length(xvals)>0){
        corr_temp[win_argsort[(length(win_argsort)-exp_D+1):length(win_argsort)]] <- matrix(xvals)
      }
      if (SDM_var == 0){
        corr_temp <- corr_temp + diff_win - mu_win
      }

      # l) If the projected period is the historical period (period 0,
      #    j=0) save the complete bias corrected series.
      #    If the projected period is not the historical period (j>0),
      #    save the value of the last year.
      if (j == 0){
        SDM_h[m,] <- corr_temp[seq(day_win,length(corr_temp),(2*day_win-1))]
      } else {
        SDM[m,j] <- corr_temp[length(corr_temp)]
      }
    }
  }

  # Reshape to a column vector
  SDM <- matrix(SDM)
  SDM_h <- matrix(SDM_h)
  
  # Concat historical and future
  SDM_series <- matrix(c(SDM_h,SDM))
  
  # If precipitation, replace no-rain values with 0
  if (SDM_var == 1){
    SDM_series[SDM_series<pp_threshold] <- 0
  }

  return(SDM_series)
}
