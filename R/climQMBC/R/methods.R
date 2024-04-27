#' Quantile Mapping
#'
#' This function performs bias correction of modeled series based on observed data by the Quantile Mapping (QM) method, as described by Cannon et al.  (2015). Correction is performed to monthly or annual precipitation or temperature data in a single location. An independent probability distribution function is assigned to each month and to each projected period based on the Kolmogorov-Smirnov test. If annual frequency is specified, this is applied to the complete period. Each projected period  considers a window whose length is equal to the number of years of the historical period and ends in the analyzed year.
#'
#' Cannon, A. J., S. R. Sobie, and T. Q. Murdock, 2015: Bias correction of GCM precipitation by quantile mapping: How well do methods preserve changes in quantiles and extremes? J. Climate, 28(17), 6938-6959, https://doi.org/10.1175/JCLI-D-14-00754.1
#'
#' @param obs A column vector of monthly or annual observed data (temperature or precipitation). If monthly frequency is specified, the length of this vector is 12 times the number of observed years [12 x y_obs, 1]. If annual frequency is specified, the length of this vector is equal to the number of observed years [y_obs, 1].
#' @param mod A column vector of monthly or annual modeled data (temperature or precipitation). If monthly frequency is specified, the length of this vector is 12 times the number of observed years [12 x y_mod, 1]. If annual frequency is specified, the length of this vector is equal to the number of observed years [y_mod, 1].
#' @param var A flag that identifies if data are temperature or precipitation. This flag tells the getDist function if it has to discard distribution functions that allow negative numbers, and if the terms in the correction equations are multiplied/divided or added/subtracted. Temperature:   var = 0; Precipitation: var = 1
#' @param frq (Optional) A string specifying if the input is annual or monthly data. If not specified, it is set monthly as default. Monthly:   frq = 'M'; Annual:    frq = 'A'
#' @param pp_threshold (Optional) A float indicating the threshold to consider physically null precipitation values.
#' @param pp_factor (Optional) A float indicating the maximum value of the random values that replace physically null precipitation values.
#'
#' @return A column vector of monthly or annual modeled data (temperature or precipitation) corrected by the QM method. If monthly frequency is specified, the length of this vector is 12 times the number of observed years [12 x y_mod, 1]. If annual frequency is specified, the length of this vector is equal to the number of observed  years [y_mod, 1].
#' @export
#'
#' @examples QM(obs,mod,var)
#' @examples QM(obs,mod,var,frq='A')
#' @examples QM(obs,mod,var,frq='M')
#' @examples QM(obs,mod,var,frq='M',pp_threshold=0.1)
#' @examples QM(obs,mod,var,frq='M',pp_factor=1/1000)
#' @examples QM(obs,mod,var,frq='M',pp_threshold=0.1,pp_factor=1/1000)
QM <- function(obs, mod, allow_negatives, frq, pp_threshold, pp_factor, day_win,
               user_pdf, pdf_obs, pdf_mod){
  
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
  
  if (frq=='D' & allow_negatives==FALSE) {
    pp_threshold_mod <- get_pp_threshold_mod(obs, mod, pp_threshold)

  } else {
    pp_threshold_mod <- pp_threshold
  }
  
  
  # 1) Format inputs and get statistics of the observed and modeled series of
  #    the historical period (formatQM function of the climQMBC package).
  format_list_obs <- formatQM(obs, allow_negatives, frq, pp_threshold, pp_factor)
  y_obs <- format_list_obs[[1]]
  obs_series <- format_list_obs[[2]]
  
  format_list_mod <- formatQM(mod, allow_negatives, frq, pp_threshold_mod, pp_factor)
  y_mod <- format_list_mod[[1]]
  mod_series <- format_list_mod[[2]]
  
  if (frq=='D') {
    obs_series_moving <- day_centered_moving_window(obs_series, day_win)
    mod_series_moving <- day_centered_moving_window(mod_series, day_win)

    obs_series_moving <- array(aperm(obs_series_moving, c(1,3,2)), c(365,(2*day_win-1)*y_obs))
    modh_series_moving <- array(aperm(mod_series_moving,c(1,3,2))[,1:y_obs,], c(365,(2*day_win-1)*y_obs))
    
    
    if (allow_negatives==FALSE) {
      obs_series_moving <- set_norain_to_nan(obs_series_moving, pp_threshold, pp_factor)
      modh_series_moving <- set_norain_to_nan(modh_series_moving, pp_threshold_mod, pp_factor)
      mod_series[mod_series<pp_threshold_mod] <- NaN
    }
    
    
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

  # 2) Assign a probability distribution function to each month for the
  #    observed and modeled data in the historical period. If annual
  #    frequency is specified, this is applied to the complete historical
  #    period (getDist function of the climQMBC package).
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
  

  # 3) Apply the cumulative distribution function of the modeled data,
  #    evaluated with the statistics of the modeled data in the historical
  #    period, to the modeled data (getCDF function of the climQMBC package).
  #    Equation 1 of Cannon et al. (2015).
  prob <- getCDF(pdf_mod,mod_series,mu_mod,std_mod,skew_mod,skewy_mod)
  prob[is.nan(mod_series)] <- NaN

  # 4) Apply the inverse cumulative distribution function of the observed
  #    data, evaluated with the statistics of the observed data in the
  #    historical period, to the probabilities obtained from 3) (getCDFinv
  #    function of the climQMBC package). Equation 1 of Cannon et al. (2015).
  QM_series <- getCDFinv(pdf_obs,prob,mu_obs,std_obs,skew_obs,skewy_obs)
  QM_series <- matrix(QM_series)
  
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
#' This function performs bias correction of modeled series based on observed data by the Detrended Quantile Mapping (DQM) method, as described by Cannon et al.  (2015). Correction is performed to monthly or annual precipitation or temperature data in a single location. An independent probability distribution function is assigned to each month and to each projected period based on the Kolmogorov-Smirnov test. If annual frequency is specified, this is applied to the complete period. Each projected period  considers a window whose length is equal to the number of years of the historical period and ends in the analyzed year. Correction of the historical period is performed by the Quantile Mapping (QM) method, as described by Cannon et al. (2015).
#'
#' Cannon, A. J., S. R. Sobie, and T. Q. Murdock, 2015: Bias correction of GCM precipitation by quantile mapping: How well do methods preserve changes in quantiles and extremes? J. Climate, 28(17), 6938-6959, https://doi.org/10.1175/JCLI-D-14-00754.1
#'
#' @param obs A column vector of monthly or annual observed data (temperature or precipitation). If monthly frequency is specified, the length of this vector is 12 times the number of observed years [12 x y_obs, 1]. If annual frequency is specified, the length of this vector is equal to the number of observed years [y_obs, 1].
#' @param mod A column vector of monthly or annual modeled data (temperature or precipitation). If monthly frequency is specified, the length of this vector is 12 times the number of observed years [12 x y_mod, 1]. If annual frequency is specified, the length of this vector is equal to the number of observed years [y_mod, 1].
#' @param var A flag that identifies if data are temperature or precipitation. This flag tells the getDist function if it has to discard distribution functions that allow negative numbers, and if the terms in the correction equations are multiplied/divided or added/subtracted. Temperature:   var = 0; Precipitation: var = 1
#' @param frq (Optional) A string specifying if the input is annual or monthly data. If not specified, it is set monthly as default. Monthly:   frq = 'M'; Annual:    frq = 'A'
#' @param pp_threshold (Optional) A float indicating the threshold to consider physically null precipitation values.
#' @param pp_factor (Optional) A float indicating the maximum value of the random values that replace physically null precipitation values.
#'
#' @return A column vector of monthly or annual modeled data (temperature or precipitation) corrected by the DQM method. If monthly frequency is specified, the length of this vector is 12 times the number of observed years [12 x y_mod, 1]. If annual frequency is specified, the length of this vector is equal to the number of observed  years [y_mod, 1].
#' @export
#'
#' @examples DQM(obs,mod,var)
#' @examples DQM(obs,mod,var,frq='A')
#' @examples DQM(obs,mod,var,frq='M')
#' @examples DQM(obs,mod,var,frq='M',pp_threshold=0.1)
#' @examples DQM(obs,mod,var,frq='M',pp_factor=1/1000)
#' @examples DQM(obs,mod,var,frq='M',pp_threshold=0.1,pp_factor=1/1000)
DQM <- function(obs, mod, mult_change, allow_negatives, frq, pp_threshold,
                pp_factor, day_win, user_pdf, pdf_obs, pdf_mod){
  
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
  
  if (frq=='D' & allow_negatives==FALSE) {
    pp_threshold_mod <- get_pp_threshold_mod(obs, mod, pp_threshold)
    
  } else {
    pp_threshold_mod <- pp_threshold
  }
  
  
  # 1) Format inputs and get statistics of the observed and modeled series of
  #    the historical period (formatQM function of the climQMBC package).
  format_list_obs <- formatQM(obs, allow_negatives, frq, pp_threshold, pp_factor)
  y_obs <- format_list_obs[[1]]
  obs_series <- format_list_obs[[2]]
  
  format_list_mod <- formatQM(mod, allow_negatives, frq, pp_threshold_mod, pp_factor)
  y_mod <- format_list_mod[[1]]
  mod_series <- format_list_mod[[2]]
  
  if (frq=='D') {
    
    obs_series_moving <- day_centered_moving_window(obs_series, day_win)
    mod_series_moving <- day_centered_moving_window(mod_series, day_win)
    
    win_series <- projected_backward_moving_window(mod_series_moving, y_obs, frq)

    obs_series_moving <- array(aperm(obs_series_moving, c(1,3,2)), c(365,(2*day_win-1)*y_obs))
    modh_series_moving <- array(aperm(mod_series_moving,c(1,3,2))[,1:y_obs,], c(365,(2*day_win-1)*y_obs))
    win_series_moving <- array(aperm(win_series,c(1,4,2,3)), c(365, (2*day_win-1)*y_obs, y_mod-y_obs))

    
    if (allow_negatives==FALSE) {
      obs_series_moving <- set_norain_to_nan(obs_series_moving, pp_threshold, pp_factor)
      modh_series_moving <- set_norain_to_nan(modh_series_moving, pp_threshold_mod, pp_factor)
      mod_series[mod_series<pp_threshold_mod] <- NaN
      win_series_moving <- set_norain_to_nan(win_series_moving, pp_threshold_mod, pp_factor)
    }

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
    win_series <- projected_backward_moving_window(mod_series, y_obs, frq)
    
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
  

  # 2) Assign a probability distribution function to each month for the
  #    observed and modeled data in the historical period. If annual
  #    frequency is specified, this is applied to the complete historical
  #    period (getDist function of the climQMBC package).
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


  mu_mod_repeated <- array(rep(mu_mod,(y_mod-y_obs)),c(dim(mod_series)[1],(y_mod-y_obs)))
  
  # 4)Compute the linear scaled values (value in square brackets in equation
  #                                     2 of Cannon et al. (2015)).
  if (mult_change == 1){ # Precipitation
    scaling_factor <- mu_win/mu_mod_repeated
    detrended_series <- mod_series[,(y_obs+1):dim(mod_series)[2]]/scaling_factor
  } else{        # Temperature
    scaling_factor <- mu_win - mu_mod_repeated
    detrended_series <- mod_series[,(y_obs+1):dim(mod_series)[2]] - scaling_factor
  }
  

  # 5) Apply the cumulative distribution function of the modeled data,
  #    evaluated with the statistics of the modeled period, to the future
  #    modeled data (getCDF function of the climQMBC package). Equation 2 of
  #    Cannon et al. (2015).
  prob <- getCDF(pdf_mod,detrended_series,mu_mod,std_mod,skew_mod,skewy_mod)
  prob[is.nan(detrended_series)] <- NaN

  # 6) Apply the inverse cumulative distribution function of the observed
  #    data, evaluated with the statistics of the observed data in the
  #    historical period, to the probabilities obtained from 5) (getCDFinv
  #    function of the climQMBC package). Equation 2 of Cannon et al. (2015).
  DQM_LS <- getCDFinv(pdf_obs,prob,mu_obs,std_obs,skew_obs,skewy_obs)

  # 7) Reimpose the trend to the values obtained in 6). Equation 2 of Cannon
  #    et al. (2015).
  if (mult_change == 1){ # Precipitation
    DQM <- DQM_LS*scaling_factor
  } else{        # Temperature
    DQM <- DQM_LS + scaling_factor
  }
  

  DQM <- matrix(DQM)

  # 8) Perform QM for the historical period.
  mod_h <- mod[1:length(obs)]
  mod_h <- matrix(mod_h)
  QM_series <- QM(obs, mod_h, allow_negatives,frq,pp_threshold, pp_factor, day_win, user_pdf, pdf_obs, pdf_mod)
  DQM_series <- matrix(c(QM_series,DQM))
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
#' This function performs bias correction of modeled series based on observed data by the Quantile Delta Mapping (QDM) method, as described by Cannon et al.  (2015). Correction is performed to monthly or annual precipitation or temperature data in a single location. An independent probability distribution function is assigned to each month and to each projected period based on the Kolmogorov-Smirnov test. If annual frequency is specified, this is applied to the complete period. Each projected period  considers a window whose length is equal to the number of years of the historical period and ends in the analyzed year. Correction of the historical period is performed by the Quantile Mapping (QM) method, as described by Cannon et al. (2015).
#'
#' Cannon, A. J., S. R. Sobie, and T. Q. Murdock, 2015: Bias correction of GCM precipitation by quantile mapping: How well do methods preserve changes in quantiles and extremes? J. Climate, 28(17), 6938-6959, https://doi.org/10.1175/JCLI-D-14-00754.1
#'
#' @param obs A column vector of monthly or annual observed data (temperature or precipitation). If monthly frequency is specified, the length of this vector is 12 times the number of observed years [12 x y_obs, 1]. If annual frequency is specified, the length of this vector is equal to the number of observed years [y_obs, 1].
#' @param mod A column vector of monthly or annual modeled data (temperature or precipitation). If monthly frequency is specified, the length of this vector is 12 times the number of observed years [12 x y_mod, 1]. If annual frequency is specified, the length of this vector is equal to the number of observed years [y_mod, 1].
#' @param var A flag that identifies if data are temperature or precipitation. This flag tells the getDist function if it has to discard distribution functions that allow negative numbers, and if the terms in the correction equations are multiplied/divided or added/subtracted. Temperature:   var = 0; Precipitation: var = 1
#' @param frq (Optional) A string specifying if the input is annual or monthly data. If not specified, it is set monthly as default. Monthly:   frq = 'M'; Annual:    frq = 'A'
#' @param pp_threshold (Optional) A float indicating the threshold to consider physically null precipitation values.
#' @param pp_factor (Optional) A float indicating the maximum value of the random values that replace physically null precipitation values.
#' @param rel_change_th (Optional) A float indicating the maximum scaling factor (Equation 4 of Cannon et al. (2015)) when the denominator is below inv_mod_th.
#' @param inv_mod_th (Optional) A float indicating the upper threshold of the denominator of the scaling factor (Equation 4 of Cannon et al. (2015)) to truncate the scaling factor. This parameter is defined as default as the pp_threshold parameter, described above.
#'
#' @return A column vector of monthly or annual modeled data (temperature or precipitation) corrected by the QDM method. If monthly frequency is specified, the length of this vector is 12 times the number of observed years [12 x y_mod, 1]. If annual frequency is specified, the length of this vector is equal to the number of observed  years [y_mod, 1].
#' @export
#'
#' @examples QDM(obs,mod,var)
#' @examples QDM(obs,mod,var,frq='A')
#' @examples QDM(obs,mod,var,frq='M')
#' @examples QDM(obs,mod,var,frq='M',pp_threshold=0.1)
#' @examples QDM(obs,mod,var,frq='M',pp_factor=1/1000)
#' @examples QDM(obs,mod,var,frq='M',pp_threshold=0.1,pp_factor=1/1000)
#' @examples QDM(obs,mod,var,frq='M',rel_change_th=4)
#' @examples QDM(obs,mod,var,frq='M',inv_mod_th=1/10)
#' @examples QDM(obs,mod,var,frq='M',rel_change_th=4,inv_mod_th=1/10)
QDM <- function(obs, mod, mult_change, allow_negatives, frq, pp_threshold,
                pp_factor ,rel_change_th, inv_mod_th, day_win, user_pdf,
                pdf_obs, pdf_mod){

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
  
  if (frq=='D' & allow_negatives==FALSE) {
    pp_threshold_mod <- get_pp_threshold_mod(obs, mod, pp_threshold)
    
  } else {
    pp_threshold_mod <- pp_threshold
  }

  # 1) Format inputs and get statistics of the observed and modeled series of
  #    the historical period (formatQM function of the climQMBC package).
  format_list_obs <- formatQM(obs, allow_negatives, frq, pp_threshold, pp_factor)
  y_obs <- format_list_obs[[1]]
  obs_series <- format_list_obs[[2]]
  
  format_list_mod <- formatQM(mod, allow_negatives, frq, pp_threshold_mod, pp_factor)
  y_mod <- format_list_mod[[1]]
  mod_series <- format_list_mod[[2]]
  
  if (frq=='D') {
    
    obs_series_moving <- day_centered_moving_window(obs_series, day_win)
    mod_series_moving <- day_centered_moving_window(mod_series, day_win)
    
    win_series <- projected_backward_moving_window(mod_series_moving, y_obs, frq)
    
    obs_series_moving <- array(aperm(obs_series_moving, c(1,3,2)), c(365,(2*day_win-1)*y_obs))
    modh_series_moving <- array(aperm(mod_series_moving,c(1,3,2))[,1:y_obs,], c(365,(2*day_win-1)*y_obs))
    win_series_moving <- array(aperm(win_series,c(1,4,2,3)), c(365, (2*day_win-1)*y_obs, y_mod-y_obs))
    
    
    if (allow_negatives==FALSE) {
      obs_series_moving <- set_norain_to_nan(obs_series_moving, pp_threshold, pp_factor)
      modh_series_moving <- set_norain_to_nan(modh_series_moving, pp_threshold_mod, pp_factor)
      mod_series[mod_series<pp_threshold_mod] <- NaN
      win_series_moving <- set_norain_to_nan(win_series_moving, pp_threshold_mod, pp_factor)
    }
    
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
    win_series <- projected_backward_moving_window(mod_series, y_obs, frq)
    
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
  

  # 2) Assign a probability distribution function to each month for the
  #    observed and modeled data in the historical period. If annual
  #    frequency is specified, this is applied to the complete historical
  #    period (getDist function of the climQMBC package).
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
  

  # 3) For each projected period:
  ks_fail_win <- 0
  pdf_win <- matrix(0,dim(mod_series)[1],dim(mod_series)[2]-y_obs)
  prob <- matrix(0,dim(mod_series)[1],dim(mod_series)[2]-y_obs)
  for (j in 1:dim(prob)[2]){

    # a) Assign a probability distribution function to each month. If
    #    annual frequency is specified, this is applied to the complete
    #    period (getDist function of the climQMBC package).
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
    
    # b) Apply the cumulative distribution function of the projected
    #    period, evaluated with the statistics of this period, to the last
    #    data of the period (getCDF function of the climQMBC package).
    #    Equation 3 of Cannon et al. (2015).
    if (frq=='D') {
      prob[,j] <- getCDF(pdf_win[,j],matrix(win_series[,day_win,j,y_obs]),mu_win[,j],std_win[,j],skew_win[,j],skewy_win[,j])
      prob[is.nan(matrix(win_series[,day_win,j,y_obs])),j] <- NaN
    } else {
      prob[,j] <- getCDF(pdf_win[,j],matrix(mod_series[,y_obs+j]),mu_win[,j],std_win[,j],skew_win[,j],skewy_win[,j])
    }
  }

  # 4) Apply the inverse cumulative distribution function:
  #    a) Of the observed data, evaluated with the statistics of the observed
  #       data in the historical period, to the probabilities obtained from
  #       3b) (getCDFinv function of the climQMBC package). Equation 5 of
  #       Cannon et al. (2015).
  inv_obs <- getCDFinv(pdf_obs,prob,mu_obs,std_obs,skew_obs,skewy_obs)

  #    b) Of the modeled data, evaluated with the statistics of the observed
  #       data in the historical period, to the probabilities obtained from
  #       3b) (getCDFinv function of the climQMBC package). Equations 4 of
  #       Cannon et al. (2015).
  inv_mod <- getCDFinv(pdf_mod,prob,mu_mod,std_mod,skew_mod,skewy_mod)

  # 5) Get the delta factor or relative change and apply it to the value
  #    obtained in 4b). Equation 4 and 6 of Cannon et al. (2015).
  if (mult_change == 1){
    delta_quantile <- mod_series[,(y_obs+1):dim(mod_series)[2]]/inv_mod
    bool_undefined <- (delta_quantile>rel_change_th) & (inv_mod<inv_mod_th)
    delta_quantile[bool_undefined] <- rel_change_th
    QDM <- inv_obs*delta_quantile
  } else{
    delta_quantile <- mod_series[,(y_obs+1):dim(mod_series)[2]] - inv_mod
    QDM <- inv_obs + delta_quantile
  }

  QDM <- matrix(QDM)

  # 6) Perform QM for the historical period.
  mod_h <- mod[1:length(obs)]
  mod_h <- matrix(mod_h)
  QM_series <- QM(obs, mod_h, allow_negatives,frq,pp_threshold, pp_factor, day_win, user_pdf, pdf_obs, pdf_mod)
  
  QDM_series <- matrix(c(QM_series,QDM))
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
#' This function performs bias correction of modeled series based on observed data by the Unbiased Quantile Mapping (UQM) method, as described by Chadwick et al.  (2021). Correction is performed to monthly or annual precipitation or temperature data in a single location. An independent probability distribution function is assigned to each month and to each projected period based on the Kolmogorov-Smirnov test. If annual frequency is specified, this is applied to the complete period. Each projected period  considers a window whose length is equal to the number of years of the historical period and ends in the analyzed year. Correction of the historical period is performed by the Quantile Mapping (QM) method, as described by Cannon et al. (2015).
#'
#' Chadwick et al. (2022) [under revision]
#'
#' @param obs A column vector of monthly or annual observed data (temperature or precipitation). If monthly frequency is specified, the length of this vector is 12 times the number of observed years [12 x y_obs, 1]. If annual frequency is specified, the length of this vector is equal to the number of observed years [y_obs, 1].
#' @param mod A column vector of monthly or annual modeled data (temperature or precipitation). If monthly frequency is specified, the length of this vector is 12 times the number of observed years [12 x y_mod, 1]. If annual frequency is specified, the length of this vector is equal to the number of observed years [y_mod, 1].
#' @param var A flag that identifies if data are temperature or precipitation. This flag tells the getDist function if it has to discard distribution functions that allow negative numbers, and if the terms in the correction equations are multiplied/divided or added/subtracted. Temperature:   var = 0; Precipitation: var = 1
#' @param frq (Optional) A string specifying if the input is annual or monthly data. If not specified, it is set monthly as default. Monthly:   frq = 'M'; Annual:    frq = 'A'
#' @param pp_threshold (Optional) A float indicating the threshold to consider physically null precipitation values.
#' @param pp_factor (Optional) A float indicating the maximum value of the random values that replace physically null precipitation values.
#'
#' @return A column vector of monthly or annual modeled data (temperature or precipitation) corrected by the UQM method. If monthly frequency is specified, the length of this vector is 12 times the number of observed years [12 x y_mod, 1]. If annual frequency is specified, the length of this vector is equal to the number of observed  years [y_mod, 1].
#' @export
#'
#' @examples UQM(obs,mod,var)
#' @examples UQM(obs,mod,var,frq='A')
#' @examples UQM(obs,mod,var,frq='M')
#' @examples UQM(obs,mod,var,frq='M',pp_threshold=0.1)
#' @examples UQM(obs,mod,var,frq='M',pp_factor=1/1000)
#' @examples UQM(obs,mod,var,frq='M',pp_threshold=0.1,pp_factor=1/1000)

UQM <- function(obs, mod, mult_change, allow_negatives, frq, pp_threshold,
                pp_factor, day_win, user_pdf, pdf_obs, pdf_mod){

  
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
  
  if (frq=='D' & allow_negatives==FALSE) {
    pp_threshold_mod <- get_pp_threshold_mod(obs, mod, pp_threshold)
    
  } else {
    pp_threshold_mod <- pp_threshold
  }
  
  
  # 1) Format inputs and get statistics of the observed and modeled series of
  #    the historical period (formatQM function of the climQMBC package).
  format_list_obs <- formatQM(obs, allow_negatives, frq, pp_threshold, pp_factor)
  y_obs <- format_list_obs[[1]]
  obs_series <- format_list_obs[[2]]
  
  format_list_mod <- formatQM(mod, allow_negatives, frq, pp_threshold_mod, pp_factor)
  y_mod <- format_list_mod[[1]]
  mod_series <- format_list_mod[[2]]
  
  if (frq=='D') {
    
    obs_series_moving <- day_centered_moving_window(obs_series, day_win)
    mod_series_moving <- day_centered_moving_window(mod_series, day_win)
    
    win_series <- projected_backward_moving_window(mod_series_moving, y_obs, frq)
    
    obs_series_moving <- array(aperm(obs_series_moving, c(1,3,2)), c(365,(2*day_win-1)*y_obs))
    modh_series_moving <- array(aperm(mod_series_moving,c(1,3,2))[,1:y_obs,], c(365,(2*day_win-1)*y_obs))
    win_series_moving <- array(aperm(win_series,c(1,4,2,3)), c(365, (2*day_win-1)*y_obs, y_mod-y_obs))
    
    
    if (allow_negatives==FALSE) {
      obs_series_moving <- set_norain_to_nan(obs_series_moving, pp_threshold, pp_factor)
      modh_series_moving <- set_norain_to_nan(modh_series_moving, pp_threshold_mod, pp_factor)
      mod_series[mod_series<pp_threshold_mod] <- NaN
      win_series_moving <- set_norain_to_nan(win_series_moving, pp_threshold_mod, pp_factor)
    }
    
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
    win_series <- projected_backward_moving_window(mod_series, y_obs, frq)
    
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


  # 2) Assign a probability distribution function to each month for the
  #    observed and modeled data in the historical period. If annual
  #    frequency is specified, this is applied to the complete historical
  #    period (getDist function of the climQMBC package).
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

  # 4) For each projected period:
  ks_fail_win <- 0
  pdf_win <- matrix(0,dim(mod_series)[1],dim(mod_series)[2]-y_obs)
  prob <- matrix(0,dim(mod_series)[1],dim(mod_series)[2]-y_obs)
  UQM <- matrix(0,dim(mod_series)[1],dim(mod_series)[2]-y_obs)
  for (j in 1:dim(prob)[2]){
    # a) Assign a probability distribution function to each month. If
    #    annual frequency is specified, this is applied to the complete
    #    period (getDist function of the climQMBC package).
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
    
    # b) Apply the cumulative distribution function of the projected
    #    period, evaluated with the statistics of this period, to the last
    #    data of the period (getCDF function of the climQMBC package).
    #    Equation 3 of Cannon et al. (2015).
    if (frq=='D') {
      prob[,j] <- getCDF(pdf_win[,j],matrix(win_series[,day_win,j,y_obs]),mu_win[,j],std_win[,j],skew_win[,j],skewy_win[,j])
    } else {
      prob[,j] <- getCDF(pdf_win[,j],matrix(mod_series[,y_obs+j]),mu_win[,j],std_win[,j],skew_win[,j],skewy_win[,j])
    }

    UQM[,j] <- getCDFinv(pdf_obs,matrix(prob[,j]),matrix(mu_projected[,j]),matrix(sigma_projected[,j]),matrix(skew_projected[,j]),matrix(skewy_projected[,j]))
  }

  UQM <- matrix(UQM)

  # 6) Perform QM for the historical period
  mod_h <- mod[1:length(obs)]
  mod_h <- matrix(mod_h)
  QM_series <- QM(obs, mod_h, allow_negatives,frq,pp_threshold, pp_factor, day_win, user_pdf, pdf_obs, pdf_mod)
  UQM_series <- matrix(c(QM_series,UQM))
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
#' This function performs bias correction of modeled series based on observed data by the Scaled Distribution Mapping (SDM) method, as described by Switanek et al. (2017). Correction is performed to monthly or annual precipitation or temperature data in a single location. The historical period of the modeled series is bias corrected as a scenario where the complete series is replaced by the bias corrected series. On the other hand, for each year after the historical period the method considers a projected period as a window whose length is equal to the number of years of the historical period and ends in the analyzed year. From the bias corrected series, only the last year is saved and all the others are discarded. This approach aims to build a continuum for the bias corrected series.
#'
#' Switanek, B. M., Troch, P., Castro, L. C., Leuprecht, A., Chang, H. I., Mukherjee, R., and Demaria, M. C. E. (2017) Scaled distribution mapping: A bias correction method that preserves raw climate model projected changes. Hydrology &amp; Earth System Sciences, 21, 2649-2666, https://doi.org/10.5194/hess-21-2649-2017.
#'
#' @param obs A column vector of monthly or annual observed data (temperature or precipitation). If monthly frequency is specified, the length of this vector is 12 times the number of observed years [12 x y_obs, 1]. If annual frequency is specified, the length of this vector is equal to the number of observed years [y_obs, 1].
#' @param mod A column vector of monthly or annual modeled data (temperature or precipitation). If monthly frequency is specified, the length of this vector is 12 times the number of observed years [12 x y_mod, 1]. If annual frequency is specified, the length of this vector is equal to the number of observed years [y_mod, 1].
#' @param var A flag that identifies if data are temperature or precipitation. This flag tells the getDist function if it has to discard distribution functions that allow negative numbers, and if the terms in the correction equations are multiplied/divided or added/subtracted. Temperature:   var = 0; Precipitation: var = 1
#' @param frq (Optional) A string specifying if the input is annual or monthly data. If not specified, it is set monthly as default. Monthly:   frq = 'M'; Annual:    frq = 'A'
#' @param pp_threshold (Optional) A float indicating the threshold to consider physically null precipitation values.
#' @param pp_factor (Optional) A float indicating the maximum value of the random values that replace physically null precipitation values.
#'
#' @return A column vector of monthly or annual modeled data (temperature or precipitation) corrected by the SDM method. If monthly frequency is specified, the length of this vector is 12 times the number of observed years [12 x y_mod, 1]. If annual frequency is specified, the length of this vector is equal to the number of observed  years [y_mod, 1].
#' @export
#'
#' @examples SDM(obs,mod,var)
#' @examples SDM(obs,mod,var,frq='A')
#' @examples SDM(obs,mod,var,frq='M')
#' @examples SDM(obs,mod,var,frq='M',pp_threshold=0.1)
#' @examples SDM(obs,mod,var,frq='M',pp_factor=1/1000)
#' @examples SDM(obs,mod,var,frq='M',pp_threshold=0.1,pp_factor=1/1000)
SDM <- function(obs ,mod, SDM_var, frq, pp_threshold, pp_factor, day_win){

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

  lower_lim <- pp_threshold
  CDF_th <- 10^-3

  # 1) Format inputs and get statistics of the observed and modeled series of
  #    the historical period (formatQM function of the climQMBC package).
  format_list_obs <- formatQM(obs, SDM_var, frq, pp_threshold, pp_factor)
  y_obs <- format_list_obs[[1]]
  obs_series <- format_list_obs[[2]]
  
  format_list_mod <- formatQM(mod, SDM_var, frq, pp_threshold, pp_factor)
  y_mod <- format_list_mod[[1]]
  mod_series <- format_list_mod[[2]]
  
  modh_series <- array(mod_series[,1:y_obs], c(dim(mod_series)[1],y_obs))
  
  
  if (frq=='D') {
    
    obs_series <- day_centered_moving_window(obs_series, day_win)
    mod_series <- day_centered_moving_window(mod_series, day_win)
    
    win_series <- projected_backward_moving_window(mod_series, y_obs, frq)
    
    obs_series <- array(aperm(obs_series, c(1,3,2)), c(365,(2*day_win-1)*y_obs))
    modh_series <- array(aperm(mod_series,c(1,3,2))[,1:y_obs,], c(365,(2*day_win-1)*y_obs))
    win_series <- array(aperm(win_series,c(1,4,2,3)), c(365, (2*day_win-1)*y_obs, y_mod-y_obs))
    
    win_series_ <- abind::abind(modh_series, win_series, along=3)
    
    
  } else {
    day_win <- 1
    win_series <- projected_backward_moving_window(mod_series, y_obs, frq)
    
    win_series_ <- abind::abind(modh_series, aperm(win_series,c(1,3,2)), along=3)
  }

  SDM  <- matrix(0,dim(mod_series)[1],y_mod-y_obs)
  SDM_h  <- matrix(0,dim(obs_series)[1],y_obs)
  for (m in 1:dim(mod_series)[1]){
    # 2) Historical period:

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

    # 3) Projected periods:
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

  SDM <- matrix(SDM)
  SDM_h <- matrix(SDM_h)
  SDM_series <- matrix(c(SDM_h,SDM))
  if (SDM_var == 1){
    SDM_series[SDM_series<pp_threshold] <- 0
  }

  return(SDM_series)
}
