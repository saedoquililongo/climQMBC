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
QM <- function(obs,mod,var,frq,pp_threshold,pp_factor){

  if(missing(pp_threshold)) {
    pp_threshold <- 1
  }

  if(missing(pp_factor)) {
    pp_factor <- 1/100
  }

  # 0) Check if annually or monthly data is specified.
  if(missing(frq)) {
    frq <- 'M'
  }

  # 1) Format inputs and get statistics of the observed and modeled series of
  #    the historical period (formatQM function of the climQMBC package).
  format_list <- formatQM(obs,mod,frq,pp_threshold,pp_factor)

  y_obs       <- format_list[[1]]
  obs_series  <- format_list[[2]]
  mod_series  <- format_list[[3]]
  mu_obs      <- format_list[[4]]
  std_obs     <- format_list[[5]]
  skew_obs    <- format_list[[6]]
  skewy_obs   <- format_list[[7]]
  mu_mod      <- format_list[[8]]
  std_mod     <- format_list[[9]]
  skew_mod    <- format_list[[10]]
  skewy_mod   <- format_list[[11]]

  # 2) Assign a probability distribution function to each month for the
  #    observed and modeled data in the historical period. If annual
  #    frequency is specified, this is applied to the complete historical
  #    period (getDist function of the climQMBC package).
  PDF_obs <- getDist(obs_series,mu_obs,std_obs,skew_obs,skewy_obs,var)
  PDF_mod <- getDist(matrix(mod_series[,1:y_obs],nrow=dim(mod_series)[1]),mu_mod,std_mod,skew_mod,skewy_mod,var)

  # 3) Apply the cumulative distribution function of the modeled data,
  #    evaluated with the statistics of the modeled data in the historical
  #    period, to the modeled data (getCDF function of the climQMBC package).
  #    Equation 1 of Cannon et al. (2015).
  Taot <- getCDF(PDF_mod,mod_series,mu_mod,std_mod,skew_mod,skewy_mod)

  # 4) Apply the inverse cumulative distribution function of the observed
  #    data, evaluated with the statistics of the observed data in the
  #    historical period, to the probabilities obtained from 3) (getCDFinv
  #    function of the climQMBC package). Equation 1 of Cannon et al. (2015).
  QM_series <- getCDFinv(PDF_obs,Taot,mu_obs,std_obs,skew_obs,skewy_obs)
  QM_series <- matrix(QM_series)
  if (var == 1){
    QM_series[QM_series<pp_threshold] <- 0
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
DQM <- function(obs,mod,var,frq,pp_threshold,pp_factor){

  if(missing(pp_threshold)) {
    pp_threshold <- 1
  }

  if(missing(pp_factor)) {
    pp_factor <- 1/100
  }

  # 0) Check if annually or monthly data is specified.
  if(missing(frq)) {
    frq <- 'M'
  }

  # 1) Format inputs and get statistics of the observed and modeled series of
  #    the historical period (formatQM function of the climQMBC package).
  format_list <- formatQM(obs,mod,frq,pp_threshold,pp_factor)

  y_obs       <- format_list[[1]]
  obs_series  <- format_list[[2]]
  mod_series  <- format_list[[3]]
  mu_obs      <- format_list[[4]]
  std_obs     <- format_list[[5]]
  skew_obs    <- format_list[[6]]
  skewy_obs   <- format_list[[7]]
  mu_mod      <- format_list[[8]]
  std_mod     <- format_list[[9]]
  skew_mod    <- format_list[[10]]
  skewy_mod   <- format_list[[11]]

  # 2) Assign a probability distribution function to each month for the
  #    observed and modeled data in the historical period. If annual
  #    frequency is specified, this is applied to the complete historical
  #    period (getDist function of the climQMBC package).
  PDF_obs <- getDist(obs_series,mu_obs,std_obs,skew_obs,skewy_obs,var)
  PDF_mod <- getDist(matrix(mod_series[,1:y_obs],nrow=dim(mod_series)[1]),mu_mod,std_mod,skew_mod,skewy_mod,var)

  # 3) Extract the long-term trend from the modeled data:
  #    a) Get the monthly mean of the historical period. If annually
  #       frequency is specified, this is applied to the complete period).
  xbarmh <- matrix(apply(matrix(mod_series[,1:y_obs],nrow=dim(mod_series)[1]),1,mean,na.rm=TRUE))

  #    b) Get the monthly mean of each projected period. If annually
  #       frequency is specified, this is applied to the complete period).
  y_mod <- dim(mod_series)[2]
  xbarmp <- matrix(0,dim(mod_series)[1],(y_mod-y_obs))

  for (j in 1:(y_mod-y_obs)){
    xbarmp[,j] <- apply(matrix(mod_series[,(j+1):(y_obs+j)],nrow=dim(mod_series)[1]),1,mean)
  }

  # 4)Compute the linear scaled values (value in square brackets in equation
  #                                     2 of Cannon et al. (2015)).
  LS <- matrix(0,dim(xbarmp)[1],dim(xbarmp)[2])
  if (var == 1){   # Precipitation
    for (m in 1:dim(mod_series)[1]){
      LS[m,] <- (mod_series[m,(y_obs+1):y_mod]*xbarmh[m])/ xbarmp[m,]
    }
  } else {          # Temperature
    for (m in 1:dim(mod_series)[1]){
      LS[m,] <- (mod_series[m,(y_obs+1):y_mod] + xbarmh[m]) - xbarmp[m,]
    }
  }

  # 5) Apply the cumulative distribution function of the modeled data,
  #    evaluated with the statistics of the modeled period, to the future
  #    modeled data (getCDF function of the climQMBC package). Equation 2 of
  #    Cannon et al. (2015).
  Taot <- getCDF(PDF_mod,LS,mu_mod,std_mod,skew_mod,skewy_mod)

  # 6) Apply the inverse cumulative distribution function of the observed
  #    data, evaluated with the statistics of the observed data in the
  #    historical period, to the probabilities obtained from 5) (getCDFinv
  #    function of the climQMBC package). Equation 2 of Cannon et al. (2015).
  DQM_LS <- getCDFinv(PDF_obs,Taot,mu_obs,std_obs,skew_obs,skewy_obs)

  # 7) Reimpose the trend to the values obtained in 6). Equation 2 of Cannon
  #    et al. (2015).
  if (var == 1){ # Precipitation
    DQM <- DQM_LS*(xbarmp/ matrix(rep(xbarmh,y_mod-y_obs),ncol=y_mod-y_obs))
  } else{        # Temperature
    DQM <- DQM_LS + (xbarmp - matrix(rep(xbarmh,y_mod-y_obs),ncol=y_mod-y_obs))
  }

  DQM <- matrix(DQM)

  # 8) Perform QM for the historical period.
  mod_h <- mod_series[,1:y_obs]
  mod_h <- matrix(mod_h)
  QM_series <- QM(obs,mod_h,var,frq)
  DQM_series <- c(QM_series,DQM)
  if (var == 1){
    DQM_series[DQM_series<pp_threshold] <- 0
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
QDM <- function(obs,mod,var,frq,pp_threshold,pp_factor,rel_change_th,inv_mod_th){

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

  # 0) Check if annually or monthly data is specified.
  if(missing(frq)) {
    frq <- 'M'
  }

  # 1) Format inputs and get statistics of the observed and modeled series of
  #    the historical period (formatQM function of the climQMBC package).
  format_list <- formatQM(obs,mod,frq,pp_threshold,pp_factor)

  y_obs       <- format_list[[1]]
  obs_series  <- format_list[[2]]
  mod_series  <- format_list[[3]]
  mu_obs      <- format_list[[4]]
  std_obs     <- format_list[[5]]
  skew_obs    <- format_list[[6]]
  skewy_obs   <- format_list[[7]]
  mu_mod      <- format_list[[8]]
  std_mod     <- format_list[[9]]
  skew_mod    <- format_list[[10]]
  skewy_mod   <- format_list[[11]]

  # 2) Assign a probability distribution function to each month for the
  #    observed and modeled data in the historical period. If annual
  #    frequency is specified, this is applied to the complete historical
  #    period (getDist function of the climQMBC package).
  PDF_obs <- getDist(obs_series,mu_obs,std_obs,skew_obs,skewy_obs,var)
  PDF_mod <- getDist(matrix(mod_series[,1:y_obs],nrow=dim(mod_series)[1]),mu_mod,std_mod,skew_mod,skewy_mod,var)

  # 3) For each projected period:
  PDF_win <- matrix(0,dim(mod_series)[1],dim(mod_series)[2]-y_obs)
  Taot <- matrix(0,dim(mod_series)[1],dim(mod_series)[2]-y_obs)
  for (j in 1:dim(Taot)[2]){
    win_series <- matrix(mod_series[,(1+j):(y_obs+j)],nrow=dim(mod_series)[1])

    mux <- apply(win_series, 1, mean, na.rm = TRUE)         # Mean
    sigmax <- apply(win_series, 1, sd, na.rm = TRUE)        # Standard deviation
    skewx <- apply(win_series,1,e1071::skewness,na.rm=TRUE,type=2) # Skewness
    Ln_win <- log(win_series)
    Ln_win[Im(Ln_win)!=0] <- 0
    Ln_win[!is.finite(Ln_win)] <- log(0.01)
    skewy <- apply(Ln_win,1,e1071::skewness,na.rm=TRUE,type=2) # Log-Skewness

    # a) Assign a probability distribution function to each month. If
    #    annual frequency is specified, this is applied to the complete
    #    period (getDist function of the climQMBC package).
    PDF_win[,j] <- getDist(win_series,mux,sigmax,skewx,skewy,var)

    # b) Apply the cumulative distribution function of the projected
    #    period, evaluated with the statistics of this period, to the last
    #    data of the period (getCDF function of the climQMBC package).
    #    Equation 3 of Cannon et al. (2015).
    Taot[,j] = getCDF(PDF_win[,j],matrix(mod_series[,y_obs+j]),mux,sigmax,skewx,skewy)
  }

  # 4) Apply the inverse cumulative distribution function:
  #    a) Of the observed data, evaluated with the statistics of the observed
  #       data in the historical period, to the probabilities obtained from
  #       3b) (getCDFinv function of the climQMBC package). Equation 5 of
  #       Cannon et al. (2015).
  inv_obs <- getCDFinv(PDF_obs,Taot,mu_obs,std_obs,skew_obs,skewy_obs)

  #    b) Of the modeled data, evaluated with the statistics of the observed
  #       data in the historical period, to the probabilities obtained from
  #       3b) (getCDFinv function of the climQMBC package). Equations 4 of
  #       Cannon et al. (2015).
  inv_mod <- getCDFinv(PDF_mod,Taot,mu_mod,std_mod,skew_mod,skewy_mod)

  # 5) Get the delta factor or relative change and apply it to the value
  #    obtained in 4b). Equation 4 and 6 of Cannon et al. (2015).
  if (var == 1){
    DM <- mod_series[,(y_obs+1):dim(mod_series)[2]]/inv_mod
    DM[(DM>rel_change_th) & (inv_mod<inv_mod_th)] <- rel_change_th
    QDM <- inv_obs*DM
  } else{
    DM <- mod_series[,(y_obs+1):dim(mod_series)[2]] - inv_mod
    QDM <- inv_obs + DM
  }

  QDM <- matrix(QDM)

  # 6) Perform QM for the historical period.
  mod_h <- mod_series[,1:y_obs]
  mod_h <- matrix(mod_h)
  QM_series <- QM(obs,mod_h,var,frq)
  QDM_series <- c(QM_series,QDM)
  if (var == 1){
    QDM_series[QDM_series<pp_threshold] <- 0
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
UQM <- function(obs,mod,var,frq,pp_threshold,pp_factor){

  if(missing(pp_threshold)) {
    pp_threshold <- 1
  }

  if(missing(pp_factor)) {
    pp_factor <- 1/100
  }

  # 0) Check if annually or monthly data is specified.
  if(missing(frq)) {
    frq <- 'M'
  }

  # 1) Format inputs and get statistics of the observed and modeled series of
  #    the historical period (formatQM function of the climQMBC package).
  format_list <- formatQM(obs,mod,frq,pp_threshold,pp_factor)

  y_obs       <- format_list[[1]]
  obs_series  <- format_list[[2]]
  mod_series  <- format_list[[3]]
  mu_obs      <- format_list[[4]]
  std_obs     <- format_list[[5]]
  skew_obs    <- format_list[[6]]
  skewy_obs   <- format_list[[7]]
  mu_mod      <- format_list[[8]]
  std_mod     <- format_list[[9]]
  skew_mod    <- format_list[[10]]
  skewy_mod   <- format_list[[11]]

  # 2) Assign a probability distribution function to each month for the
  #    observed and modeled data in the historical period. If annual
  #    frequency is specified, this is applied to the complete historical
  #    period (getDist function of the climQMBC package).
  PDF_obs <- getDist(obs_series,mu_obs,std_obs,skew_obs,skewy_obs,var)

  # 3) For each projected period, get the delta factor (delta) and time
  #    dependent (aster) statistics (mean, standard deviation, skewness, and
  #    log skewness). Equations X to X in Chadwick et al. (2021).
  xbarmt <- matrix(0,dim(mod_series)[1],dim(mod_series)[2]-y_obs)
  xhatmt <- matrix(0,dim(xbarmt)[1],dim(xbarmt)[2])
  skwmt <- matrix(0,dim(xbarmt)[1],dim(xbarmt)[2])
  Lskwmt <- matrix(0,dim(xbarmt)[1],dim(xbarmt)[2])

  Dmu    <- matrix(0,dim(xbarmt)[1],dim(xbarmt)[2])
  Dsigma <- matrix(0,dim(xbarmt)[1],dim(xbarmt)[2])
  Dskw   <- matrix(0,dim(xbarmt)[1],dim(xbarmt)[2])
  DLskw  <- matrix(0,dim(xbarmt)[1],dim(xbarmt)[2])
  muAster <- matrix(0,dim(xbarmt)[1],dim(xbarmt)[2])
  sigmaAster <- matrix(0,dim(xbarmt)[1],dim(xbarmt)[2])
  skwAster  <- matrix(0,dim(xbarmt)[1],dim(xbarmt)[2])
  LskwAster <- matrix(0,dim(xbarmt)[1],dim(xbarmt)[2])

  if (var==1){ # Precipitation
    for (i in 1:dim(xbarmt)[1]){
      for (j in 1:dim(xbarmt)[2]){
        xbarmt[i,j] <- mean(mod_series[i,(j+1):(y_obs+j)],na.rm=TRUE)
        xhatmt[i,j] <- sd(mod_series[i,(j+1):(y_obs+j)],na.rm=TRUE)

        skwmt[i,j]  <- e1071::skewness(mod_series[i,(j+1):(y_obs+j)],na.rm=TRUE,type=2)
        Lnskwmt     <- log(mod_series[i,(j+1):(y_obs+j)])
        Lnskwmt[Im(Lnskwmt)!=0] <- 0
        Lnskwmt[!is.finite(Lnskwmt)] <- log(0.01)
        Lskwmt[i,j] <- e1071::skewness(Lnskwmt,na.rm=TRUE,type=2)

        Dmu[i,j]    <- (xbarmt[i,j]- mu_mod[i])/ mu_mod[i]
        Dsigma[i,j] <- (xhatmt[i,j]- std_mod[i])/ std_mod[i]
        Dskw[i,j]   <- (skwmt[i,j] - skew_mod[i])/ skew_mod[i]
        DLskw[i,j]  <- (Lskwmt[i,j]- skewy_mod[i])/ skewy_mod[i]

        muAster[i,j] <- mu_obs[i]*(1+ Dmu[i,j])
        sigmaAster[i,j] <- std_obs[i]*(1+ Dsigma[i,j])
        skwAster[i,j]  <- skew_obs[i]*(1+ Dskw[i,j])
        LskwAster[i,j] <- skewy_obs[i]*(1+ DLskw[i,j])
      }
    }
  } else { # Temperature
    for (i in 1:dim(xbarmt)[1]){
      for (j in 1:dim(xbarmt)[2]){
        xbarmt[i,j] <- mean(mod_series[i,(j+1):(y_obs+j)],na.rm=TRUE)
        xhatmt[i,j] <- sd(mod_series[i,(j+1):(y_obs+j)],na.rm=TRUE)

        skwmt[i,j]  <- e1071::skewness(mod_series[i,(j+1):(y_obs+j)],na.rm=TRUE,type=2)
        Lnskwmt     <- log(mod_series[i,(j+1):(y_obs+j)])
        Lnskwmt[Im(Lnskwmt)!=0] <- 0
        Lnskwmt[!is.finite(Lnskwmt)] <- log(0.01)
        Lskwmt[i,j] <- e1071::skewness(Lnskwmt,na.rm=TRUE,type=2)

        Dmu[i,j]   = xbarmt[i,j] - mu_mod[i]
        Dsigma[i,j] = xhatmt[i,j] - std_mod[i]
        Dskw[i,j]   = skwmt[i,j] - skew_mod[i]
        DLskw[i,j] = Lskwmt[i,j]- skewy_mod[i]

        muAster[i,j]= mu_obs[i]+(Dmu[i,j])
        sigmaAster[i,j]= std_obs[i]+(Dsigma[i,j])
        skwAster[i,j]  = skew_obs[i]+(Dskw[i,j])
        LskwAster[i,j] = skewy_obs[i]+(DLskw[i,j])
      }
    }
  }

  # 4) For each projected period:
  PDF_win <- matrix(0,dim(mod_series)[1],dim(mod_series)[2]-y_obs)
  Taot <- matrix(0,dim(mod_series)[1],dim(mod_series)[2]-y_obs)
  UQM <- matrix(0,dim(mod_series)[1],dim(mod_series)[2]-y_obs)
  for (j in 1:dim(Taot)[2]){
    win_series <- matrix(mod_series[,(1+j):(y_obs+j)],nrow=dim(mod_series)[1])

    mux <- apply(win_series, 1, mean, na.rm = TRUE)         # Mean
    sigmax <- apply(win_series, 1, sd, na.rm = TRUE)        # Standard deviation
    skewx <- apply(win_series,1,e1071::skewness,na.rm=TRUE,type=2) # Skewness
    Ln_win <- log(win_series)
    Ln_win[Im(Ln_win)!=0] <- 0
    Ln_win[!is.finite(Ln_win)] <- log(0.01)
    skewy <- apply(Ln_win,1,e1071::skewness,na.rm=TRUE,type=2) # Log-Skewness

    # a) Assign a probability distribution function to each month. If
    #    annual frequency is specified, this is applied to the complete
    #    period (getDist function of the climQMBC package).
    PDF_win[,j] <- getDist(win_series,mux,sigmax,skewx,skewy,var)

    # b) Apply the cumulative distribution function of the projected
    #    period, evaluated with the statistics of this period, to the last
    #    data of the period (getCDF function of the climQMBC package).
    #    Equation X of Chadwick et al. (2021).
    Taot[,j] <- getCDF(PDF_win[,j],matrix(mod_series[,y_obs+j]),mux,sigmax,skewx,skewy)
  }

  # 5) Apply the inverse cumulative distribution function of the observed
  #    data, evaluated with the time dependent statistics, to the values
  #    obtained in 4b) (getCDFinv function of the climQMBC package). Equation
  #    X of Chadwick et al. (2021).
  for (yr in 1:dim(Taot)[2]){
    UQM[,yr] <- getCDFinv(PDF_obs,matrix(Taot[,yr]),matrix(muAster[,yr]),matrix(sigmaAster[,yr]),matrix(skwAster[,yr]),matrix(LskwAster[,yr]))
  }

  UQM <- matrix(UQM)

  # 6) Perform QM for the historical period
  mod_h <- mod_series[,1:y_obs]
  mod_h <- matrix(mod_h)
  QM_series <- QM(obs,mod_h,var,frq)
  UQM_series <- c(QM_series,UQM)
  if (var == 1){
    UQM_series[UQM_series<pp_threshold] <- 0
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
SDM <- function(obs,mod,var,frq,pp_threshold,pp_factor){

  if(missing(pp_threshold)) {
    pp_threshold <- 1
  }

  if(missing(pp_factor)) {
    pp_factor <- 1/100
  }

  lower_lim <- pp_threshold
  CDF_th <- 10^-3

  # 0) Check if annually or monthly data is specified.
  if(missing(frq)) {
    frq <- 'M'
  }

  # 1) Format inputs and get statistics of the observed and modeled series of
  #    the historical period (formatQM function of the climQMBC package).
  format_list <- formatQM(obs,mod,frq,pp_threshold,pp_factor)

  y_obs       <- format_list[[1]]
  obs_series  <- format_list[[2]]
  mod_series  <- format_list[[3]]

  SDM  <- matrix(0,dim(mod_series)[1],dim(mod_series)[2]-y_obs)
  SDM_h  <- matrix(0,dim(obs_series)[1],y_obs)
  for (m in 1:dim(mod_series)[1]){
    # 2) Historical period:

    # a) [Switanek et al. (2017), step 1)]
    #     For temperature, get the detrended modeled and observed series in
    #     the historical period.
    #     For precipitation, get the rainday values and its frequency for
    #     the modeled and observed series in the historical period.
    if (var == 0){
      D_obs <- pracma::detrend(obs_series[m,])
      D_mod <- pracma::detrend(mod_series[m,1:y_obs])

      mu_obs <- mean(obs_series[m,])
      mu_mod <- mean(mod_series[m,1:y_obs])
    } else{
      D_obs <- sort(obs_series[m,obs_series[m,]>lower_lim])
      D_mod <- sort(mod_series[m,1:y_obs][mod_series[m,1:y_obs]>lower_lim])

      freq_obs <- length(D_obs)/dim(obs_series)[2]
      freq_mod <- length(D_mod)/length(mod_series[m,1:y_obs])

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
    if (var == 0){
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
      fit_obs <- fitdistrplus::fitdist(D_obs, distr = "gamma", method = "mle")
      fit_mod <- fitdistrplus::fitdist(D_mod, distr = "gamma", method = "mle")
      CDF_obs <- pgamma(D_obs,shape=fit_obs[[1]][1],rate=fit_obs[[1]][2])
      CDF_mod <- pgamma(D_mod,shape=fit_mod[[1]][1],rate=fit_mod[[1]][2])
      CDF_obs[CDF_obs>1-CDF_th] <- 1-CDF_th
      CDF_mod[CDF_mod>1-CDF_th] <- 1-CDF_th
    }

    # 3) Projected periods:
    for (j in 0:(dim(mod_series)[2]-y_obs)){
      # c) Initialize correction array.
      corr_temp <- matrix(0,y_obs,1)

      # d) Define projected window.
      win_series <- mod_series[m,(1+j):(y_obs+j)]

      # e) [Switanek et al. (2017), step 1)]
      #    For temperature, get the detrended series of the projected
      #    period.
      #    For precipitation, get the rainday values, its frequency, and
      #    expected raindays for the projected period.
      #    Get the index of the sorted detrended or rainday values.
      if (var == 0){
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
      if (var == 0){
        mu_winD <- mean(D_win,na.rm=TRUE)
        sigma_winD <- sd(D_win,na.rm=TRUE)

        diff_win <- win_series - D_win

        CDF_win <- pnorm(sort(D_win),mu_winD,sigma_winD)
        CDF_win[CDF_win<CDF_th] <- CDF_th
        CDF_win[CDF_win>1-CDF_th] <- 1-CDF_th
      } else {
        fit_win <- fitdistrplus::fitdist(D_win, distr = "gamma", method = "mle")
        CDF_win <- pgamma(D_win,shape=fit_win[[1]][1],rate=fit_win[[1]][2])
        CDF_win[CDF_win>1-CDF_th] <- 1-CDF_th
      }

      # g) [Switanek et al. (2017), step 3)]
      #    Get the scaling between the model projected period and
      #    historical period distributions.
      if (var == 0){
        SF <- (sigma_obsD/sigma_modD)*(qnorm(CDF_win,mu_winD,sigma_winD) - qnorm(CDF_win,mu_modD,sigma_modD))
      } else {
        SF <- qgamma(CDF_win,shape=fit_win[[1]][1],rate=fit_win[[1]][2])/qgamma(CDF_win,shape=fit_mod[[1]][1],rate=fit_mod[[1]][2])
      }

      # h) Interpolate observed and modeled CDF of the historical period
      #    to the length of the projected period.
      obs_cdf_intpol <- pracma::interp1(pracma::linspace(1,length(D_obs),length(D_obs)),CDF_obs,pracma::linspace(1,length(D_obs),length(D_win)))
      mod_cdf_intpol <- pracma::interp1(pracma::linspace(1,length(D_mod),length(D_mod)),CDF_mod,pracma::linspace(1,length(D_mod),length(D_win)))

      # i) [Switanek et al. (2017), steps 4 and 5)]
      #    Get recurrence intervals and its scaled version for the
      #    observed, historical modeled and projected period modeled
      #    CDFs.
      if (var == 0){
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
      if (var == 0){
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
      corr_temp[win_argsort[(length(win_argsort)-exp_D+1):length(win_argsort)]] <- matrix(xvals)
      if (var == 0){
        corr_temp <- corr_temp + diff_win - mu_win
      }

      # l) If the projected period is the historical period (period 0,
      #    j=0) save the complete bias corrected series.
      #    If the projected period is not the historical period (j>0),
      #    save the value of the last year.
      if (j == 0){
        SDM_h[m,] <- corr_temp
      } else {
        SDM[m,j] <- corr_temp[length(corr_temp)]
      }
    }
  }

  SDM <- matrix(SDM)
  SDM_h <- matrix(SDM_h)
  SDM_series <- c(SDM_h,SDM)
  if (var == 1){
    SDM_series[SDM_series<pp_threshold] <- 0
  }

  return(SDM_series)
}
