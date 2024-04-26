get_pp_threshold_mod <- function(obs, mod, pp_threshold){
  
  obs_rainday_hist <- sum(obs>pp_threshold)
  mod_sort_descending <- sort(mod[1:length(obs)], decreasing=TRUE)
  days_kept <- min(obs_rainday_hist, length(obs))
  
  pp_threshold_mod <- mod_sort_descending[days_kept+1]
  
  return(pp_threshold_mod)
}


#' formatQM
#'
#'This function formats the inputs and gets basic statistics for the different Quantile Mapping (QM, DQM, QDM, UQM and SDM) methods available in the climQMBC package. If monthly data is specified, the input series will be reshaped to a matrix of 12 rows and several columns equal to the number of years of each series. If annual data is specified, the input is reshaped to a row vector with same entries as the input series. For precipitation, physically null values (values below pp_threshold) are replaced by random positive values below pp_factor.
#'
#' @param obs A column vector of monthly or annual observed data (temperature or precipitation). If monthly frequency is specified, the length of this vector is 12 times the number of observed years [12 x y_obs, 1]. If annual frequency is specified, the length of this vector is equal to the number of observed years [y_obs, 1].
#' @param mod A column vector of monthly or annual modeled data (temperature or precipitation). If monthly frequency is specified, the length of this vector is 12 times the number of observed years [12 x y_mod, 1]. If annual frequency is specified, the length of this vector is equal to the number of observed years [y_mod, 1].
#' @param var A flag that identifies if data are temperature or precipitation. This flag tells the getDist function if it has to discard distribution functions that allow negative numbers, and if the terms in the correction equations are multiplied/divided or added/subtracted. Temperature:   var = 0; Precipitation: var = 1
#' @param frq A string specifying if the input is annual or monthly data. If not specified, it is set monthly as default. Monthly:   frq = 'M'; Annual:    frq = 'A'
#' @param pp_threshold A float indicating the threshold to consider physically null precipitation values.
#' @param pp_factor A float indicating the maximum value of the random values that replace physically null precipitation values.
#'
#' @return y_obs:          Number of observed years.
#' @return obs_series:     A column vector of monthly or annual observed data (temperature or precipitation). If monthly frequency is specified, the length of this vector is 12 times the number of observed years [12, y_obs]. If annual frequency is specified, the length of this vector is equal to the number of observed years [1, y_obs].
#' @return mod_series:     A column vector of monthly or annual modeled data (temperature or precipitation). If monthly frequency is specified, the length of this vector is 12 times the number of observed years [12, y_mod]. If annual frequency is specified, the length of this vector is equal to the number of observed years [1, y_mod].
#' @return mu_obs:         If monthly frequency is specified, a column vector of monthly mean of observed data [12,1]. If annual frequency is specified, the mean of the observed data (float).
#' @return mu_mod:         If monthly frequency is specified, a column vector of monthly mean of modeled data of the historical period [12,1]. If annual frequency is specified, the mean of the modeled data of the historical period(float).
#' @return sigma_obs:      If monthly frequency is specified, a column vector of monthly standard deviation of observed data [12,1]. If annual frequency is specified, the standard deviation of the observed data (float).
#' @return sigma_mod:      If monthly frequency is specified, a column vector of monthly standard deviation of modeled data of the historical period [12,1]. If annual frequency is specified, the standard deviation of the modeled data of the historical period(float).
#' @return skew_obs:       If monthly frequency is specified, a column vector of monthly skewness of observed data [12,1]. If annual frequency is specified, the skewness of the observed data (float).
#' @return skew_mod:       If monthly frequency is specified, a column vector of monthly skewness of modeled data of the historical period [12,1]. If annual frequency is specified, the skewness of the modeled data of the historical period (float).
#' @return skewy_obs:      If monthly frequency is specified, a column vector of monthly skewness of the logarithm of observed data [12,1]. If annual frequency is specified, the skewness of the logarithm of the observed data (float).
#' @return skewy_mod:      If monthly frequency is specified, a column vector of monthly skewness of the logarithm of modeled data of the historical period [12,1]. If annual frequency is specified, the skewness of the logarithm of the modeled data of the historical period(float).
#' @export
#'
#' @examples formatQM(series_, allow_negatives, frq, pp_threshold, pp_factor)
formatQM <- function(series_, allow_negatives, frq, pp_threshold, pp_factor){
  
  series <- series_

  # 0) Check if annually or monthly data is specified.
  if (frq=='D') {
    I <- 365
  } else if (frq=='M') {
    I <- 12
  } else if (frq=='A') {
    I <- 1
  } else {
    I <- 1
  }

  # 1) If variable is precipitation, replace low values with random values.
  if (allow_negatives==0) {
    bool_low <- series<pp_threshold
    series[bool_low] <- runif(sum(bool_low))*pp_factor*pp_threshold
  }

  # 2) Get number of years of the observed period.
  years <- length(series)/I

  # 3) If monthly data is specified, reshape the input series to a matrix of
  #    12 rows and several columns equal to the number of years of each
  #    series. If annually data is specified, reshape the input to a row
  #    vector with same entries as the input series.
  series <- matrix(series, nrow=I, ncol=length(series)/I)

  return(list(years, series))
}

#' getStats
#'
#'This function computes the mean, standard deviation, skewness and skewness of the logarithmic values, for each sub-period within the year, according to the frequency initially defined.
#'
#' @param series A matrix of monthly or annual data. If monthly frequency is specified, the number of rows of the matrix is 12 and the number of columns is the number of years [12, years]. If annual frequency is specified, the number of rows of the matrix is 1 and the number of columns is the number of years [1, years].
#'
#' @return mu:         If monthly frequency is specified, a column vector of monthly mean of observed data [12,1]. If annual frequency is specified, the mean of the observed data (float).
#' @return sigma:      If monthly frequency is specified, a column vector of monthly standard deviation of observed data [12,1]. If annual frequency is specified, the standard deviation of the observed data (float).
#' @return skew:       If monthly frequency is specified, a column vector of monthly skewness of observed data [12,1]. If annual frequency is specified, the skewness of the observed data (float).
#' @return skewy:     If monthly frequency is specified, a column vector of monthly skewness of observed data [12,1]. If annual frequency is specified, the skewness of the observed data (float).
#' @export
#'
#' @examples getStats(series)
getStats <- function(series, frq){

  if (frq=='D') {
    if (length(dim(series))==3) {
      dim_stats <- c(1,3)
    } else {
      dim_stats <- 1
    }
  } else {
    if (length(dim(series))==3) {
      dim_stats <- c(1,2)
    } else {
      dim_stats <- 1
    }
  }

  
  mu <- apply(series, dim_stats, mean, na.rm = TRUE)         # Mean
  sigma <- apply(series, dim_stats, sd, na.rm = TRUE)        # Standard deviation
  skew <- apply(series,dim_stats,e1071::skewness,na.rm=TRUE,type=1) # Skewness
  series_log <- log(series)
  series_log[!is.finite(series_log)] <- log(runif(sum(!is.finite(series_log)))*0.01)
  skewy <- apply(series_log,dim_stats,e1071::skewness,na.rm=TRUE,type=1) # Log-Skewness

  return(list(mu,sigma,skew,skewy))
}


day_centered_moving_window <- function(series, day_win){
  
  series_moving <- rbind(series[(dim(series)[1]-day_win+1):dim(series)[1],],rep(1,day_win*2) %x% series,series[1:day_win,])
  series_moving <- array(series_moving, c(dim(series)[1]+1, day_win*2, dim(series)[2]))[1:dim(series)[1],2:(day_win*2),]
  series_moving <- array(series_moving, c(dim(series)[1], day_win*2-1, dim(series)[2]))
  
  return(series_moving)
}



projected_backward_moving_window <- function(series, projected_win, frq){
  
  y_mod <- dim(series)[length(dim(series))]
  
  if (frq=='D') {
    day_win <- as.integer(dim(series)[2]+1)/2
    
    win_series <- abind::abind(array(rep(1,projected_win), c(1,1,projected_win)) %x% series,
                               array(0, c(dim(series)[1],day_win*2-1,projected_win)), along=3)
    win_series <- array(win_series, c(dim(series)[1],day_win*2-1,y_mod+1,projected_win))[,,2:(y_mod-projected_win+1),]
    win_series <- array(win_series, c(dim(series)[1],day_win*2-1,y_mod-projected_win,projected_win))
    
  } else {
    win_series <- cbind(array(rep(series,projected_win), c(dim(eries)[1],dim(mod_series)[2]*projected_win)),
                        array(0, c(dim(mod_series)[1],projected_win)))
    win_series <- array(win_series, c(dim(mod_series)[1],y_mod+1,projected_win))[,2:(y_mod-projected_win+1),]
    win_series <- array(win_series, c(dim(mod_series)[1],y_mod-projected_win,projected_win))
  }
  
  return(win_series)
}


set_norain_to_nan <- function(series_moving, pp_threshold, pp_factor, min_rainday){
  
  if(missing(min_rainday)) {
    min_rainday <- 30
  }
  
  if (length(dim(series_moving))==2) {
    rainday_count <- apply(series_moving>pp_threshold,1,sum)
    for (per in 1:length(rainday_count)){
      bool_low <- series_moving[per,]<pp_threshold
      replace_values_nans <-runif(sum(bool_low))*pp_factor*pp_threshold
      replace_values_nans[max(1,min_rainday-rainday_count[per]):length(replace_values_nans)] <- NaN
      series_moving[per,bool_low] <- replace_values_nans
    }
  } else {
    rainday_count <- apply(series_moving>pp_threshold,c(1,3),sum)
    for (per1 in 1:dim(rainday_count)[1]) {
      for (per2 in 1:dim(rainday_count)[2]) { 
        bool_low <- series_moving[per1,,per2]<pp_threshold
        replace_values_nans <-runif(sum(bool_low))*pp_factor*pp_threshold
        replace_values_nans[max(1,min_rainday-rainday_count[per1,per2]):length(replace_values_nans)] <- NaN
        series_moving[per1,bool_low, per2] <- replace_values_nans
      }
    }
  }
  
  return(series_moving)
}

#' Get probability distribution function for each month of the period
#'
#' This function assigns an independent probability distribution function to each row of the input series by comparing the empirical probability distribution function with seven distributions based on the Kolmogorov-Smirnov (KS) test. If the series consider monthly data, it will have 12 rows and each row will represent a month. For annual data the series will have only one row. Only strictly positive distributions are considered for precipitation and strictly positive distributions are discarded if the series has negative values.
#'
#' The available distributions are: 1) Normal distribution; 2) Log-Normal distribution; 3) Gamma 2 parameters distribution; 4) Gamma 3 parameters distribution (Pearson 3 parameters distribution); 5) Log-Gamma 3 parameters distribution (Log-Pearson 3 parameters distribution); 6) Gumbel distribution; 7) Exponential distribution
#'
#' For precipitation, only 2), 3) and 5) are considered (1, 4, 6, and 7 are discarded). For series with negative values, only 1), 3), 4), 6), and 7) are considered (2, 3 and 5 are discarded).
#'
#' @param series A matrix of monthly or annual data (temperature or precipitation). If the series consider monthly data, it will have 12 rows and each row will represent a month. For annual data the series will have only one row.
#' @param mu A column vector of mean values of the series. [12,1] if the series consider monthly data and [1,1] if the series consider annual data.
#' @param sigma A column vector of standard deviation of the series. [12,1] if the series consider monthly data and [1,1] if the series consider annual data.
#' @param skew A column vector of skewness of the series. [12,1] if the series consider monthly data and [1,1] if the series consider annual data.
#' @param skewy A column vector of skewness of the logarithm of the series. [12,1] if the series consider monthly data and [1,1] if the series consider annual data.
#' @param var A flag that identifies if data are temperature or precipitation. This flag tells the getDist function if it has to discard distribution functions that allow negative numbers. Temperature:   var = 0; Precipitation: var = 1
#'
#' @return PDF: A column vector with an ID for the resulting distribution from the KS test. [12,1] if the series consider monthly data and [1,1] if the series consider annual data. The ID is related to the numeration of the distribution listed in the description of this function. This ID is used in the getCDF and getCDFinv functions of the climQMBC package.
#' @export
#'
#' @examples getDist(series,mu,sigma,skew,skewy,var)
getDist <- function(series, allow_negatives, mu, sigma, skew, skewy){

  # 1) Get the number of years to compute the empirical distribution in step
  #    3) and get the number of rows of the input series.
  n <- dim(series)[1]

  # 2) Initialize column vectors for the statistics needed for the available
  #    probability distribution functions.
  PDF    <- matrix(0,n,1)
  sigmay <- matrix(0,n,1)
  muy    <- matrix(0,n,1)
  A      <- matrix(0,n,1)
  B      <- matrix(0,n,1)
  Alp    <- matrix(0,n,1)
  Bet    <- matrix(0,n,1)
  Gam    <- matrix(0,n,1)
  Alpy   <- matrix(0,n,1)
  Bety   <- matrix(0,n,1)
  Gamy   <- matrix(0,n,1)
  a      <- matrix(0,n,1)
  u      <- matrix(0,n,1)
  
  ks_fail <- 0

  # 3) Perform the Kolmogorov-Smirnov test for each row.
  for (m in 1:n){
    series_sub <- series[m,]
    series_sub <- series_sub[!is.nan(series_sub)]
      
      
    y_series <- length(series_sub)

    # a) Get empirical distribution.
    sortdata <- sort(series_sub)
    probEmp  <- seq(from = 1/(y_series+1), to = y_series/(y_series+1), by = 1/(y_series+1))

    # b) Compare distributions.
    # i) Normal distribution.
    normal   <- pnorm(sortdata,mu[m],sigma[m])
    KSnormal <- max(abs(normal-probEmp))

    # ii) Log Normal distribution.
    if (any(series_sub < 0)){
      KSlognormal <- 1
    }else{
      sigmay[m] <- sqrt(log(1+(sigma[m]/mu[m])^2))
      muy[m]    <- log(mu[m])-(sigmay[m]^2)/2
      lognormal <- plnorm(sortdata,muy[m],sigmay[m])
      KSlognormal = max(abs(lognormal-probEmp))
    }

    # iii) Gamma 2 parameters distribution.
    if (any(series_sub < 0)){
      KSgammaII <- 1
    }else{
      A[m]    <- (sigma[m]^2)/mu[m]
      B[m]    <- (mu[m]/sigma[m])^2
      GammaII <- pgamma(sortdata,shape=B[m],scale=A[m])
      KSgammaII <- max(abs(GammaII-probEmp))
    }

    #% iv) Gamma 3 parameters distribution.
    #     (Pearson 3 parameters distribution)
    Bet[m]   <- (2/skew[m])^2
    Alp[m]   <- sigma[m]/sqrt(Bet[m])
    Gam[m]   <- mu[m]-(Alp[m]*Bet[m])
    GammaIII <- pgamma((sortdata-Gam[m]),shape=Bet[m],scale=Alp[m])
    KSgammaIII <- max(abs(GammaIII-probEmp))

    # v) Log-Gamma 3 parameters distribution.
    #    (Log-Pearson 3 parameters distribution)
    if (any(series_sub < 0)){
      KSLpIII <- 1
    }else{
      sigmay[m] <- sqrt(log(1+(sigma[m]/mu[m])^2))
      muy[m]    <- log(mu[m])-(sigmay[m]^2)/2
      Bety[m] <- (2/skewy[m])^2
      Alpy[m] <- sigmay[m]/sqrt(Bety[m])
      Gamy[m] <- muy[m]-(Alpy[m]*Bety[m])
      Lnsortdata = log(sortdata)
      Lnsortdata[!is.finite(Lnsortdata)] <- log(0.01)
      LpIII <- pgamma((Lnsortdata-Gamy[m]),shape=Bety[m],scale=Alpy[m])
      KSLpIII <- max(abs(LpIII-probEmp))
    }

    # vi) Gumbel distribution.
    Sn    <- pi/sqrt(6)
    yn    <- 0.5772
    a[m] <- Sn/sigma[m]
    u[m] <- mu[m]-(yn/a[m])
    gumbel <- exp(-exp(-a[m]*(sortdata-u[m])))
    KSgumbel <- max(abs(gumbel-probEmp))

    # vii) Exponential distribution.
    gamexp <- mu[m]-sigma[m]
    exponential <- pmax(1-exp(-1/sigma[m]*(sortdata-gamexp)),0)
    KSexponential <- max(abs(exponential-probEmp))
    
    # c) If variable is precipitation, set KS=1 to distributions that allow
    #    negative values (this will discard those distributions).
    if (allow_negatives==0){
      KSnormal <- 1
      KSgammaIII <- 1
      KSgumbel <- 1
      KSexponential <- 1
    }

    # d) The distribution with lower KS value is considered for each month.
    KS_vals <- c(KSnormal,KSlognormal,KSgammaII,KSgammaIII,KSLpIII,KSgumbel,KSexponential)
    bestPDF <- which.min(KS_vals)
    min_KS <- min(KS_vals)
    
    ks_crit <- 1.3581/sqrt(y_series)
    ks_fail <- ks_fail + (sign(min_KS-ks_crit)+1)/2

    PDF[m] <- bestPDF
  }

  return(list(PDF, ks_fail))
}

#' Get probability of a set of values
#'
#' This function evaluates each row of the series in the respective cumulative distribution function assigned by the Kolmogorov-Smirnov (KS) test in the getDist function of the climQMBC package.
#'
#' The available distributions are: 1) Normal distribution; 2) Log-Normal distribution; 3) Gamma 2 parameters distribution; 4) Gamma 3 parameters distribution (Pearson 3 parameters distribution); 5) Log-Gamma 3 parameters distribution (Log-Pearson 3 parameters distribution); 6) Gumbel distribution; 7) Exponential distribution
#'
#' @param PDF A column vector with an ID for the resulting distribution from the KS test. [12,1] if the series consider monthly data and [1,1] if the series consider annual data. The ID is related to the numeration of the distribution listed in the description of this function.
#' @param series A matrix of monthly or annual data (temperature or precipitation). If the series consider monthly data, it will have 12 rows and each row will represent a month. For annual data the series will have only one row.
#' @param mu A column vector of mean values of the series. [12,1] if the series consider monthly data and [1,1] if the series consider annual data.
#' @param sigma A column vector of standard deviation of the series. [12,1] if the series consider monthly data and [1,1] if the series consider annual data.
#' @param skew A column vector of skewness of the series. [12,1] if the series consider monthly data and [1,1] if the series consider annual data.
#' @param skewy A column vector of skewness of the logarithm of the series. [12,1] if the series consider monthly data and [1,1] if the series consider annual data.
#'
#' @return Taot: A column vector with the non-exceedance probability for each row of the input series. [12,1] if the series consider monthly data and [1,1] if the series consider annual data.
#' @export
#'
#' @examples getCDF(PDF,series,mu,sigma,skew,skewy)
getCDF <- function(PDF,series,mu,sigma,skew,skewy){

  # 1) Get the number of rows and years of the series.
  n_m <- dim(series)[1]
  n_y <- dim(series)[2]

  # 2) Compute the cumulative distribution function to the values of each
  #    row, based on the distribution assigned in the getDist function of the
  #    climQMBC package.
  Taot <- matrix(0,n_m,n_y)
  for (m in 1:n_m){
    if (PDF[m] == 1){ # i) Normal distribution.
      Taot[m,] <- pnorm(series[m,],mu[m],sigma[m])

    } else if (PDF[m] == 2){ # ii) Log-Normal distribution.
      sigmay <- sqrt(log(1+(sigma[m]/mu[m])^2))
      muy <- log(mu[m])-(sigmay^2)/2
      Taot[m,] <- plnorm(series[m,],muy,sigmay)

    } else if (PDF[m] == 3){ # iii) Gamma 2 parameters distribution.
      A <- (sigma[m]^2)/mu[m]
      B <- (mu[m]/sigma[m])^2
      Taot[m,] <- pgamma(series[m,],shape=B,scale=A)

    } else if (PDF[m] == 4){ # iv) Gamma 3 parameters distribution.
      Bet  <- (2/skew[m])^2
      Alp  <- sigma[m]/sqrt(Bet)
      Gam  <- mu[m]-(Alp*Bet)
      Taot[m,]  <- pgamma((series[m,]-Gam),shape=Bet,scale=Alp)

    } else if (PDF[m] == 5){ # v) Log-Gamma 3 parameters distribution.
      Bety   <- (2/skewy[m])^2
      sigmay <- sqrt(log(1+(sigma[m]/mu[m])^2))
      Alpy   <- sigmay/sqrt(Bety)
      muy    <- log(mu[m])-(sigmay^2)/2
      Gamy   <- muy - (Alpy*Bety)
      Lnsortdata <- log(series[m,])
      Lnsortdata[Im(Lnsortdata)!=0] <- log(0.01)
      Lnsortdata[!is.finite(Lnsortdata)] <- log(0.01)
      Taot[m,]  <- pgamma((Lnsortdata-Gamy),shape=Bety,scale=Alpy)

    } else if (PDF[m] == 6){ # vi) Gumbel distribution.
      Sn <- pi/sqrt(6)
      yn <- 0.5772
      a <- Sn/sigma[m]
      u <- mu[m]-(yn/a)
      Taot[m,] <- exp(-exp(-a*(series[m,]-u)))

    } else if (PDF[m] == 7){ # vii) Exponential distribution.
      gamexp <- mu[m] - sigma[m]
      Taot[m,] <- pmax(1-exp(-1/sigma[m]*(series[m,]-gamexp)),0)
    }
  }

  # 3) Tail probabilities are set to (0 +  threshold) and (1 - threshold) to
  #    avoid numerical errors.
  th <- 10^-3
  Taot[Taot>1-th] <- 1-th
  Taot[Taot<th] <- th

  return(Taot)
}


#' Get the value associated to a certain probability
#'
#' This function evaluates the probability, Taot, in the respective inverse cumulative distribution function assigned by the Kolmogorov-Smirnov (KS) test in the getDist function of the climQMBC package.
#'
#' The available distributions are: 1) Normal distribution; 2) Log-Normal distribution; 3) Gamma 2 parameters distribution; 4) Gamma 3 parameters distribution (Pearson 3 parameters distribution); 5) Log-Gamma 3 parameters distribution (Log-Pearson 3 parameters distribution); 6) Gumbel distribution; 7) Exponential distribution
#'
#' @param PDF A column vector with an ID for the resulting distribution from the KS test. [12,1] if the series consider monthly data and [1,1] if the series consider annual data. The ID is related to the numeration of the distribution listed in the description of this function.
#' @param Taot A column vector with the non-exceedance probability for each row of the input series. [12,1] if the series consider monthly data and [1,1] if the series consider annual data.
#' @param mu A column vector of mean values of the series. [12,1] if the series consider monthly data and [1,1] if the series consider annual data.
#' @param sigma A column vector of standard deviation of the series. [12,1] if the series consider monthly data and [1,1] if the series consider annual data.
#' @param skew A column vector of skewness of the series. [12,1] if the series consider monthly data and [1,1] if the series consider annual data.
#' @param skewy A column vector of skewness of the logarithm of the series. [12,1] if the series consider monthly data and [1,1] if the series consider annual data.
#'
#' @return  xhat: A column vector with the values obtained when the inverse cumulative distribution function is applied. [12,1] if the series consider monthly data and [1,1] if the series consider annual data.
#' @export
#'
#' @examples getCDFinv(PDF,Taot,mu,sigma,skew,skewy)
getCDFinv <- function(PDF,Taot,mu,sigma,skew,skewy){

  # 1) Get the number of rows and years of the series.
  n_m <- dim(Taot)[1]
  n_y <- dim(Taot)[2]

  # 2) Compute the inverse cumulative distribution function to the values of
  #    each row, based on the distribution assigned in the getDist function
  #    of the climQMBC package.
  xhat <- matrix(0,n_m,n_y)
  for (m in 1:n_m){
    if (PDF[m] == 1){ # i) Normal distribution.
      xhat[m,] <- qnorm(Taot[m,],mu[m],sigma[m])

    } else if (PDF[m] == 2){ # ii) Log-Normal distribution.
      sigmay <- sqrt(log(1+(sigma[m]/mu[m])^2))
      muy <- log(mu[m])-(sigmay^2)/2
      xhat[m,] <- qlnorm(Taot[m,],muy,sigmay)

    } else if (PDF[m] == 3){ # iii) Gamma 2 parameters distribution.
      A <- (sigma[m]^2)/mu[m]
      B <- (mu[m]/sigma[m])^2
      xhat[m,] <- qgamma(Taot[m,],shape=B,scale=A)

    } else if (PDF[m] == 4){ # iv) Gamma 3 parameters distribution.
      Bet  <- (2/skew[m])^2
      Alp  <- sigma[m]/sqrt(Bet)
      Gam  <- mu[m]-(Alp*Bet)
      xhat[m,]  <- qgamma(Taot[m,],shape=Bet,scale=Alp) + Gam

    } else if (PDF[m] == 5){ # v) Log-Gamma 3 parameters distribution.
      Bety   <- (2/skewy[m])^2
      sigmay <- sqrt(log(1+(sigma[m]/mu[m])^2))
      Alpy   <- sigmay/sqrt(Bety)
      muy    <- log(mu[m])-(sigmay^2)/2
      Gamy   <- muy - (Alpy*Bety)
      xhat[m,]  <- exp(qgamma(Taot[m,],shape=Bety,scale=Alpy) + Gamy)

    } else if (PDF[m] == 6){ # vi) Gumbel distribution.
      Sn <- pi/sqrt(6)
      yn <- 0.5772
      a <- Sn/sigma[m]
      u <- mu[m]-(yn/a)
      xhat[m,] <- u-log(-log(Taot[m,]))/a

    } else if (PDF[m] == 7){ # vii) Exponential distribution.
      gamexp <- mu[m] - sigma[m]
      xhat[m,] <- gamexp-(sigma[m]*log(1-Taot[m,]))
    }
  }

  return(xhat)
}
