#' get_pp_threshold_mod
#'
#' This functions gets a no-rain value threshold for the modeled series by mathcing the no-rain days of the modeled series in the historical period with the observed series.
#'
#' @param obs A column vector of daily observed data, without  considering leap  days. The length of the column vector  should by a multiple of 365. [ndata_obs, 1]
#' @param mod A column vector of daily modeled data, without  considering leap  days. The length of the column vector should by a multiple of 365. [ndata_mod, 1]
#' @param pp_threshold A float indicating the threshold to consider no-rain values in the observed data.
#'
#' @return pp_threshold_mod: A float indicating the threshold to consider no-rain values in the modeled data.
#' @export
#'
#' @examples get_pp_threshold_mod(obs, mod, pp_threshold)
get_pp_threshold_mod <- function(obs, mod, pp_threshold){
  
  # Get the rain days of the observed series
  obs_rainday_hist <- sum(obs>pp_threshold)
  
  # Define the pp_threshold of the modeled series by matching the rain days
  # of the observed series
  mod_sort_descending <- sort(mod[1:length(obs)], decreasing=TRUE)
  days_kept <- min(obs_rainday_hist, length(obs))
  
  if (days_kept!=obs_rainday_hist){
    pp_threshold_mod <- mod_sort_descending[days_kept+1]
  } else {
    pp_threshold_mod <- 0
  }

  return(pp_threshold_mod)
}


#' formatQM
#'
#' This function formats a time series from a column vector to a matrix. The number of rows depends on the frequency of the data (1, 12 and 365 for annual, monthyl and daily data, respectively). If negative values are not allowed, determined by allow_negatives=1, no-rain values (values below pp_threshold) are replaced by random positive values below the product pp_threshold*pp_factor.
#'
#' @param series_ A column vector of daily data, without considering leap days. The length of the column vector should by a multiple of 365. [ndata_obs, 1]
#' @param allow_negatives A flag that identifies if data allows negative values and also to replace no-rain values with random small values (Chadwick et al., 2023) to avoid numerical problems with the probability distribution functions. allow_negatives = 1 or True: Allow negatives ; allow_negatives = 0 or False: Do not allow negative
#' @param frq A string specifying if the input frequency is daily, monthly or annual. Daily:     frq = 'D' ; Monthly:   frq = 'M' ; Annual:    frq = 'A'
#' @param pp_threshold A float indicating the threshold to consider no-rain values.
#' @param pp_factor A float which multiplied to pp_threshold indicates the maximum value of the random values that replace no-rain values.
#'
#' @return years: Number of years of the series.
#' @return series_matrix: A matrix of daily, monthly or annual observed data, without considering leap days. For daily, monthly and annual frequency, the number of rows will be 365, 12 and 1, respectively. The number of columns will be equal to the number of years. [1, 12 or 365, years]
#' @export
#'
#' @examples formatQM(series_, allow_negatives, frq, pp_threshold, pp_factor)
formatQM <- function(series_, allow_negatives, frq, pp_threshold, pp_factor){
  
  # This prevents modyfing the original input passed by reference
  ## Must look for a more elegant way to skip this step
  series <- series_

  # Check frequency and define the number of rows
  if (frq=='D') {
    I <- 365
  } else if (frq=='M') {
    I <- 12
  } else if (frq=='A') {
    I <- 1
  } else {
    I <- 1
  }

  # If variable is precipitation, determined by allow_negatives=1, replace low
  # values with random values.
  if (allow_negatives==0) {
    bool_low <- series<pp_threshold
    series[bool_low] <- runif(sum(bool_low))*pp_factor*pp_threshold
  }

  # Get number of years of the series.
  years <- length(series)/I

  # Reshape to get a matrix of shape [I , years]
  series <- matrix(series, nrow=I, ncol=length(series)/I)

  return(list(years, series))
}


#' getStats
#'
#' This function computes the mean, standard deviation, skewness and skewness of the logarithmic values, for each sub-period within the year, according to the frequency initially defined. Statistics are computed along axis 1. NOTE: If input series is a 2D array, outputs will be a column vector. If input series is a 3D array, outputs will be a 2D array.
#'
#' @param series An array of daily, monthly or annual data, without considering leap days. Possible input dimensions: 2D array: [sub-periods, years] or [sub-periods + days window, years] ; 3D array: [sub_periods, projected periods window, years] or [sub-periods + days window, projected periods window, years]
#'
#' @return mu: Mean values of each sub-period (and projected period if input is a 3D array).
#' @return sigma: Standar deviation values of each sub-period (and projected  period if input is a 3D array).
#' @return skew: Skewness values of each sub-period (and projected period if input is a 3D array).
#' @return skewy: Skewness values of the logarithm of theseries of each sub-period (and projected period if input is a 3D array).
#' @export
#'
#' @examples getStats(series, frq)
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

  # Get the mean, standard deviation, skewness and skewness of the
  # logarithmic values of each year sub-period of the series.
  mu <- apply(series, dim_stats, mean, na.rm = TRUE)         # Mean
  sigma <- apply(series, dim_stats, sd, na.rm = TRUE)        # Standard deviation
  skew <- apply(series,dim_stats,e1071::skewness,na.rm=TRUE,type=1) # Skewness
  series_log <- log(series)
  series_log[!is.finite(series_log)] <- log(runif(sum(!is.finite(series_log)))*0.01)
  skewy <- apply(series_log,dim_stats,e1071::skewness,na.rm=TRUE,type=1) # Log-Skewness

  return(list(mu,sigma,skew,skewy))
}


#' day_centered_moving_window
#'
#' This function transform a 2D matrix of dimension [days in year, years] to a 3D array, where the new dimension is a centered moving window for each day of the year.
#'
#' @param series_matrix A matrix of daily data, without considering leap days. The number of rows is 365 and number of columns is the  number of years of the series. [365, years]
#' @param day_win An integer indicating how many days to consider backwards and forward to get the statistics of each calendar day.The length of the window will be (2*win_day-1). For example, day_win=15 -> window of 29.
#' 
#' @return series_moving: A 3D array of daily data, without considering leap days, with a dimension that considers a centered moving window for each day of the year. [365, 2*win-1, years]
#' @export
#'
#' @examples day_centered_moving_window(series, day_win)
day_centered_moving_window <- function(series, day_win){
  
  series_moving <- rbind(series[(dim(series)[1]-day_win+1):dim(series)[1],],rep(1,day_win*2) %x% series,series[1:day_win,])
  series_moving <- array(series_moving, c(dim(series)[1]+1, day_win*2, dim(series)[2]))[1:dim(series)[1],2:(day_win*2),]
  
  # Avoid dropping dimensions of length 1
  series_moving <- array(series_moving, c(dim(series)[1], day_win*2-1, dim(series)[2]))
  
  return(series_moving)
}


#' projected_backward_moving_window
#'
#' This function transform a 2D or 3D array of dimensions [days in year, year] or [days in year, centered window for each day, years] to a 3D or 4D array adding a new dimension to represent a backward moving window for each  projected period.
#'
#' @param series A 2D or 3D array of daily, monthly or annual observed data, without considering leap days. For daily, monthly and annual frequency, the number of rows will be 365, 12 and 1, respectively. If daily data, it is a 3D array that has a  dimension with a centered moving window for each day of the year. [1, 12 or 365, years] or [365, day centered window, years]
#' @param projected_win An integer the length of the backward moving window.
#' @param frq A string specifying if the input frequency is daily, monthly or annual. Daily:     frq = 'D' ; Monthly:   frq = 'M' ; Annual:    frq = 'A'
#' 
#' @return win_series: An array of daily, monthly or annual future data, without considering leap days, with a dimension associated to a backward moving window for each projected period. Possible dimensions: 3D array: [sub-periods, projected periods window, years] ; 4D array: [sub-periods, days window, projected periods window, years]
#' @export
#'
#' @examples projected_backward_moving_window(series, projected_win, frq)
projected_backward_moving_window <- function(series, projected_win, frq){
  
  y_mod <- dim(series)[length(dim(series))]
  
  if (frq=='D') {
    day_win <- as.integer(dim(series)[2]+1)/2
    
    win_series <- abind::abind(array(rep(1,projected_win), c(1,1,projected_win)) %x% series,
                               array(0, c(dim(series)[1],day_win*2-1,projected_win)), along=3)
    win_series <- array(win_series, c(dim(series)[1],day_win*2-1,y_mod+1,projected_win))[,,2:(y_mod-projected_win+1),]
    
    # Avoid dropping dimensions of length 1
    win_series <- array(win_series, c(dim(series)[1],day_win*2-1,y_mod-projected_win,projected_win))
    
  } else {
    win_series <- cbind(array(rep(series,projected_win), c(dim(series)[1],dim(series)[2]*projected_win)),
                        array(0, c(dim(series)[1],projected_win)))
    win_series <- array(win_series, c(dim(series)[1],y_mod+1,projected_win))[,2:(y_mod-projected_win+1),]
    
    # Avoid dropping dimensions of length 1
    win_series <- array(win_series, c(dim(series)[1],y_mod-projected_win,projected_win))
  }
  
  return(win_series)
}


#' projected_backward_moving_window
#'
#' This function replace no-rain values with nans, leaving a minimum amout of values (min_rainday) to fit a distribution. If fewer than min_rainday values are above pp_threshold, random small values will be added.
#'
#' @param series_moving A 2D or 3D array of daily, monthly or annual observed data, without considering leap days. For daily, monthly and annual frequency, the number of rows will be 365, 12 and 1, respectively. If daily data, it is a 3D array that has a  dimension with a centered moving window for each day of the year. [1, 12 or 365, years] or [365, day centered window, years]
#' @param pp_factor A float which multiplied to pp_threshold indicates the maximum value of the random values that replace no-rain values.
#' @param pp_threshold A float indicating the threshold to consider no-rain values.
#' @param min_rainday (Optional) Minimum amount of values to keep for each sub-period and projected period, to ensure a minimum amount to fit a  distribution. Default is 30
#' 
#' @return series_moving: The input, but with no-rain values replaced with nans.
#' @export
#'
#' @examples set_norain_to_nan(series_moving, pp_threshold, pp_factor)
#' @examples set_norain_to_nan(series_moving, pp_threshold, pp_factor, min_rainday=20)
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


#' getDist
#' 
#' This function assigns an independent probability distribution function to each row of the input series by comparing the empirical probability distribution function with seven distributions based on the Kolmogorov-Smirnov (KS) test.
#'
#' The available distributions are: 1) Normal ; 2) Log-Normal ; 3) Gamma 2 parameters ; 4) Gamma 3 parameters ; (Pearson 3 parameters) ; 5) Log-Gamma 3 parameters ; 6) Gumbel ; 7) Exponential
#'
#' For allow_negatives=0, only 2), 3) and 5) are considered (1, 4, 6, and 7 are discarded). For series with negative values, only 1), 3), 4), 6), and 7) are considered (2, 3 and 5 are discarded).
#'
#' @param series An array of daily, monthly or annual data, without considering leap days. Possible input dimensions: 2D array: [sub-periods, years] or [sub-periods + days window, years]
#' @param allow_negatives A flag that identifies if data allows negative values and also to replace no-rain values with random small values (Chadwick et al., 2023) to avoid numerical problems with the probability distribution functions. allow_negatives = 1 or True: Allow negatives ; allow_negatives = 0 or False: Do not allow negative
#' @param mu A vector with mean values of each sub-period.
#' @param sigma A vector with standard deviation values of each sub-period.
#' @param skew A vector with skewness values of each sub-period.
#' @param skewy A vector with skewness values of the logarithm of the series of each sub-period.
#' 
#' @return pdf: A vector with an ID for the resulting distribution from the KS test. The ID is related to  the numeration of the distribution listed in the description of this function. This ID is used in the getCDF and getCDFinv  functions of the climQMBC package.
#' @export
#'
#' @examples getDist(series, allow_negatives, mu, sigma, skew, skewy)
getDist <- function(series, allow_negatives, mu, sigma, skew, skewy){

  # Get the number of rows of the input series.
  nrows <- dim(series)[1]

  # Initialize column vectors for the statistics needed for the available
  # probability distribution functions.
  pdf    <- matrix(0,nrows,1)
  
  # Initialize a counter of failures of the KS test
  ks_fail <- 0

  # Perform the Kolmogorov-Smirnov test for sub-period or row.
  for (sp in 1:nrows){
    series_sub <- series[sp,]
    
    # Get the number of years with valid data (for allow_negatives=0, the 
    # series might have nan values to ignore them)
    series_sub <- series_sub[!is.nan(series_sub)]
    y_series <- length(series_sub)

    # a) Get empirical distribution.
    sortdata <- sort(series_sub)
    probEmp  <- seq(from = 1/(y_series+1), to = y_series/(y_series+1), by = 1/(y_series+1))

    # b) Compare distributions.
    # i) Normal distribution.
    normal   <- pnorm(sortdata,mu[sp],sigma[sp])
    KSnormal <- max(abs(normal-probEmp))

    # ii) Log Normal distribution.
    if (any(series_sub < 0)){
      KSlognormal <- 1
    }else{
      sigmay <- sqrt(log(1+(sigma[sp]/mu[sp])^2))
      muy <- log(mu[sp])-(sigmay^2)/2
      lognormal <- plnorm(sortdata,muy,sigmay)
      KSlognormal = max(abs(lognormal-probEmp))
    }

    # iii) Gamma 2 parameters distribution.
    if (any(series_sub < 0)){
      KSgammaII <- 1
    }else{
      A <- (sigma[sp]^2)/mu[sp]
      B <- (mu[sp]/sigma[sp])^2
      GammaII <- pgamma(sortdata,shape=B,scale=A)
      KSgammaII <- max(abs(GammaII-probEmp))
    }

    #% iv) Gamma 3 parameters distribution.
    #     (Pearson 3 parameters distribution)
    Bet <- (2/skew[sp])^2
    Alp <- sigma[sp]/sqrt(Bet)
    Gam <- mu[sp]-(Alp*Bet)
    GammaIII <- pgamma((sortdata-Gam),shape=Bet,scale=Alp)
    KSgammaIII <- max(abs(GammaIII-probEmp))

    # v) Log-Gamma 3 parameters distribution.
    #    (Log-Pearson 3 parameters distribution)
    if (any(series_sub < 0)){
      KSLpIII <- 1
    }else{
      sigmay <- sqrt(log(1+(sigma[sp]/mu[sp])^2))
      muy <- log(mu[sp])-(sigmay^2)/2
      Bety <- (2/skewy[sp])^2
      Alpy <- sigmay/sqrt(Bety)
      Gamy <- muy-(Alpy*Bety)
      Lnsortdata = log(sortdata)
      Lnsortdata[!is.finite(Lnsortdata)] <- log(runif(sum(!is.finite(Lnsortdata)))*0.01)
      LpIII <- pgamma((Lnsortdata-Gamy),shape=Bety,scale=Alpy)
      KSLpIII <- max(abs(LpIII-probEmp))
    }

    # vi) Gumbel distribution.
    Sn <- pi/sqrt(6)
    yn <- 0.5772
    a <- Sn/sigma[sp]
    u <- mu[sp]-(yn/a)
    gumbel <- exp(-exp(-a*(sortdata-u)))
    KSgumbel <- max(abs(gumbel-probEmp))

    # vii) Exponential distribution.
    gamexp <- mu[sp]-sigma[sp]
    exponential <- pmax(1-exp(-1/sigma[sp]*(sortdata-gamexp)),0)
    KSexponential <- max(abs(exponential-probEmp))
    
    # KSLpIII <- 1
    
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

    pdf[sp] <- bestPDF
  }
  
  return(list(pdf, ks_fail))
}


#' getCDF
#' 
#' This function evaluates each row of the series in the respective cumulative distribution function assigned by the Kolmogorov-Smirnov (KS) test in the getDist function of the climQMBC package.
#'
#' The available distributions are: 1) Normal distribution; 2) Log-Normal distribution; 3) Gamma 2 parameters distribution; 4) Gamma 3 parameters distribution (Pearson 3 parameters distribution); 5) Log-Gamma 3 parameters distribution (Log-Pearson 3 parameters distribution); 6) Gumbel distribution; 7) Exponential distribution
#'
#' @param pdf A vector with an ID for the resulting distribution from the KS test. The ID is related to  the numeration of the distribution listed in the description of this function.
#' @param series An array of daily, monthly or annual data, without considering leap days. Possible input dimensions: 2D array: [sub-periods, years] ; [sub-periods + days window, years]
#' @param mu A vector with mean values of each sub-period.
#' @param sigma A vector with standard deviation values of each sub-period.
#' @param skew A vector with skewness values of each sub-period.
#' @param skewy A vector with skewness values of the logarithm of the series of each sub-period.
#'
#' @return prob: A matrix with the non-exceedance probability for value row of the input series.
#' @export
#'
#' @examples getCDF(pdf,series,mu,sigma,skew,skewy)
getCDF <- function(pdf,series,mu,sigma,skew,skewy){

  # Get the number of rows and columns of the series.
  nrows <- dim(series)[1]
  ncols <- dim(series)[2]

  # Compute the cumulative distribution function to the values of each 
  # row, based on the distribution assigned in the getDist function of the
  # climQMBC package.
  prob <- matrix(0,nrows,ncols)
  for (sp in 1:nrows){
    series_sub <- series[sp,]
    if (pdf[sp] == 1){ # i) Normal distribution.
      prob[sp,] <- pnorm(series_sub,mu[sp],sigma[sp])

    } else if (pdf[sp] == 2){ # ii) Log-Normal distribution.
      sigmay <- sqrt(log(1+(sigma[sp]/mu[sp])^2))
      muy <- log(mu[sp])-(sigmay^2)/2
      prob[sp,] <- plnorm(series_sub,muy,sigmay)

    } else if (pdf[sp] == 3){ # iii) Gamma 2 parameters distribution.
      A <- (sigma[sp]^2)/mu[sp]
      B <- (mu[sp]/sigma[sp])^2
      prob[sp,] <- pgamma(series_sub,shape=B,scale=A)

    } else if (pdf[sp] == 4){ # iv) Gamma 3 parameters distribution.
      Bet  <- (2/skew[sp])^2
      Alp  <- sigma[sp]/sqrt(Bet)
      Gam  <- mu[sp]-(Alp*Bet)
      prob[sp,]  <- pgamma((series_sub-Gam),shape=Bet,scale=Alp)

    } else if (pdf[sp] == 5){ # v) Log-Gamma 3 parameters distribution.
      Bety   <- (2/skewy[sp])^2
      sigmay <- sqrt(log(1+(sigma[sp]/mu[sp])^2))
      Alpy   <- sigmay/sqrt(Bety)
      muy    <- log(mu[sp])-(sigmay^2)/2
      Gamy   <- muy - (Alpy*Bety)
      Lnsortdata <- log(series_sub)
      Lnsortdata[!is.finite(Lnsortdata)] <- log(runif(sum(!is.finite(Lnsortdata)))*0.01)
      prob[sp,]  <- pgamma((Lnsortdata-Gamy),shape=Bety,scale=Alpy)

    } else if (pdf[sp] == 6){ # vi) Gumbel distribution.
      Sn <- pi/sqrt(6)
      yn <- 0.5772
      a <- Sn/sigma[sp]
      u <- mu[sp]-(yn/a)
      prob[sp,] <- exp(-exp(-a*(series_sub-u)))

    } else if (pdf[sp] == 7){ # vii) Exponential distribution.
      gamexp <- mu[sp] - sigma[sp]
      prob[sp,] <- pmax(1-exp(-1/sigma[sp]*(series_sub-gamexp)),0)
    }
  }

  # Tail probabilities are set to (0 +  threshold) and (1 - threshold) to
  # avoid numerical errors.
  th <- 10^-3
  prob[prob>1-th] <- 1-th
  prob[prob<th] <- th

  return(prob)
}


#' getCDFinv
#' 
#' This function evaluates the probability, prob, in the respective inverse cumulative distribution function assigned by the Kolmogorov-Smirnov (KS) test in the getDist function of the climQMBC package.
#' 
#' The available distributions are: 1) Normal distribution; 2) Log-Normal distribution; 3) Gamma 2 parameters distribution; 4) Gamma 3 parameters distribution (Pearson 3 parameters distribution); 5) Log-Gamma 3 parameters distribution (Log-Pearson 3 parameters distribution); 6) Gumbel distribution; 7) Exponential distribution
#'
#' @param pdf A vector with an ID for the resulting distribution from the KS test. The ID is related to  the numeration of the distribution listed in the description of this function.
#' @param prob A matrix with the non-exceedance probability for value row of the input series.
#' @param mu A vector with mean values of each sub-period.
#' @param sigma A vector with standard deviation values of each sub-period.
#' @param skew A vector with skewness values of each sub-period.
#' @param skewy A vector with skewness values of the logarithm of the series of each sub-period.
#'
#' @return  xhat: A matrix with the values obtained when the inverse cumulative distribution function is applied to prob.
#' @export
#'
#' @examples getCDFinv(PDF,Taot,mu,sigma,skew,skewy)
getCDFinv <- function(pdf,prob,mu,sigma,skew,skewy){

  # Get the number of rows and columns of the series.
  rows <- dim(prob)[1]
  cols <- dim(prob)[2]

  # Compute the cumulative distribution function to the values of each 
  # row, based on the distribution assigned in the getDist function of the
  # climQMBC package.
  xhat <- matrix(0,rows,cols)
  for (sp in 1:rows){
    if (pdf[sp] == 1){ # i) Normal distribution.
      xhat[sp,] <- qnorm(prob[sp,],mu[sp],sigma[sp])

    } else if (pdf[sp] == 2){ # ii) Log-Normal distribution.
      sigmay <- sqrt(log(1+(sigma[sp]/mu[sp])^2))
      muy <- log(mu[sp])-(sigmay^2)/2
      xhat[sp,] <- qlnorm(prob[sp,],muy,sigmay)

    } else if (pdf[sp] == 3){ # iii) Gamma 2 parameters distribution.
      A <- (sigma[sp]^2)/mu[sp]
      B <- (mu[sp]/sigma[sp])^2
      xhat[sp,] <- qgamma(prob[sp,],shape=B,scale=A)

    } else if (pdf[sp] == 4){ # iv) Gamma 3 parameters distribution.
      Bet  <- (2/skew[sp])^2
      Alp  <- sigma[sp]/sqrt(Bet)
      Gam  <- mu[sp]-(Alp*Bet)
      xhat[sp,]  <- qgamma(prob[sp,],shape=Bet,scale=Alp) + Gam

    } else if (pdf[sp] == 5){ # v) Log-Gamma 3 parameters distribution.
      Bety   <- (2/skewy[sp])^2
      sigmay <- sqrt(log(1+(sigma[sp]/mu[sp])^2))
      Alpy   <- sigmay/sqrt(Bety)
      muy    <- log(mu[sp])-(sigmay^2)/2
      Gamy   <- muy - (Alpy*Bety)
      xhat[sp,]  <- exp(qgamma(prob[sp,],shape=Bety,scale=Alpy) + Gamy)

    } else if (pdf[sp] == 6){ # vi) Gumbel distribution.
      Sn <- pi/sqrt(6)
      yn <- 0.5772
      a <- Sn/sigma[sp]
      u <- mu[sp]-(yn/a)
      xhat[sp,] <- u-log(-log(prob[sp,]))/a

    } else if (pdf[sp] == 7){ # vii) Exponential distribution.
      gamexp <- mu[sp] - sigma[sp]
      xhat[sp,] <- gamexp-(sigma[sp]*log(1-prob[sp,]))
    }
  }

  return(xhat)
}
