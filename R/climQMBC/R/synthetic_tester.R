#' get_synthetic_gamma
#'
#' Gets a random set of yn values of mean mu and standard deviation sigma based on a gamma distribution
#'
#' @param mu Mean of the series
#' @param sigma Standard deviation of the series
#' @param yn Number of years of the sample
#'
#' @return An array with the sampled values
#' @export
#'
#' @examples get_synthetic_gamma(mu, sigma, yn)
get_synthetic_gamma <- function(mu, sigma, yn){
  k <- (mu/sigma)^2
  th <- (sigma^2)/mu
  return(rgamma(yn,shape=k,scale=th))
}


#' get_synthetic_normal
#'
#' Gets a random set of yn values of mean mu and standard deviation sigma based on a normal distribution
#'
#' @param mu Mean of the series
#' @param sigma Standard deviation of the series
#' @param yn Number of years of the sample
#'
#' @return An array with the sampled values
#' @export
#'
#' @examples get_synthetic_normal(mu, sigma, yn)
get_synthetic_normal <- function(mu, sigma, yn){
  return(rnorm(yn,mu,sigma))
}


#' add_bias
#'
#' Adds bias based on a equi-probability transform
#'
#' @param series A column vector with the original series
#' @param mu_bias The relative or absolute bias to be added in the mean
#' @param sigma_bias The relative or absolute bias to be added in the standard deviation
#' @param mult_change (Optional) A flag that indicates if projected changes should be computed as multiplicative (fut = hist*delta) or  additive (fut = hist + delta) changes. mult_change = 1 or True: Multiplicative (default) ; mult_change = 0 or False: Additive
#'
#' @return An array with the biased series
#' @export
#'
#' @examples add_bias(eries, mu_bias, sigma_bias, mult_change)
add_bias <- function(series, mu_bias, sigma_bias, mult_change){
  mu <- mean(series)
  sigma <- sd(series)

  if (mult_change == 1){
    mu_scaled <- mu*(1 + mu_bias)
    sigma_scaled <- sigma*(1 + sigma_bias)
  } else {
    mu_scaled <- mu + mu_bias
    sigma_scaled <- sigma + sigma_bias
  }

  biased_series <- mu_scaled + (series - mu)*sigma_scaled/sigma
  return(biased_series)
}


#' get_tester_stats
#'
#' This function computes the tester statistics, including the NSE, KGE and ratios (or difference) in the mean, std. dev., skewness and several percentiles.
#'
#' @param bc_series A column vector with the biased series
#' @param obj_series A column vector with the objective series
#' @param mult_change (Optional) A flag that indicates if projected changes should be computed as multiplicative (fut = hist*delta) or  additive (fut = hist + delta) changes. mult_change = 1 or True: Multiplicative (default) ; mult_change = 0 or False: Additive
#'
#' @return A list with the tester statistics
#' @export
#'
#' @examples get_tester_stats(bc_series, obj_series, mult_change)
get_tester_stats <- function(bc_series, obj_series, mult_change){
  bc_mean <- mean(bc_series)
  obj_mean <- mean(obj_series)

  a <- sum((bc_series - bc_mean)*(obj_series - obj_mean))
  b <- sqrt(sum((bc_series - bc_mean)^2))
  c <- sqrt(sum((obj_series - obj_mean)^2))

  r2 <- (a/(b*c))^2

  if (mult_change==1){
    mu_ratio <- mean(bc_series)/mean(obj_series)
    sigma_ratio <- sd(bc_series)/sd(obj_series)
    sk_ratio <- e1071::skewness(bc_series,type=1)/e1071::skewness(obj_series, type=1)
    p05_ratio <- quantile(bc_series,0.05,names=FALSE)/quantile(obj_series,0.05,names=FALSE)
    p10_ratio <- quantile(bc_series,0.10,names=FALSE)/quantile(obj_series,0.10,names=FALSE)
    p25_ratio <- quantile(bc_series,0.25,names=FALSE)/quantile(obj_series,0.25,names=FALSE)
    p50_ratio <- quantile(bc_series,0.50,names=FALSE)/quantile(obj_series,0.50,names=FALSE)
    p75_ratio <- quantile(bc_series,0.75,names=FALSE)/quantile(obj_series,0.75,names=FALSE)
    p90_ratio <- quantile(bc_series,0.90,names=FALSE)/quantile(obj_series,0.90,names=FALSE)
    p95_ratio <- quantile(bc_series,0.95,names=FALSE)/quantile(obj_series,0.95,names=FALSE)

    kge <- 1-sqrt((1-r2)^2 + (1-mu_ratio)^2 + (1-sigma_ratio)^2)
    nse <- 1-sum((bc_series-obj_series)^2)/sum((mean(obj_series)-obj_series)^2)

  } else {
    mu_ratio <- mean(bc_series)-mean(obj_series)
    sigma_ratio <- sd(bc_series)-sd(obj_series)
    sk_ratio <- e1071::skewness(bc_series,type=1)-e1071::skewness(obj_series, type=1)
    p05_ratio <- quantile(bc_series,0.05,names=FALSE)-quantile(obj_series,0.05,names=FALSE)
    p10_ratio <- quantile(bc_series,0.10,names=FALSE)-quantile(obj_series,0.10,names=FALSE)
    p25_ratio <- quantile(bc_series,0.25,names=FALSE)-quantile(obj_series,0.25,names=FALSE)
    p50_ratio <- quantile(bc_series,0.50,names=FALSE)-quantile(obj_series,0.50,names=FALSE)
    p75_ratio <- quantile(bc_series,0.75,names=FALSE)-quantile(obj_series,0.75,names=FALSE)
    p90_ratio <- quantile(bc_series,0.90,names=FALSE)-quantile(obj_series,0.90,names=FALSE)
    p95_ratio <- quantile(bc_series,0.95,names=FALSE)-quantile(obj_series,0.95,names=FALSE)

    kge <- sqrt((1-r2)^2 + (1-mu_ratio)^2 + (1-sigma_ratio)^2)
    nse <- sum((bc_series-obj_series)^2)/sum((mean(obj_series)-obj_series)^2)
  }
  return(list(nse,kge,mu_ratio,sigma_ratio,sk_ratio,
              p05_ratio,p10_ratio,p25_ratio,p50_ratio,
              p75_ratio,p90_ratio,p95_ratio))
}


#' synthetic_tester
#'
#' This function This function runs the synthetic tester function to compare the performance of different methods based on a synthetic example rather than the original data, to compare the theoretical performance.
#'
#' @param obs A column vector of monthly observed data.
#' @param mod A column vector of monthly modeled or GCM data.
#' @param allow_negatives A flag that identifies if data allows negative values and also to replace no-rain values with random small  values (Chadwick et al., 2023) to avoid numerical problems with the probability distribution functions. allow_negatives = 1 or True: Allow negatives (default) ; allow_negatives = 0 or False: Do not allow negative
#' @param mult_change A flag that indicates if projected changes should be computed as multiplicative (fut = hist*delta) or additive (fut = hist + delta) changes. mult_change = 1 or True: Multiplicative (default) ; mult_change = 0 or False: Additive
#' @param SDM_var A flag that identifies if data are temperature or precipitation. Temperature:   var = 0 ; Precipitation: var = 1
#' @param yn Number of years of the sample
#'
#' @return A table with the tester statistics
#' @export
#'
#' @examples synthetic_tester(obs, mod, allow_negatives, mult_change, SDM_var, yn)
synthetic_tester <- function(obs, mod, allow_negatives, mult_change, SDM_var, yn){
  if(missing(yn)) {
    yn <- 10000
  }

  # Get series statistics
  mu_obs_h <- mean(obs)
  sigma_obs_h <- sd(obs)

  mu_mod_h <- mean(mod[1:length(obs)])
  sigma_mod_h <- sd(mod[1:length(obs)])

  mu_mod_f <- mean(mod[length(obs):length(mod)])
  sigma_mod_f <- sd(mod[length(obs):length(mod)])

  # Get bias and future values
  if (mult_change == 1){
    mu_bias <- mu_mod_h/mu_obs_h - 1
    sigma_bias <- sigma_mod_h/sigma_obs_h - 1

    mu_change <- mu_mod_f/mu_mod_h - 1
    sigma_change <- sigma_mod_f/sigma_mod_h - 1

    mu_obj_f <- mu_obs_h*(1+mu_change)
    sigma_obj_f <- sigma_obs_h*(1+sigma_change)
  } else {
    mu_bias <- mu_mod_h-mu_obs_h
    sigma_bias <- sigma_mod_h-sigma_obs_h

    mu_change <- mu_mod_f-mu_mod_h
    sigma_change <- sigma_mod_f-sigma_mod_h

    mu_obj_f <- mu_obs_h + mu_change
    sigma_obj_f <- sigma_obs_h + sigma_change
  }


  # Get synthetic series
  if (allow_negatives == 0){
    obs_h <- get_synthetic_gamma(mu_obs_h,sigma_obs_h,yn)
    obs_f <- get_synthetic_gamma(mu_obj_f,sigma_obj_f,yn)
  } else {
    obs_h <- get_synthetic_normal(mu_obs_h,sigma_obs_h,yn)
    obs_f <- get_synthetic_normal(mu_obj_f,sigma_obj_f,yn)
  }

  mod_h <- add_bias(obs_h,mu_bias,sigma_bias,mult_change)
  mod_f <- add_bias(obs_f,mu_bias,sigma_bias,mult_change)


  # Concat historical and future period
  ## Future period is concated twice so that only the second future is
  ## analyzed, without moving windows
  obs <- matrix(obs_h)
  mod <- t(rbind(c(mod_h,mod_f,mod_f)))


  # Apply bias correction methods
  QM_h_series <- QM(obs_h,mod_h,allow_negatives=allow_negatives,frq='A')
  SDM_h_series <- SDM(obs_h,mod_h,SDM_var=SDM_var,frq='A')
  QM_f_series <- QM(obs,mod,allow_negatives=allow_negatives,frq='A')[(2*yn+1):length(mod)]
  DQM_f_series <- DQM(obs,mod,allow_negatives=allow_negatives,mult_change=mult_change,frq='A')[(2*yn+1):length(mod)]
  QDM_f_series <- QDM(obs,mod,allow_negatives=allow_negatives,mult_change=mult_change,frq='A')[(2*yn+1):length(mod)]
  UQM_f_series <- UQM(obs,mod,allow_negatives=allow_negatives,mult_change=mult_change,frq='A')[(2*yn+1):length(mod)]
  SDM_f_series <- SDM(obs,mod,SDM_var=SDM_var,frq='A')[(2*yn+1):length(mod)]

  li_stats_qm_h <- get_tester_stats(QM_h_series,obs_h,mult_change)
  li_stats_sdm_h <- get_tester_stats(SDM_h_series,obs_h,mult_change)
  li_stats_qm_f <- get_tester_stats(QM_f_series,obs_f,mult_change)
  li_stats_dqm_f <- get_tester_stats(DQM_f_series,obs_f,mult_change)
  li_stats_qdm_f <- get_tester_stats(QDM_f_series,obs_f,mult_change)
  li_stats_uqm_f <- get_tester_stats(UQM_f_series,obs_f,mult_change)
  li_stats_sdm_f <- get_tester_stats(SDM_f_series,obs_f,mult_change)

  tester_stats <- rbind(li_stats_qm_h,li_stats_qm_f,li_stats_qm_f,li_stats_dqm_f,li_stats_qdm_f,li_stats_uqm_f,li_stats_sdm_f)

  rownames(tester_stats) <- c('QM_h','SDM_h','QM_f','DQM_f','QDM_f','UQM_f','SDM_f')
  if (mult_change == 1){
    colnames(tester_stats) <- c('NSE','KGE','Mean','Std. Dev','Skew','P05','P10','P25','P50','P75','P90','P95')
  } else {
    colnames(tester_stats) <- c('1-NSE','1-KGE','Mean','Std. Dev','Skew','P05','P10','P25','P50','P75','P90','P95')
  }



  return(tester_stats)
}
