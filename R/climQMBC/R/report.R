#' climQMBC example report
#'
#' This function generates two report of the performance of the different methods (QM, DQM, QDM, UQM and SDM) available in the climQMBC package. These reports are based on the mean and standard deviation of the series in the historical and future projected periods.
#'
#' The first report is a summary table with the overall performance of the QM method in the historical period and of the different methods in future projected periods. The methods performance is addressed by comparing its variations in the mean and standard deviation with the ones of the modeled data.
#'
#' The second report consist of three figures. Figure 1 shows the cumulative distribution of the observed, modeled, and corrected series in the historical and future period. This figure also shows the complete time series. Figure 2 and 3 shows the monthly mean and standard deviation, respectively, of each series in the historical and future projected periods. In these two figures, in the future projected periods, the observed series is replaced by an objective series which is computed as the observed series (monthly mean or standard deviation), scaled by the variation between the projected and the historical period of the modeled series. Each projected period is centered in a moving window whose length is equal to the length of the historical period.
#' 
#' NOTE: This routine considers that obs and mod series start in the same day/month/year and are continuous until the end day/month/year.
#'
#' @param obs A column vector of monthly observed data. The length of the column vector should by a multiple of 12. [ndata_obs, 1]
#' @param mod A column vector of monthly modeled or GCM data.The length of the column vector should by a multiple of 12. [ndata_obs, 1]
#' @param var A flag that identifies if data are temperature or precipitation. Temperature:   var = 0 ; Precipitation: var = 1
#' @param mult_change (Optional) A flag that indicates if projected changes should be computed as multiplicative (fut = hist*delta) or additive (fut = hist + delta) changes. mult_change = 1 or True: Multiplicative (default) ; mult_change = 0 or False: Additive
#' @param allow_negatives (Optional) A flag that identifies if data allows negative values and also to replace no-rain values with random small  values (Chadwick et al., 2023) to avoid numerical problems with the probability distribution functions. allow_negatives = 1 or True: Allow negatives (default) ; allow_negatives = 0 or False: Do not allow negative
#' @param fun (Optional) A list of strings with the desired bias correction methods to be reported. If this input is not recieved by the function, all bias correction methods available in the climQMBC package will be reported. The methods supported are: a) 'QM' : Quantile Mapping; b) 'DQM': Detrended Quantile Mapping; c) 'QDM': Quantile Delta Mapping; d) 'UQM': Unbiased Quantile Mapping; Default: fun = ['QM','DQM','QDM','UQM','SDM']
#' @param y_init (Optional) First year of the observed and modeled series (integer). Default: y_init = 0
#' @param y_wind (Optional) A list of integers with the year of the center of the projected periods to be reported. Default: y_wind = [int(y_obs+y_obs/2), int(y_mod-y_obs/2)] This value sets a first projected period just after the end of historical period, and a second projected period just before the end of the modeled series.
#' @param user_pdf (Optional) A flag indicating if the user will define the probability distribution functions (pdf) for the observed and modeled series. The distributions will be the same for all periods and sub-periods. user_pdf = 1 or True: User defines the pdf ; user_pdf = 0 or False: pdf defined by the Kolmogorov -Smirnov test (default). NOTE: The available distributions for pdf_obs and pdf_mod are:  1) Normal ; 2) Log-Normal ; 3) Gamma 2 parameters ; 4) Gamma 3 parameters ; 5) Log-Gamma 3 parameters ; 6) Gumbel ; 7) Exponential
#' @param pdf_obs (Optional) An integer indicating the probability distribution function (pdf) to be used for the observed data. The pdf will be the same for all periods and sub-periods. Default: FALSE
#' @param pdf_mod (Optional) An integer indicating the probability distribution  function (pdf) to be used for the modeled data. The pdf will be the same for all periods and sub-periods. Default: FALSE
#' 
#' @return A list with the five methods implemented in the climQMBC ppackage: QM, DQM, QDM, UQM, and SDM.
#' @export
#'
#' @examples report(obs, mod, SDM_var)
#' @examples report(obs, mod, SDM_var, mult_change, allow_negatives)
#' @examples report(obs, mod, SDM_var,fun=['QDM','UQM','SDM'],y_init = 1979,y_wind = [2035,2060,2080])
report <- function(obs, mod ,SDM_var, mult_change, allow_negatives, fun, y_init, y_wind, user_pdf, pdf_obs, pdf_mod){

  # 0) Get the number of observed and modeled years
  n_obs <- length(obs)
  y_obs <- n_obs/12

  n_mod <- length(mod)
  y_mod = n_mod/12

  # 1) Set non-declared arguments
  if(missing(mult_change)) {
    mult_change <- 1
  }
  if(missing(allow_negatives)) {
    allow_negatives <- 1
  }
  if(missing(fun)) {
    fun <- c('QM','DQM','QDM','UQM','SDM')
  }
  if(missing(y_init)) {
    y_init <- 0
  }
  if(missing(y_wind)) {
    y_wind <- c(floor(y_obs+y_obs/2),floor(y_mod-y_obs/2))
    w_label <- c('First period','Last period')
    rep_head <- c('Variation              Hist.        First        Last')
  } else {
    rep_head <- c('Variation              Hist')
    w_label <- c()
    for (w in 1:length(y_wind)){
      rep_head <- paste(rep_head,'       ',pracma::num2str(y_wind[w],0))
      w_label = c(w_label,pracma::num2str(y_wind[w],0))
    }
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

  y_wind <- y_wind - y_init

  # 2) Apply QM methods
  QM_series <- QM(obs, mod, allow_negatives=allow_negatives, frq='M', user_pdf=user_pdf, pdf_obs=pdf_obs, pdf_mod=pdf_mod)
  DQM_series <- DQM(obs, mod, mult_change=mult_change, allow_negatives=allow_negatives, frq='M', user_pdf=user_pdf, pdf_obs=pdf_obs, pdf_mod=pdf_mod)
  QDM_series <- QDM(obs, mod, mult_change=mult_change, allow_negatives=allow_negatives, frq='M', user_pdf=user_pdf, pdf_obs=pdf_obs, pdf_mod=pdf_mod)
  UQM_series <- UQM(obs, mod, mult_change=mult_change, allow_negatives=allow_negatives, frq='M', user_pdf=user_pdf, pdf_obs=pdf_obs, pdf_mod=pdf_mod)
  SDM_series <- SDM(obs, mod, SDM_var=SDM_var, frq='M')
  
  # 3) Get observed, modeled and bias corrected statistics
  # a) Get obs, mod, and QMs as [12 x n] matrix
  obs_M <- matrix(obs,nrow = 12, ncol = y_obs)
  mod_M <- matrix(mod,nrow = 12, ncol = y_mod)
  QM_M <- matrix(QM_series,nrow = 12, ncol = y_mod)
  DQM_M <- matrix(DQM_series,nrow = 12, ncol = y_mod)
  QDM_M <- matrix(QDM_series,nrow = 12, ncol = y_mod)
  UQM_M <- matrix(UQM_series,nrow = 12, ncol = y_mod)
  SDM_M <- matrix(SDM_series,nrow = 12, ncol = y_mod)

  # b) Mean
  mu_obs <- mean(obs_M)
  mu_mod_h <- mean(mod_M[,1:y_obs])
  mu_mod_f <- mean(mod_M[,(y_obs+1):y_mod])
  mu_QM_h <- mean(QM_M[,1:y_obs])
  mu_QM_f <- mean(QM_M[,(y_obs+1):y_mod])
  mu_DQM_f <- mean(DQM_M[,(y_obs+1):y_mod])
  mu_QDM_f <- mean(QDM_M[,(y_obs+1):y_mod])
  mu_UQM_f <- mean(UQM_M[,(y_obs+1):y_mod])
  mu_SDM_f <- mean(SDM_M[,(y_obs+1):y_mod])

  # c) Monthly mean
  mu_obs_M <- apply(obs_M,1,mean)
  mu_mod_Mh <- apply(mod_M[,1:y_obs],1,mean)
  mu_mod_Mf <- apply(mod_M[,(y_obs+1):y_mod],1,mean)
  mu_QM_Mh <- apply(QM_M[,1:y_obs],1,mean)
  mu_QM_Mf <- apply(QM_M[,(y_obs+1):y_mod],1,mean)
  mu_DQM_Mf <- apply(DQM_M[,(y_obs+1):y_mod],1,mean)
  mu_QDM_Mf <- apply(QDM_M[,(y_obs+1):y_mod],1,mean)
  mu_UQM_Mf <- apply(UQM_M[,(y_obs+1):y_mod],1,mean)
  mu_SDM_Mf <- apply(SDM_M[,(y_obs+1):y_mod],1,mean)

  # d) Standard deviation
  s_obs <- sd(obs_M)
  s_mod_h <- sd(mod_M[,1:y_obs])
  s_mod_f <- sd(mod_M[,(y_obs+1):y_mod])
  s_QM_h <- sd(QM_M[,1:y_obs])
  s_QM_f <- sd(QM_M[,(y_obs+1):y_mod])
  s_DQM_f <- sd(DQM_M[,(y_obs+1):y_mod])
  s_QDM_f <- sd(QDM_M[,(y_obs+1):y_mod])
  s_UQM_f <- sd(UQM_M[,(y_obs+1):y_mod])
  s_SDM_f <- sd(SDM_M[,(y_obs+1):y_mod])

  # e) Monthly standard deviation
  s_obs_M <- apply(obs_M,1,sd)
  s_mod_Mh <- apply(mod_M[,1:y_obs],1,sd)
  s_mod_Mf <- apply(mod_M[,(y_obs+1):y_mod],1,sd)
  s_QM_Mh <- apply(QM_M[,1:y_obs],1,sd)
  s_QM_Mf <- apply(QM_M[,(y_obs+1):y_mod],1,sd)
  s_DQM_Mf <- apply(DQM_M[,(y_obs+1):y_mod],1,sd)
  s_QDM_Mf <- apply(QDM_M[,(y_obs+1):y_mod],1,sd)
  s_UQM_Mf <- apply(UQM_M[,(y_obs+1):y_mod],1,sd)
  s_SDM_Mf <- apply(SDM_M[,(y_obs+1):y_mod],1,sd)

  # 4) Get observed, modeled and bias corrected statistics of the projected
  #    periods
  # a) Projected period mean
  mu_mod_w <- matrix(0,length(y_wind),1)
  mu_QM_w <- matrix(0,length(y_wind),1)
  mu_DQM_w <- matrix(0,length(y_wind),1)
  mu_QDM_w <- matrix(0,length(y_wind),1)
  mu_UQM_w <- matrix(0,length(y_wind),1)
  mu_SDM_w <- matrix(0,length(y_wind),1)

  # b) Projected period monthly mean
  mu_mod_Mw <- matrix(0,12,length(y_wind))
  mu_QM_Mw <- matrix(0,12,length(y_wind))
  mu_DQM_Mw <- matrix(0,12,length(y_wind))
  mu_QDM_Mw <- matrix(0,12,length(y_wind))
  mu_UQM_Mw <- matrix(0,12,length(y_wind))
  mu_SDM_Mw <- matrix(0,12,length(y_wind))

  # c) Projected period standard deviation
  s_mod_w <- matrix(0,length(y_wind),1)
  s_QM_w <- matrix(0,length(y_wind),1)
  s_DQM_w <- matrix(0,length(y_wind),1)
  s_QDM_w <- matrix(0,length(y_wind),1)
  s_UQM_w <- matrix(0,length(y_wind),1)
  s_SDM_w <- matrix(0,length(y_wind),1)

  # d) Projected period monthly standard deviation
  s_mod_Mw <- matrix(0,12,length(y_wind))
  s_QM_Mw <- matrix(0,12,length(y_wind))
  s_DQM_Mw <- matrix(0,12,length(y_wind))
  s_QDM_Mw <- matrix(0,12,length(y_wind))
  s_UQM_Mw <- matrix(0,12,length(y_wind))
  s_SDM_Mw <- matrix(0,12,length(y_wind))

  # Compute projected period statistics
  for (w in 1:length(y_wind)){
    # Projected period mean
    mu_mod_w[w] <- mean(mod_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)])
    mu_QM_w[w] <- mean(QM_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)])
    mu_DQM_w[w] <- mean(DQM_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)])
    mu_QDM_w[w] <- mean(QDM_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)])
    mu_UQM_w[w] <- mean(UQM_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)])
    mu_SDM_w[w] <- mean(SDM_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)])

    # Projected period monthly mean
    mu_mod_Mw[,w] <- apply(mod_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)],1,mean)
    mu_QM_Mw[,w] <- apply(QM_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)],1,mean)
    mu_DQM_Mw[,w] <- apply(DQM_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)],1,mean)
    mu_QDM_Mw[,w] <- apply(QDM_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)],1,mean)
    mu_UQM_Mw[,w] <- apply(UQM_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)],1,mean)
    mu_SDM_Mw[,w] <- apply(SDM_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)],1,mean)

    # Projected period standard deviation
    s_mod_w[w] <- sd(mod_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)])
    s_QM_w[w] <- sd(QM_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)])
    s_DQM_w[w] <- sd(DQM_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)])
    s_QDM_w[w] <- sd(QDM_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)])
    s_UQM_w[w] <- sd(UQM_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)])
    s_SDM_w[w] <- sd(SDM_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)])

    # Projected period monthly standard deviation
    s_mod_Mw[,w] <- apply(mod_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)],1,sd)
    s_QM_Mw[,w] <- apply(QM_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)],1,sd)
    s_DQM_Mw[,w] <- apply(DQM_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)],1,sd)
    s_QDM_Mw[,w] <- apply(QDM_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)],1,sd)
    s_UQM_Mw[,w] <- apply(UQM_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)],1,sd)
    s_SDM_Mw[,w] <- apply(SDM_M[,floor(y_wind[w]-y_obs/2):floor(y_wind[w]+y_obs/2)],1,sd)
  }

  # 5) Get delta mean and delta standard deviation
  dm_mod_w <- matrix(0,length(y_wind),1)
  dm_QM_w <- matrix(0,length(y_wind),1)
  dm_DQM_w <- matrix(0,length(y_wind),1)
  dm_QDM_w <- matrix(0,length(y_wind),1)
  dm_UQM_w <- matrix(0,length(y_wind),1)
  dm_SDM_w <- matrix(0,length(y_wind),1)

  dm_mod_Mw <- matrix(0,12,length(y_wind))
  dm_QM_Mw <- matrix(0,12,length(y_wind))
  dm_DQM_Mw <- matrix(0,12,length(y_wind))
  dm_QDM_Mw <- matrix(0,12,length(y_wind))
  dm_UQM_Mw <- matrix(0,12,length(y_wind))
  dm_SDM_Mw <- matrix(0,12,length(y_wind))

  ds_mod_w <- matrix(0,length(y_wind),1)
  ds_QM_w <- matrix(0,length(y_wind),1)
  ds_DQM_w <- matrix(0,length(y_wind),1)
  ds_QDM_w <- matrix(0,length(y_wind),1)
  ds_UQM_w <- matrix(0,length(y_wind),1)
  ds_SDM_w <- matrix(0,length(y_wind),1)

  ds_mod_Mw <- matrix(0,12,length(y_wind))
  ds_QM_Mw <- matrix(0,12,length(y_wind))
  ds_DQM_Mw <- matrix(0,12,length(y_wind))
  ds_QDM_Mw <- matrix(0,12,length(y_wind))
  ds_UQM_Mw <- matrix(0,12,length(y_wind))
  ds_SDM_Mw <- matrix(0,12,length(y_wind))

  if (mult_change==1){
    dm_obs <- mu_QM_h/mu_obs
    dm_mod <- mu_mod_f/mu_mod_h
    dm_QM <- mu_QM_f/mu_QM_h
    dm_DQM <- mu_DQM_f/mu_QM_h
    dm_QDM <- mu_QDM_f/mu_QM_h
    dm_UQM <- mu_UQM_f/mu_QM_h
    dm_SDM <- mu_SDM_f/mu_QM_h

    ds_obs <- s_QM_h/s_obs
    ds_mod <- s_mod_f/s_mod_h
    ds_QM <- s_QM_f/s_QM_h
    ds_DQM <- s_DQM_f/s_QM_h
    ds_QDM <- s_QDM_f/s_QM_h
    ds_UQM <- s_UQM_f/s_QM_h
    ds_SDM <- s_SDM_f/s_QM_h

    # Scale factor for the objective series in the future projected series.
    dm_mod_M <- mu_mod_Mf/mu_mod_Mh
    ds_mod_M <- s_mod_Mf/s_mod_Mh

    for (w in 1:length(y_wind)){
      dm_mod_w[w] <- mu_mod_w[w]/mu_mod_h
      dm_QM_w[w] <- mu_QM_w[w]/mu_QM_h
      dm_DQM_w[w] <- mu_DQM_w[w]/mu_QM_h
      dm_QDM_w[w] <- mu_QDM_w[w]/mu_QM_h
      dm_UQM_w[w] <- mu_UQM_w[w]/mu_QM_h
      dm_SDM_w[w] <- mu_SDM_w[w]/mu_QM_h

      ds_mod_w[w] <- s_mod_w[w]/s_mod_h
      ds_QM_w[w] <- s_QM_w[w]/s_QM_h
      ds_DQM_w[w] <- s_DQM_w[w]/s_QM_h
      ds_QDM_w[w] <- s_QDM_w[w]/s_QM_h
      ds_UQM_w[w] <- s_UQM_w[w]/s_QM_h
      ds_SDM_w[w] <- s_SDM_w[w]/s_QM_h

      dm_mod_Mw[,w] <- mu_mod_Mw[,w]/mu_mod_Mh
      dm_QM_Mw[,w] <- mu_QM_Mw[,w]/mu_QM_Mh
      dm_DQM_Mw[,w] <- mu_DQM_Mw[,w]/mu_QM_Mh
      dm_QDM_Mw[,w] <- mu_QDM_Mw[,w]/mu_QM_Mh
      dm_UQM_Mw[,w] <- mu_UQM_Mw[,w]/mu_QM_Mh
      dm_SDM_Mw[,w] <- mu_SDM_Mw[,w]/mu_QM_Mh

      ds_mod_Mw[,w] <- s_mod_Mw[,w]/s_mod_Mh
      ds_QM_Mw[,w] <- s_QM_Mw[,w]/s_QM_Mh
      ds_DQM_Mw[,w] <- s_DQM_Mw[,w]/s_QM_Mh
      ds_QDM_Mw[,w] <- s_QDM_Mw[,w]/s_QM_Mh
      ds_UQM_Mw[,w] <- s_UQM_Mw[,w]/s_QM_Mh
      ds_SDM_Mw[,w] <- s_SDM_Mw[,w]/s_QM_Mh
    }
  } else {
    dm_obs <- mu_QM_h-mu_obs
    dm_mod <- mu_mod_f-mu_mod_h
    dm_QM <- mu_QM_f-mu_QM_h
    dm_DQM <- mu_DQM_f-mu_QM_h
    dm_QDM <- mu_QDM_f-mu_QM_h
    dm_UQM <- mu_UQM_f-mu_QM_h
    dm_SDM <- mu_SDM_f-mu_QM_h

    ds_obs <- s_QM_h-s_obs
    ds_mod <- s_mod_f-s_mod_h
    ds_QM <- s_QM_f-s_QM_h
    ds_DQM <- s_DQM_f-s_QM_h
    ds_QDM <- s_QDM_f-s_QM_h
    ds_UQM <- s_UQM_f-s_QM_h
    ds_SDM <- s_SDM_f-s_QM_h

    # Scale factor for the objective series in the future projected series.
    dm_mod_M <- mu_mod_Mf - mu_mod_Mh
    ds_mod_M <- s_mod_Mf - s_mod_Mh

    for (w in 1:length(y_wind)){
      dm_mod_w[w] <- mu_mod_w[w]-mu_mod_h
      dm_QM_w[w] <- mu_QM_w[w]-mu_QM_h
      dm_DQM_w[w] <- mu_DQM_w[w]-mu_QM_h
      dm_QDM_w[w] <- mu_QDM_w[w]-mu_QM_h
      dm_UQM_w[w] <- mu_UQM_w[w]-mu_QM_h
      dm_SDM_w[w] <- mu_SDM_w[w]-mu_QM_h

      ds_mod_w[w] <- s_mod_w[w]-s_mod_h
      ds_QM_w[w] <- s_QM_w[w]-s_QM_h
      ds_DQM_w[w] <- s_DQM_w[w]-s_QM_h
      ds_QDM_w[w] <- s_QDM_w[w]-s_QM_h
      ds_UQM_w[w] <- s_UQM_w[w]-s_QM_h
      ds_SDM_w[w] <- s_SDM_w[w]-s_QM_h

      dm_mod_Mw[,w] <- mu_mod_Mw[,w]-mu_mod_Mh
      dm_QM_Mw[,w] <- mu_QM_Mw[,w]-mu_QM_Mh
      dm_DQM_Mw[,w] <- mu_DQM_Mw[,w]-mu_QM_Mh
      dm_QDM_Mw[,w] <- mu_QDM_Mw[,w]-mu_QM_Mh
      dm_UQM_Mw[,w] <-  mu_UQM_Mw[,w]-mu_QM_Mh
      dm_SDM_Mw[,w] <- mu_SDM_Mw[,w]-mu_QM_Mh

      ds_mod_Mw[,w] <- s_mod_Mw[,w]-s_mod_Mh
      ds_QM_Mw[,w] <- s_QM_Mw[,w]-s_QM_Mh
      ds_DQM_Mw[,w] <- s_DQM_Mw[,w]-s_QM_Mh
      ds_QDM_Mw[,w] <- s_QDM_Mw[,w]-s_QM_Mh
      ds_UQM_Mw[,w] <- s_UQM_Mw[,w]-s_QM_Mh
      ds_SDM_Mw[,w] <- s_SDM_Mw[,w]-s_QM_Mh
    }
  }

  # 6) Display report
  # a) Report description

 pracma::disp('Description of this report:')
 pracma::disp(' ')
 pracma::disp('Description of the table')
 pracma::disp('The first line shows the difference (relative or absolute, according')
 pracma::disp('to the variable analyzed) between the QM method and the observed data')
 pracma::disp('in the historical period. For precipitation, a value of 1 means that')
 pracma::disp('the mean precipitation of the observed data and QM series are the same')
 pracma::disp('for the historical period. For temperature, this is achieved with a')
 pracma::disp('value of 0. The second line shows the headers of the periods reported.')
 pracma::disp('Remember that the moving window is centered in the period and its')
 pracma::disp('length is equal to the length of the historical period. The first ')
 pracma::disp('column (Hist.) shows the difference between the future and the')
 pracma::disp('historical period. The next columns show the difference between the')
 pracma::disp('corresponding period and the historical period. The third and the')
 pracma::disp('following lines show the difference in each period for the')
 pracma::disp('corresponding series. As a general rule, it is desirable that the')
 pracma::disp('difference between the future periods and the historical period of the')
 pracma::disp('bias corrected series match the differences of the modeled series.')
 pracma::disp(' ')
 pracma::disp('Description of the figures')
 pracma::disp('Figure 1: (left) shows the cumulative distribution functions of each')
 pracma::disp('series split between historical and future period. (right) shows the')
 pracma::disp('observed, modeled, and corrected series.')
 pracma::disp('Figure 2: shows the monthly mean of the historical and future period.')
 pracma::disp('Additionally, monthly mean of the projected periods is displayed.')
 pracma::disp('Figure 3: shows the monthly standard deviation of the historical and')
 pracma::disp('future period. Additionally, monthly mean of the projected periods is')
 pracma::disp('displayed.')
 pracma::disp('* In Fig. 2 and 3, the obj (red) line represents objective values.')
 pracma::disp('This objective series is computed as the observed series (monthly mean')
 pracma::disp('or standard deviation), scaled by the variation between the projected')
 pracma::disp('and the historical period of the modeled series. Each projected period')
 pracma::disp('is centered in  a moving window whose length is equal to the length of')
 pracma::disp('the historical period.')
 pracma::disp(' ')

  # Format each line of the table report
  n = 3 # Number of decimals to which round the values
  sp = c('-',' ',' ') # A small trick to prevent negative values to use more
                     # space than positive values
  rep_mu_M <- paste('Modeled          :', ' ', sp[sign(dm_mod)+2], pracma::num2str(abs(round(dm_mod,n))))
  rep_mu_QM <- paste('QM               :', ' ', sp[sign(dm_QM)+2], pracma::num2str(abs(round(dm_QM,n))))
  rep_mu_DQM <- paste('DQM              :', ' ', sp[sign(dm_DQM)+2], pracma::num2str(abs(round(dm_DQM,n))))
  rep_mu_QDM <- paste('QDM              :', ' ', sp[sign(dm_QDM)+2], pracma::num2str(abs(round(dm_QDM,n))))
  rep_mu_UQM <- paste('UQM              :', ' ', sp[sign(dm_UQM)+2], pracma::num2str(abs(round(dm_UQM,n))))
  rep_mu_SDM <- paste('SDM              :', ' ', sp[sign(dm_SDM)+2], pracma::num2str(abs(round(dm_SDM,n))))

  rep_s_M <- paste('Modeled          :', ' ', sp[sign(ds_mod)+2], pracma::num2str(abs(round(ds_mod,n))))
  rep_s_QM <- paste('QM               :', ' ', sp[sign(ds_QM)+2], pracma::num2str(abs(round(ds_QM,n))))
  rep_s_DQM <- paste('DQM              :', ' ', sp[sign(ds_DQM)+2], pracma::num2str(abs(round(ds_DQM,n))))
  rep_s_QDM <- paste('QDM              :', ' ', sp[sign(ds_QDM)+2], pracma::num2str(abs(round(ds_QDM,n))))
  rep_s_UQM <- paste('UQM              :', ' ', sp[sign(ds_UQM)+2], pracma::num2str(abs(round(ds_UQM,n))))
  rep_s_SDM <- paste('SDM              :', ' ', sp[sign(ds_SDM)+2], pracma::num2str(abs(round(ds_SDM,n))))

  for (w in 1:length(y_wind)){
    rep_mu_M <- paste(rep_mu_M, '  ; ', sp[sign(dm_mod_w[w])+2], pracma::num2str(abs(round(dm_mod_w[w],n))))
    rep_mu_QM <- paste(rep_mu_QM, '  ; ', sp[sign(dm_QM_w[w])+2], pracma::num2str(abs(round(dm_QM_w[w],n))))
    rep_mu_DQM <- paste(rep_mu_DQM, '  ; ', sp[sign(dm_DQM_w[w])+2], pracma::num2str(abs(round(dm_DQM_w[w],n))))
    rep_mu_QDM <- paste(rep_mu_QDM, '  ; ', sp[sign(dm_QDM_w[w])+2], pracma::num2str(abs(round(dm_QDM_w[w],n))))
    rep_mu_UQM <- paste(rep_mu_UQM, '  ; ', sp[sign(dm_UQM_w[w])+2], pracma::num2str(abs(round(dm_UQM_w[w],n))))
    rep_mu_SDM <- paste(rep_mu_SDM, '  ; ', sp[sign(dm_SDM_w[w])+2], pracma::num2str(abs(round(dm_SDM_w[w],n))))

    rep_s_M <- paste(rep_s_M, '  ; ', sp[sign(ds_mod_w[w])+2], pracma::num2str(abs(round(ds_mod_w[w],n))))
    rep_s_QM <- paste(rep_s_QM, '  ; ', sp[sign(ds_QM_w[w])+2], pracma::num2str(abs(round(ds_QM_w[w],n))))
    rep_s_DQM <- paste(rep_s_DQM, '  ; ', sp[sign(ds_DQM_w[w])+2], pracma::num2str(abs(round(ds_DQM_w[w],n))))
    rep_s_QDM <- paste(rep_s_QDM, '  ; ', sp[sign(ds_QDM_w[w])+2], pracma::num2str(abs(round(ds_QDM_w[w],n))))
    rep_s_UQM <- paste(rep_s_UQM, '  ; ', sp[sign(ds_UQM_w[w])+2], pracma::num2str(abs(round(ds_UQM_w[w],n))))
    rep_s_SDM <- paste(rep_s_SDM, '  ; ', sp[sign(ds_SDM_w[w])+2], pracma::num2str(abs(round(ds_SDM_w[w],n))))
  }

  # b) Table
 pracma::disp('Mean')
 pracma::disp(paste('QM performance in historical period  :', ' ', sp[sign(dm_obs)+2], pracma::num2str(abs(round(dm_obs,n)))))
 pracma::disp(rep_head)
 pracma::disp(rep_mu_M)
  if ('QM' %in% fun){
   pracma::disp(rep_mu_QM)
  }
  if ('DQM' %in% fun){
   pracma::disp(rep_mu_DQM)
  }
  if ('QDM' %in% fun){
   pracma::disp(rep_mu_QDM)
  }
  if ('UQM' %in% fun){
   pracma::disp(rep_mu_UQM)
  }
  if ('SDM' %in% fun){
   pracma::disp(rep_mu_SDM)
  }
 pracma::disp(' ')

 pracma::disp('Standard deviation')
 pracma::disp(paste('QM performance in historical period  :', ' ', sp[sign(ds_obs)+2], pracma::num2str(abs(round(ds_obs,n)))))
 pracma::disp(rep_head)
 pracma::disp(rep_s_M)
  if ('QM' %in% fun){
   pracma::disp(rep_s_QM)
  }
  if ('DQM' %in% fun){
   pracma::disp(rep_s_DQM)
  }
  if ('QDM' %in% fun){
   pracma::disp(rep_s_QDM)
  }
  if ('UQM' %in% fun){
   pracma::disp(rep_s_UQM)
  }
  if ('SDM' %in% fun){
   pracma::disp(rep_s_SDM)
  }
 pracma::disp(' ')

  # c) Figures
  series_mn <- min(c(obs,mod,QM_series,DQM_series,QDM_series,UQM_series,SDM_series))
  series_mx <- max(c(obs,mod,QM_series,DQM_series,QDM_series,UQM_series,SDM_series))

  mu_mn <- min(c(mu_mod_Mw,mu_QM_Mw,mu_DQM_Mw,mu_QDM_Mw,mu_UQM_Mw,mu_SDM_Mw))
  mu_mx <- max(c(mu_mod_Mw,mu_QM_Mw,mu_DQM_Mw,mu_QDM_Mw,mu_UQM_Mw,mu_SDM_Mw))

  s_mn <- min(c(s_mod_Mw,s_QM_Mw,s_DQM_Mw,s_QDM_Mw,s_UQM_Mw,s_SDM_Mw))
  s_mx <- max(c(s_mod_Mw,s_QM_Mw,s_DQM_Mw,s_QDM_Mw,s_UQM_Mw,s_SDM_Mw))

  # Figure 1: Cumulative probabilities and time series
  # Get empirical cumulative distribution functions
  xo <- sort(obs)
  fo <- ecdf(xo)(xo)
  xm <- sort(mod[1:length(obs)])
  fm <- ecdf(xm)(xm)
  xm_ <- sort(mod[length(obs)+1:length(mod)])
  fm_ <- ecdf(xm_)(xm_)
  xq <- sort(QM_series[1:length(obs)])
  fq <- ecdf(xq)(xq)
  xq1 <- sort(QM_series[length(obs)+1:length(mod)])
  fq1 <- ecdf(xq1)(xq1)
  xq2 <- sort(DQM_series[length(obs)+1:length(mod)])
  fq2 <- ecdf(xq2)(xq2)
  xq3 <- sort(QDM_series[length(obs)+1:length(mod)])
  fq3 <- ecdf(xq3)(xq3)
  xq4 <- sort(UQM_series[length(obs)+1:length(mod)])
  fq4 <- ecdf(xq4)(xq4)
  xq5 <- sort(SDM_series[length(obs)+1:length(mod)])
  fq5 <- ecdf(xq5)(xq5)

  lgnd <- c('Obs_h','Mod_h','Mod_f','QM_h')
  clr <- c('red','blue','blue','black')
  mstyl <- c(1,1,2,1)
  # Plot empirical cumulative distribution function
  dev.new()
  par(mfrow=c(1,2))
  plot(xo,fo,col='red',type='l',xlab=c('Temperature (C)','Precipitation (mm)')[mult_change+1],ylab='Probability',main='Empirical cumulative distribution function',ylim=c(0,1),xlim=c(series_mn,series_mx))
  grid()
  par(new=TRUE)
  lines(xm,fm,col='blue')
  par(new=TRUE)
  lines(xm_,fm_,col='blue',lty=2)
  par(new=TRUE)
  lines(xq,fq,col='black')
  if ('QM' %in% fun){
    par(new=TRUE)
    lines(xq1,fq1,col='black',lty=2)
    lgnd <- c(lgnd,'QM_f')
    clr <- c(clr,'black')
    mstyl <- c(mstyl,2)
  }
  if ('DQM' %in% fun){
    par(new=TRUE)
    lines(xq2,fq2,col='green',lty=2)
    lgnd <- c(lgnd,'DQM')
    clr <- c(clr,'green')
    mstyl <- c(mstyl,2)
  }
  if ('QDM' %in% fun){
    par(new=TRUE)
    lines(xq3,fq3,col='yellow',lty=2)
    lgnd <- c(lgnd,'QDM')
    clr <- c(clr,'yellow')
    mstyl <- c(mstyl,2)
  }
  if ('UQM' %in% fun){
    par(new=TRUE)
    lines(xq4,fq4,col='cyan',lty=2)
    lgnd <- c(lgnd,'UQM')
    clr <- c(clr,'cyan')
    mstyl <- c(mstyl,2)
  }
  if ('SDM' %in% fun){
    par(new=TRUE)
    lines(xq5,fq5,col='magenta',lty=2)
    lgnd <- c(lgnd,'SDM')
    clr <- c(clr,'magenta')
    mstyl <- c(mstyl,2)
  }
  legend('bottomright',legend=lgnd,lty=mstyl,col=clr,cex=0.75)

  lgnd <- c('Mod')
  clr <- c('blue')
  # Plot series
  plot(mod,col='blue',type='l',xlab='Month since starting date',ylab=c('Temperature (°C)','Precipitation (mm)')[mult_change+1],main=paste(c('Temperature','Precipitation')[mult_change+1],'time series'), ylim=c(series_mn,series_mx))
  grid()
  if ('QM' %in% fun){
    par(new=TRUE)
    lines(QM_series,col='black')
    lgnd <- c(lgnd,'QM')
    clr <- c(clr,'black')
  }
  if ('DQM' %in% fun){
    par(new=TRUE)
    lines(DQM_series,col='green')
    lgnd <- c(lgnd,'DQM')
    clr <- c(clr,'green')
  }
  if ('QDM' %in% fun){
    par(new=TRUE)
    lines(QDM_series,col='yellow')
    lgnd <- c(lgnd,'QDM')
    clr <- c(clr,'yellow')
  }
  if ('UQM' %in% fun){
    par(new=TRUE)
    lines(UQM_series,col='cyan')
    lgnd <- c(lgnd,'UQM')
    clr <- c(clr,'cyan')
  }
  if ('SDM' %in% fun){
    par(new=TRUE)
    lines(SDM_series,col='magenta')
    lgnd <- c(lgnd,'SDM')
    clr <- c(clr,'magenta')
  }
  par(new=TRUE)
  lines(obs,col='red')
  lgnd <- c(lgnd,'Obs.')
  clr <- c(clr,'red')
  legend('topright',legend=lgnd,lty=1,col=clr,cex=0.75,bg='white')


  # Figure 2: Monthly means
  # Historical period
  dev.new()
  par(mfrow=c(1,length(y_wind)+2))
  # Plot series
  plot(mu_mod_Mh,col='blue',type='l',xlab='Month',ylab='',xlim=c(1,12),ylim=c(mu_mn,mu_mx))
  title(main='Monthly mean of the historical period')
  grid()
  par(new=TRUE)
  lines(mu_QM_Mh,col='black')
  par(new=TRUE)
  lines(mu_obs_M,col='red')
  legend('top',legend=c('Mod_h','QM_h','Obs'),lty=1,col=c('blue','black','red'),bg='white')

  # Future period
  lgnd <- c('Mod')
  clr <- c('blue')
  plot(mu_mod_Mf,col='blue',type='l',xlab='Month',ylab='',xlim=c(1,12),ylim=c(mu_mn,mu_mx))
  title(main='Monthly mean of the future period')
  grid()
  if ('QM' %in% fun){
    par(new=TRUE)
    lines(mu_QM_Mf,col='black')
    lgnd <- c(lgnd,'QM')
    clr <- c(clr,'black')
  }
  if ('DQM' %in% fun){
    par(new=TRUE)
    lines(mu_DQM_Mf,col='green')
    lgnd <- c(lgnd,'DQM')
    clr <- c(clr,'green')
  }
  if ('QDM' %in% fun){
    par(new=TRUE)
    lines(mu_QDM_Mf,col='yellow')
    lgnd <- c(lgnd,'QDM')
    clr <- c(clr,'yellow')
  }
  if ('UQM' %in% fun){
    par(new=TRUE)
    lines(mu_UQM_Mf,col='cyan')
    lgnd <- c(lgnd,'UQM')
    clr <- c(clr,'cyan')
  }
  if ('SDM' %in% fun){
    par(new=TRUE)
    lines(mu_SDM_Mf,col='magenta')
    lgnd <- c(lgnd,'SDM')
    clr <- c(clr,'magenta')
  }
  par(new=TRUE)
  if (mult_change==1){
    lines(mu_obs_M*dm_mod_M,col='red')
  } else {
    lines(mu_obs_M+dm_mod_M,col='red')
  }
  lgnd <- c(lgnd,'Obj.')
  clr <- c(clr,'red')
  legend('topright',legend=lgnd,lty=1,col=clr,bg='white')

  # Projected periods
  for (w in 1:length(y_wind)){
    lgnd <- c('Mod')
    clr <- c('blue')
    plot(mu_mod_Mw[,w],col='blue',type='l',xlab='Month',ylab='',xlim=c(1,12),ylim=c(mu_mn,mu_mx))
    title(main=w_label[w])
    grid()
    if ('QM' %in% fun){
      par(new=TRUE)
      lines(mu_QM_Mw[,w],col='black')
      lgnd <- c(lgnd,'QM')
      clr <- c(clr,'black')
    }
    if ('DQM' %in% fun){
      par(new=TRUE)
      lines(mu_DQM_Mw[,w],col='green')
      lgnd <- c(lgnd,'DQM')
      clr <- c(clr,'green')
    }
    if ('QDM' %in% fun){
      par(new=TRUE)
      lines(mu_QDM_Mw[,w],col='yellow')
      lgnd <- c(lgnd,'QDM')
      clr <- c(clr,'yellow')
    }
    if ('UQM' %in% fun){
      par(new=TRUE)
      lines(mu_UQM_Mw[,w],col='cyan')
      lgnd <- c(lgnd,'UQM')
      clr <- c(clr,'cyan')
    }
    if ('SDM' %in% fun){
      par(new=TRUE)
      lines(mu_SDM_Mw[,w],col='magenta')
      lgnd <- c(lgnd,'SDM')
      clr <- c(clr,'magenta')
    }
    par(new=TRUE)
    if (mult_change==1){
      lines(mu_obs_M*dm_mod_Mw[,w],col='red')
    } else {
      lines(mu_obs_M+dm_mod_Mw[,w],col='red')
    }
    lgnd <- c(lgnd,'Obj.')
    clr <- c(clr,'red')
    legend('topright',legend=lgnd,lty=1,col=clr,bg='white')
  }


  # Figure 3: Monthly standard deviation
  # Historical period
  dev.new()
  par(mfrow=c(1,length(y_wind)+2))
  # Plot series
  plot(s_mod_Mh,col='blue',type='l',xlab='Month',ylab='',xlim=c(1,12),ylim=c(s_mn,s_mx))
  title(main='Monthly std. dev. of the historical period')
  grid()
  par(new=TRUE)
  lines(s_QM_Mh,col='black')
  par(new=TRUE)
  lines(s_obs_M,col='red')
  legend('top',legend=c('Mod_h','QM_h','Obs'),lty=1,col=c('blue','black','red'),bg='white')

  # Future period
  lgnd <- c('Mod')
  clr <- c('blue')
  plot(s_mod_Mf,col='blue',type='l',xlab='Month',ylab='',xlim=c(1,12),ylim=c(s_mn,s_mx))
  title(main='Monthly std. dev. of the future period')
  grid()
  if ('QM' %in% fun){
    par(new=TRUE)
    lines(s_QM_Mf,col='black')
    lgnd <- c(lgnd,'QM')
    clr <- c(clr,'black')
  }
  if ('DQM' %in% fun){
    par(new=TRUE)
    lines(s_DQM_Mf,col='green')
    lgnd <- c(lgnd,'DQM')
    clr <- c(clr,'green')
  }
  if ('QDM' %in% fun){
    par(new=TRUE)
    lines(s_QDM_Mf,col='yellow')
    lgnd <- c(lgnd,'QDM')
    clr <- c(clr,'yellow')
  }
  if ('UQM' %in% fun){
    par(new=TRUE)
    lines(s_UQM_Mf,col='cyan')
    lgnd <- c(lgnd,'UQM')
    clr <- c(clr,'cyan')
  }
  if ('SDM' %in% fun){
    par(new=TRUE)
    lines(s_SDM_Mf,col='magenta')
    lgnd <- c(lgnd,'SDM')
    clr <- c(clr,'magenta')
  }
  par(new=TRUE)
  if (mult_change==1){
    lines(s_obs_M*ds_mod_M,col='red')
  } else {
    lines(s_obs_M+ds_mod_M,col='red')
  }
  lgnd <- c(lgnd,'Obj.')
  clr <- c(clr,'red')
  legend('topright',legend=lgnd,lty=1,col=clr,bg='white')

  # Projected periods
  for (w in 1:length(y_wind)){
    lgnd <- c('Mod')
    clr <- c('blue')
    plot(s_mod_Mw[,w],col='blue',type='l',xlab='Month',ylab='',xlim=c(1,12),ylim=c(s_mn,s_mx))
    title(main=w_label[w])
    grid()
    if ('QM' %in% fun){
      par(new=TRUE)
      lines(s_QM_Mw[,w],col='black')
      lgnd <- c(lgnd,'QM')
      clr <- c(clr,'black')
    }
    if ('DQM' %in% fun){
      par(new=TRUE)
      lines(s_DQM_Mw[,w],col='green')
      lgnd <- c(lgnd,'DQM')
      clr <- c(clr,'green')
    }
    if ('QDM' %in% fun){
      par(new=TRUE)
      lines(s_QDM_Mw[,w],col='yellow')
      lgnd <- c(lgnd,'QDM')
      clr <- c(clr,'yellow')
    }
    if ('UQM' %in% fun){
      par(new=TRUE)
      lines(s_UQM_Mw[,w],col='cyan')
      lgnd <- c(lgnd,'UQM')
      clr <- c(clr,'cyan')
    }
    if ('SDM' %in% fun){
      par(new=TRUE)
      lines(s_SDM_Mw[,w],col='magenta')
      lgnd <- c(lgnd,'SDM')
      clr <- c(clr,'magenta')
    }
    par(new=TRUE)
    if (mult_change==1){
      lines(s_obs_M*ds_mod_Mw[,w],col='red')
    } else {
      lines(s_obs_M+ds_mod_Mw[,w],col='red')
    }
    lgnd <- c(lgnd,'Obj.')
    clr <- c(clr,'red')
    legend('topright',legend=lgnd,lty=1,col=clr,bg='white')
  }


  return(list(QM_series,DQM_series,QDM_series,UQM_series,SDM_series))
}
