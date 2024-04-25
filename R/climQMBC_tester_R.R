# We hope this script will help the user understand how to use the main
# functionalities of the climQMBC package. Alongside this script, in a folder
# called Sample_data, the user will find sample datasets for monthly
# precipitation (mm) and mean temperature (C). The datasets include modeled and
# observed data. The files name structure is xxx_yyy.csv, where xxx can be mod
# or obs for modeled or observed data, respectively, and yyy can be tmp or pp
# for temperature or precipitation data, respectively. The historical period of
# the sample datasets begin in 1979 and has a length of 36 years (1979 to 2014,
# including both years), and the modeled period begins in 1979 and has a length
# of 122 years (1979 to 2100, including both years).
#
# This script is divided into three examples. The first two examples can be used
# to evaluate the performance of each bias correction method available in the
# climQMBC package.
# Example 1 shows how to use the report function with the minimum number of
# inputs. The five methods available in the climQMBC package will be reported,
# and the future projected windows displayed are the first period after the
# historical period and the last period before the end of the modeled period.
# Remember that the projected periods length is equal to the length of the
# historical period.
# Example 2 shows how to use the report function for specific bias correction
# methods and projected periods. The Quantile Delta Mapping (QDM), Unbiased
# Quantile Mapping (UQM), and Scaled Distribution Mapping (SDM) methods will be
# reported. The report will analyze the projected periods centered in 2035 2060
# and 2080.
# Example 3 shows how each bias correction method available in the climQMBC
# package should be called. The outputs of each function are columns vector with
#  monthly corrected data.
#
# Feel free to uncomment each example, modify the periods and try your own
# datasets.
#
#
# Written by Sebastian Aedo Quililongo (1*)
#            Cristian Chadwick         (2)
#            Fernando Gonzalez-Leiva   (3)
#            Jorge Gironas             (3, 4)
# 
#  (1) Stockholm Environment Institute, Latin America Centre, Bogota, Colombia
#  (2) Faculty of Engineering and Sciences, Universidad Adolfo Ibanez, Santiago,
#      Chile
#  (3) Department of Hydraulics and Environmental Engineering, Pontificia
#      Universidad Catolica de Chile, Santiago, Chile
#  (4) Centro de Cambio Global UC, Pontificia Universidad Catolica de Chile,
#      Santiago, Chile
# 
# *Maintainer contact: sebastian.aedo.q@gmail.com
# Revision: 0, updated Dec 2021

## Recordar a usuarios que deben descargar libreria abind

# # If the package is not installed, save the .tar.gz file in the same directory
# # as this script. The following lines will install the package and import it.
# install.packages(paste(getwd(),'/','climQMBC_0.1.1.tar.gz',sep=''),repos=NULL,type='source')


# library(climQMBC)
source(paste(getwd(),'/climQMBC/R/methods.R',sep=''))
source(paste(getwd(),'/climQMBC/R/utils.R',sep=''))

# Casos mensuales y anuales
# # # Load observed and model data. Remember that for temperature, var = 0, and
# # # for precipitation, var = 1.
# # var <- 1
# # 
# # if (var==1){
# #   agg <- sum
# #   var_txt <- 'pp' 
# # } else {
# #   var_txt <- 'tmax' 
# #   agg <- mean
# # }
# 
var_txt <- 'pr'

obs <- read.csv(paste(getwd(),'/../Sample_data/obs_',var_txt,'_M.csv',sep=''))
obs <- matrix(obs[,var_txt])

mod <- read.csv(paste(getwd(),'/../Sample_data/mod_',var_txt,'_M.csv',sep=''))
mod <- matrix(mod[,var_txt])

# frq <- 'D' #'M' for monthly data; 'A' for annually data
# QM_series <- QM(obs,mod,allow_negatives=0,frq=frq)
# DQM_series <- DQM(obs,mod,mult_change=1,allow_negatives=0,frq=frq)
# QDM_series <- QDM(obs,mod,mult_change=1,allow_negatives=0,frq=frq)
# UQM_series <- UQM(obs,mod,mult_change=1,allow_negatives=0,frq=frq)
# SDM_series <- SDM(obs,mod,SDM_var=1,frq=frq)





## Example 1
#   Example 1 shows how to use the report function with the minimum number
#   of inputs. The five methods available in the climQMBC package will be
#   reported, and the future projected windows displayed are the first
#   period after the historical period and the last period before the end
#   of the modeled period. Remember that the projected periods length is
#   equal to the length of the historical period.
# 
# rep_series <- report(obs,mod,var)
# QM_series <- rep_series[[1]]
# DQM_series <- rep_series[[2]]
# QDM_series <- rep_series[[3]]
# UQM_series <- rep_series[[4]]
# SDM_series <- rep_series[[5]]


## Example 2
#   Example 2 shows how to use the report function for specific bias
#   correction methods and projected periods. The Quantile Delta Mapping
#   (QDM), Unbiased Quantile Mapping (UQM), and Scaled Distribution Mapping
#   (SDM) methods will be reported. The report will analyze the projected
#   periods centered in 2035 2060 and 2080.

# rep_series <- report(obs,mod,var,fun=c('QDM','UQM','SDM'),y_init=1979,y_wind=c(2035,2060,2080))
# QM_series <- rep_series[[1]]
# DQM_series <- rep_series[[2]]
# QDM_series <- rep_series[[3]]
# UQM_series <- rep_series[[4]]
# SDM_series <- rep_series[[5]]


## Example 3
#   Example 3 shows how each bias correction method available in the
#   climQMBC package should be called. The outputs of each function are
#   columns vector with monthly corrected data.
#  
# frq <- 'M' #'M' for monthly data; 'A' for annually data
# QM_series <- QM(obs,mod,var,frq)
# DQM_series <- DQM(obs,mod,var,frq)
# QDM_series <- QDM(obs,mod,var,frq)
# UQM_series <- UQM(obs,mod,var,frq)
# SDM_series <- SDM(obs,mod,var,frq)




###########
frq <- 'D'
variable <- 'pr'
allow_negatives <- 0
SDM_var <- 1

day_win <- 15
pp_threshold <- 1
pp_factor <- 1/10000

# Casos diarios
obs <- read.csv(paste(getwd(),'/../Sample_data/obs_pr_D.csv',sep=''))
obs <- matrix(obs$pr)

mod <- read.csv(paste(getwd(),'/../Sample_data/mod_pr_D.csv',sep=''))
mod <- matrix(mod$pr)
# 
# 
frq <- 'D' #'M' for monthly data; 'A' for annually data
qm_series <- QM(obs,mod,allow_negatives=allow_negatives, frq=frq, pp_threshold=pp_threshold, pp_factor=pp_factor,day_win=day_win, user_pdf=TRUE, pdf_obs=5, pdf_mod=5)
# DQM_series <- DQM(obs,mod,mult_change=1,allow_negatives=0,frq=frq,pp_threshold=1, pp_factor=1/(100*100), win=15)
# QDM_series <- QDM(obs,mod,mult_change=1,allow_negatives=0,frq=frq,pp_threshold=1, pp_factor=1/(100*100), win=15)
# UQM_series <- UQM(obs,mod,mult_change=1,allow_negatives=0,frq=frq,pp_threshold=1, pp_factor=1/(100*100), win=15)
# # SDM_series <- SDM(obs,mod,SDM_var=1,frq=frq)
# 
plot(qm_series)

write.csv(qm_series, "../../A.csv")