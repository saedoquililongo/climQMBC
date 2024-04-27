# climQMBC for R
A package with multiple bias correction methods for climatic variables, including the QM, DQM, QDM, UQM, and SDM methods.

Authors: Sebastian Aedo Quililongo (1*), Cristian Chadwick (2), Fernando Gonzalez-Leiva (3), and Jorge Gironas (3, 4)

(1) Stockholm Environment Institute, Latin America Centre, Bogota, Colombia.\
(2) Faculty of Engineering and Sciences, Universidad Adolfo Ibanez, Santiago, Chile. \
(3) Department of Hydraulics and Environmental Engineering, Pontificia Universidad Catolica de Chile, Santiago, Chile.\
(4) Centro de Cambio Global UC, Pontificia Universidad Catolica de Chile, Santiago, Chile. 

*Maintainer contact: sebastian.aedo.q@gmail.com\
Version: 1.0, updated Apr 2024\
Written in R version 4.2.2.

## Dependencies
Imports needed:
- abind (>=1.4-5)
- e1071 (>=1.7-6)
- fitdistrplus (>=1.1-6)
- pracma (>=2.3.3)

## Installation
As this package is not available in CRAN, you must make sure to previosuly install the required dependencies. Also, to install packages with windows, tha package Rtools is needed (https://cran.r-project.org/bin/windows/Rtools/).

- Save the .tar.gz file and identify the path to the directory (avoid using complex paths with dots "." which could raise an error indicating that R can not find the package)
- To install the package, run the following line in R, where < dir > is the path to the directory where the .tar.gz file is located:
```R
install.packages("<dir>\\climQMBC_1.0.tar.gz",repos=NULL,type="source")
```
- Verify that the package was succesfully installed by running the following line in R:
```R
library(climQMBC)
```

We suggest the user to look at the climQMBC_tester_R.R script. This script show examples of the capacities of the package and suggestions of its use.


## Version history
### Version 1.0 (current version on github)
- **General:** Bias correction of daily data is now available, considering years of 365 days and a moving window centered in the day to be corrected. The length of the moving window is defined by the user by the argument `day_win` (Default value is 1). If `allow_negatives=0` (like precipitation) only rainday values will be corrected and used to define the most appropiate probability distribution function and scaling factors or projected changes. Additionally, for `allow_negatives=0` a `pp_threshold_mod` will be calculated internally to match the number of raindays in the historical period between the observations and modeled series.
- **General:** Improved the calculation of moving windows and projected changes, making the functions faster.
- **All but SDM:** The `var` argument is replaced with `allow_negatives` and `mult_change`. `allow_negatives` is an optional argument that takes values of 1 or 0 to define if the data allow negative values or not, respectively. If `allow_negatives=0`, probability distribution functions that allow negative values are dropped and values below `pp_threshold` are replaced with random values between 0 and `pp_threshold*pp_factor` (ideal for precipitation). `mult_change` is an optional argument that takes values of 1 or 0 to define if scaling factors or projected changes are incorporated multiplicatively or additively, respectively. Default values are `mult_change=1` and `allow_negatives=1`. As a suggestion, for precipitation use `mult_change=1` and `allow_negatives=0`; for temperature in Celsius use `mult_change=0` and `allow_negatives=1`; for temperature in Kelvin use `mult_change=1` and `allow_negatives=1`.
- **All but SDM:** Three optional arguments to allow the user to define the probability distribution functions to be used instead of using the KS-test were added. `user_pdf` is a flag indicating if the user will define the probability distribution functions for the observed and modeled series. The distributions will be the same for all periods and sub-periods. `user_pdf=1` indicates that the user will define the pdf. `pdf_obs` and `pdf_mod` is an integer indicating the probability distribution function (pdf) to be used for the observed data. The available distributions for `pdf_obs` and `pdf_mod` are: 0) Normal, 1) Log-Normal, 2) Gamma 2 parameters, 3) Gamma 3 parameters, 4) Log-Gamma 3 parameters, 5) Gumbel and 6) Exponential.
- **All but SDM:** A warning will be print when the probability distribution function defined by the KS-test fails the test (being the better does not mean it is good). Only one warning per use of the function.
- **SDM:** The `var` argument is replaced by `SDM_var`.
- **SDM:** Fixed bugs related to dry periods.
- **General:** Updated the sample data with daily and monthly data for the same location.


### Version 0.1.1
- **(R) formatQM:** Fixed missing `var` argument in formatQM function and in each method when formatQM is called.


### Version 0.1.0
- **General:** Incorporated the capacity to replace low precipitation values (values below `pp_threshold`) with random values between 0 and `pp_factor`. This capacity also solves the problem when the moving window has only cero values (dry months or arid zones), which causes numerical problems with the invers Cumulative Distribution Functions. This can be disabled by defining `pp_threshold = 0`. `pp_threshold` and `pp_factor` are optional arguments defined with values of `1` and `1/100`, respectively.
- **QDM:** Incorporated the capacity to truncate the scaling factor for precipitations in the QDM method to avoid overestimations in dry periods caused by numerical indefinitions of the scaling factor. The scaling factor is truncated to a value defined with the parameter `rel_change_th` when the denominator of the scaling factor is below the parameter `inv_mod_th`. This can be disabled by defining `inv_mod_th = 0`. `rel_change_th` and `inv_mod_th` are optional arguments defined with values of `2` and `pp_threshold`, respectively.
- **(R) Report function:** Fixed the xlims and ylims of the CDF and monthly series plots, respectively.


### Version 0.0.0
First release of the climQMBC package.