# climQMBC
A package with multiple bias correction methods for climatic variables, including the QM, DQM, QDM, UQM, and SDM methods.

Authors: Sebastian Aedo Quililongo (1*), Cristian Chadwick (2), Fernando Gonzalez-Leiva (3), and Jorge Gironas (3, 4)

(1) Stockholm Environment Institute, Latin America Centre, Bogota, Colombia.\
(2) Faculty of Engineering and Sciences, Universidad Adolfo Ibanez, Santiago, Chile. \
(3) Department of Hydraulics and Environmental Engineering, Pontificia Universidad Catolica de Chile, Santiago, Chile.\
(4) Centro de Cambio Global UC, Pontificia Universidad Catolica de Chile, Santiago, Chile. 

*Maintainer contact: sebastian.aedo.q@gmail.com\
Version: 1.0.0, updated Apr 2024\


## Introduction
Bias Correction is the process of eliminating GCM or RCM biases based on observed values. Several bias correction bethods have been developed to correct the bias of GCM or RCM simulations; however, based on the formulation of these methods, they preserve different changes of the modeled series and perform different under certain conditions.

The climQMBC package makes available five different Quantile Mapping methods for bias correction of climatic variables. These are the Quantile Mapping (QM), Detrended Quantile Mapping (DQM), Quantile Delta Mapping (DQM), Unbiased Quantile Mapping (UQM), and Scaled Distribution Mapping (SDM) methods. The package also provides a **_report_** function to compare the performance or evaluate the tradeoffs between each method, based on the conservation of changes in the mean and standard deviation of the GMCs.

The climQMBC package considers bias correction of monthly and annual series of variables whose changes are relative (precipitation) or absolute (temperature).

Bias correction of modeled series (GCMs outputs) is performed at a single location based on the observed values. As a general approach, the QM methods looks for relationships between the modeled series and observed data and applies it to the modeled series to obtain an unbiased climate series in the historical and future periods. The climQMBC package applies these relations as a continuum, where for each year it considers a window with length equal to the number of years of the historical period and ends in the analyzed year to correct that specific year. This process is repeated for each year of the future period. Note that this continuum methodology does not apply for the standard QM method.

For each period, the package selects the most appropiate probability ditribution function based in the Kolmogorov-Smirnov test. If the series frequency is daily or monthly, a distribution function is assigned for each day or month independently. Distributions are fitted by the Method of Moments. The available distributions are:
- Normal distribution
- Log-normal distribution
- Gamma 2 parameters distribution
- Gamma 3 parameters distribution
- Log-Gamma distribution
- Gumbel distribution
- Exponential distribution

For precipitation values, only strictly positive distributions are considered. For temperature values, strictly positive distributions are discarded if there are negative values in the analyzed period.

This distribution adjustment does not apply for the SDM method as it considers the gamma distribution for precipitation and normal distribution for temperatures, fitted by the Maximum likelihood Estimates Method.

The climQMBC package has the capacity to replace low precipitation values with random values. Considering a default threshold of 1 mm, values below this threshold are replaced with random values between 0 and 1/100. These parameters can be modified by the user and even the user can disable this capacity.

For daily bias correction, the user can define a moving window centered on the day to be corrected to fit the corresponding probability distribution function and the corresponding statistics. Additionally, for precipitation it adjust the number of rain-days by setting a threshold on the modeled series to fit the observed number of rain-days and performs the bias correction process to the rain values. If the number of rain-values to fit a probability distribution function is below 30, it will consider some no-rain values to fit te ditribution.

## Getting started
### Installation or importing

Within the [Python](https://github.com/saedoquililongo/climQMBC/tree/main/Python), [Matlab](https://github.com/saedoquililongo/climQMBC/tree/main/Matlab), and [R](https://github.com/saedoquililongo/climQMBC/tree/main/R) directories you will find the installation and importing procedure for each programming language.


### Quick start

Within the [Python](https://github.com/saedoquililongo/climQMBC/tree/main/Python), [Matlab](https://github.com/saedoquililongo/climQMBC/tree/main/Matlab), and [R](https://github.com/saedoquililongo/climQMBC/tree/main/R) directories you will find scripts with examples of the capacities of the package and suggestions of its use.

- Python: [climQMBC_tester_Python.py](https://github.com/saedoquililongo/climQMBC/blob/main/Python/climQMBC_tester_Python.py)
- Matlab: [climQMBC_tester_Matlab.m](https://github.com/saedoquililongo/climQMBC/blob/main/Matlab/climQMBC_tester_Matlab.m)
- R: [climQMBC_tester_R.R](https://github.com/saedoquililongo/climQMBC/blob/main/R/climQMBC_tester_R.R)

Additionally, the [Example_notebook.ipynb](https://github.com/saedoquililongo/climQMBC/blob/main/Example_notebook.ipynb) is a notebook showing the usage of the package in Python and examples of the results of the methods and functions implemented.


## References
#### climQMBC:
[Aedo, S., Chadwick, C., González-Leiva, F., and Gironás, J. (2021). climQMBC a new package to Bias Correct Climatic Variables While Preserving raw GCM Changes in the Mean and Standard Deviation For R, Python and Matlab. American Geophysical Union fall meeting, 2021.](https://agu2021fallmeeting-agu.ipostersessions.com/Default.aspx?s=52-1C-3B-41-27-7C-34-E2-DE-3F-55-24-7B-0C-34-48)

#### QM, DQM and QDM:
[Cannon, A. J., S. R. Sobie, and T. Q. Murdock. (2015). Bias correction of GCM precipitation by quantile mapping: How well do methods preserve changes in quantiles and extremes? J. Climate, 28(17), 6938-6959, https://doi.org/10.1175/JCLI-D-14-00754.1](https://doi.org/10.1175/JCLI-D-14-00754.1)

#### UQM:
[Chadwick, C., Gironás, J., González-Leiva, F., and Aedo, S. (2023). Bias adjustment to preserve changes in variability: the unbiased mapping of GCM changes. Hydrological Sciences Journal, https://doi.org/10.1080/02626667.2023.2201450](https://doi.org/10.1080/02626667.2023.2201450)

[Chadwick, C., Gironás, J., González-Leiva, F., and Aedo, S. (2021). Bias Correction and the Importance of Preserving the Changes in Variability: Unbiased Quantile Mapping (UQM) a Novel Bias Correction. American Geophysical Union fall meeting, 2021.](https://agu2021fallmeeting-agu.ipostersessions.com/default.aspx?s=48-67-54-07-35-60-D9-5B-8D-0C-9C-6C-1C-1A-92-EE)

#### SDM:
[Switanek, B. M., Troch, P., Castro, L. C., Leuprecht, A., Chang, H. I., Mukherjee, R., and Demaria, M. C. E. (2017). Scaled distribution mapping: A bias correction method that preserves raw climate model projected changes. Hydrology &amp; Earth System Sciences, 21, 2649-2666, https://doi.org/10.5194/hess-21-2649-2017.](https://doi.org/10.5194/hess-21-2649-2017)


## Version history
### Version 1.0.0 (current version on github)
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
