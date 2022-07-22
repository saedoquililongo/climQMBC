# climQMBC for R
A package with multiple bias correction methods for climatic variables, including the QM, DQM, QDM, UQM, and SDM methods.

Authors: Sebastian Aedo Quililongo (1*), Cristian Chadwick (2), Fernando Gonzalez-Leiva (3), and Jorge Gironas (3)

(1) Centro de Cambio Global UC, Pontificia Universidad Catolica de Chile, Santiago, Chile.\
(2) Faculty of Engineering and Sciences, Universidad Adolfo Ibanez, Santiago, Chile. \
(3) Department of Hydraulics and Environmental Engineering, Pontificia Universidad Catolica de Chile, Santiago, Chile. \

*Maintainer contact: slaedo@uc.cl\
Version: 0.1.0, updated Jul 2022\
Written in R version 3.6.1

## Dependencies
Imports needed:
- e1071 (>=1.7-6)
- fitdistrplus (>=1.1-6)
- pracma (>=2.3.3)

## Installation
As this package is not available in CRAN, you must make sure to previosuly install the required dependencies. Also, to install packages with windows, tha package Rtools is needed (https://cran.r-project.org/bin/windows/Rtools/).

- Save the .tar.gz file and identify the path to the directory (avoid using complex paths with dots "." which could raise an error indicating that R can not find the package)
- To install the package, run the following line in R, where < dir > is the path to the directory where the .tar.gz file is located:
```R
install.packages("<dir>\\climQMBC_0.1.0.tar.gz",repos=NULL,type="source")
```
- Verify that the package was succesfully installed by running the following line in R:
```R
library(climQMBC)
```

We suggest the user to look at the climQMBC_tester_R.R script. This script show examples of the capacities of the package and suggestions of its use.


## Version history
### Version 0.1.0 (current version on github)
- **General:** Incorporated the capacity to replace low precipitation values (values below `pp_threshold`) with random values between 0 and `pp_factor`. This capacity also solves the problem when the moving window has only cero values (dry months or arid zones), which causes numerical problems with the invers Cumulative Distribution Functions. This can be disabled by defining `pp_threshold = 0`. `pp_threshold` and `pp_factor` are optional arguments defined with values of `1` and `1/100`, respectively.
- **QDM:** Incorporated the capacity to truncate the scaling factor for precipitations in the QDM method to avoid overestimations in dry periods caused by numerical indefinitions of the scaling factor. The scaling factor is truncated to a value defined with the parameter `rel_change_th` when the denominator of the scaling factor is below the parameter `inv_mod_th`. This can be disabled by defining `inv_mod_th = 0`. `rel_change_th` and `inv_mod_th` are optional arguments defined with values of `2` and `pp_threshold`, respectively.
- **(R) Report function:** Fixed the xlims and ylims of the CDF and monthly series plots, respectively.


### Version 0.0.0
First release of the climQMBC package.