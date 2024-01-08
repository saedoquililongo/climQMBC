# climQMBC for Python
A package with multiple bias correction methods for climatic variables, including the QM, DQM, QDM, UQM, and SDM methods.

Authors: Sebastian Aedo Quililongo (1*), Cristian Chadwick (2), Fernando Gonzalez-Leiva (3), and Jorge Gironas (3, 4)

(1) Stockholm Environment Institute, Latin America Centre, Bogota, Colombia.\
(2) Faculty of Engineering and Sciences, Universidad Adolfo Ibanez, Santiago, Chile. \
(3) Department of Hydraulics and Environmental Engineering, Pontificia Universidad Catolica de Chile, Santiago, Chile.\
(4) Centro de Cambio Global UC, Pontificia Universidad Catolica de Chile, Santiago, Chile. 

*Maintainer contact: sebastian.aedo.q@gmail.com\
Version: 0.1.0, updated Jun 2022\
Written in Python version 3.7.9

## Dependencies
Imports needed:
- numpy (>=1.20.3)
- scipy (>=1.6.3)
- pandas (>=1.2.4) [For the `report` function]

## Importing
As this package is not available to install via pip, you must make sure to previosuly install the required dependencies. To import the climQMBC package we suggest the following procedure:

- Download the py_.zip file and extract the climQMBC folder located inside the zip file.
- Add these lines at the beginning of your code to add the path to the python scripts, where < dir > is the path to the climQMBC folder (Note that if the climQMBC folder is located in the same path as your script, you could skip this step):
```Python
import sys
sys.path.append('<dir>')
```

- To import the methods, you could run:
```Python
from climQMBC.methods import QM, DQM, QDM, UQM, SDM
```
- To import the report function, you could run:
```Python
from climQMBC.report import report
```

We suggest the user to look at the climQMBC_tester_Python.py script. This script show examples of the capacities of the package and suggestions of its use.

## Version history
### Version 0.1.1 (current version on github; no changes in Python since version 0.1.0)
- **(R) formatQM:** Fixed missing `var` argument in formatQM function and in each method when formatQM is called.


### Version 0.1.0
- **General:** Incorporated the capacity to replace low precipitation values (values below `pp_threshold`) with random values between 0 and `pp_factor`. This capacity also solves the problem when the moving window has only cero values (dry months or arid zones), which causes numerical problems with the invers Cumulative Distribution Functions. This can be disabled by defining `pp_threshold = 0`. `pp_threshold` and `pp_factor` are optional arguments defined with values of `1` and `1/100`, respectively.
- **QDM:** Incorporated the capacity to truncate the scaling factor for precipitations in the QDM method to avoid overestimations in dry periods caused by numerical indefinitions of the scaling factor. The scaling factor is truncated to a value defined with the parameter `rel_change_th` when the denominator of the scaling factor is below the parameter `inv_mod_th`. This can be disabled by defining `inv_mod_th = 0`. `rel_change_th` and `inv_mod_th` are optional arguments defined with values of `2` and `pp_threshold`, respectively.
- **(R) Report function:** Fixed the xlims and ylims of the CDF and monthly series plots, respectively.


### Version 0.0.0
First release of the climQMBC package.