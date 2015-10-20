# Hierarchical Inference Testing

[![CRAN](http://www.r-pkg.org/badges/version/hit)](http://cran.r-project.org/package=hit)


Master branch bild and test status for Linux (Mac)
[![Travis-CI Build Status](https://travis-ci.org/QTCAT/hit.png?branch=master)]
(https://travis-ci.org/QTCAT/hit) 
and for Windows 
 [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/kttq4x98q6hra6ct/branch/master?svg=true)](https://ci.appveyor.com/project/jrklasen/hit)
.



## Description:
Hierarchical inference testing (HIT) for linear models with correlated 
covariates. HIT is furthermore applicable to high-dimensional settings. For 
details see:

**Mandozzi, J. and Buehlmann, P. (2015)**. *Hierarchical testing in the 
high-dimensional setting with correlated variables*. Journal of the American 
Statistical Association. Preprint 
[arXiv:1312.5556](http://arxiv.org/abs/1312.5556).

## Install:
The package can be installed from CRAN,

```R
install.packages("hit")

```

or via [`devtools`](https://github.com/hadley/devtools#updating-to-the-latest-version-of-devtools), if you haven't yet `devtools` installed, you have to do so first.

```R
# install.packages("devtools")
devtools::install_github("QTCAT/hit")
```

## Example:
The `hit`-function example gives an overview of the functionality of the 
package and can be accessed once the package is loaded.

```R
library(hit)
example(hit)
```

--------------------------------------------------------------------------------
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
&copy; 2015 JR Klasen
