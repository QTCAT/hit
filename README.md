# Hierarchical Inference Testing

The current built and test status for Linux 
[![Travis-CI Build Status](https://travis-ci.org/jrklasen/hit.png?branch=master)]
(https://travis-ci.org/jrklasen/hit) and for Windows 
[![AppVeyor Build Status]
  (https://ci.appveyor.com/api/projects/status/github/jrklasen/hit?branch=master&svg=true)]
(https://ci.appveyor.com/project/jrklasen/hit).

## Description:
Hierarchical inference testing (HIT) for linear models with correlated covariates. 
HIT is furthermore applicable to high-dimensional settings. For details see:

**Mandozzi, J. and Buehlmann, P. (2013)**. *Hierarchical testing in the high-dimensional 
setting with correlated variables*. To appear in the Journal of the American Statistical 
Association. Preprint [arXiv:1312.5556](http://arxiv.org/abs/1312.5556).

## Install:

The package can be installed from an R console via `devtools`. If you haven't yet 
`devtools` installed, you have to do so first.

```R
# install.packages("devtools")
devtools::install_github("jrklasen/hit")
```

If you want to install a specific release e.g. `v0.0-3`, use: 

```R
devtools::install_github("jrklasen/hit", ref="v0.0-3")
```

## Example:

The `hit`-function example gives an overview of the functionality of the package 
and can be accessed once the package is loaded (please be aware, that running the 
example can take a few seconds).

```R
library(hit)
example(hit)
```

[![License]
  (http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)]
(http://www.gnu.org/licenses/gpl-2.0.html)
&copy; 2015 JR Klasen
