# Hierarchical Inference Testing (HIT)

Master branch build and test passing status at Linux:
[![Travis-CI Build Status](https://travis-ci.org/jrklasen/hit.png?branch=master)](https://travis-ci.org/jrklasen/hit), and at Windows:
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/jrklasen/hit?branch=master&svg=true)](https://ci.appveyor.com/project/jrklasen/hit).

## Install

The package can be installed from an R console via `devtools`. If you haven't 
yet `devtools` installed, you have to do so first.

    # install.packages("devtools")
    devtools::install_github("jrklasen/hit")

If you want to install a specific release e.g. `v0.0-2`, use: 
`devtools::install_github("jrklasen/hit", ref="v0.0-2")`

## Example

The `hit`-function example gives an overview of the functionality of the 
package and can be accessed once the package is loaded (please be aware, 
that running the example can take a few seconds).

    library(hit)
    example(hit)

## Algorithm
The algorithm is described in:

**Mandozzi, J. and Buehlmann, P. (2013)**. *Hierarchical testing in the 
high-dimensional setting with correlated variables.* To appear in the Journal 
of the American Statistical Association. Preprint [arXiv:1312.5556](http://arxiv.org/abs/1312.5556).

