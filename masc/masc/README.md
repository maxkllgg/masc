
 masc: Matching and Synthetic Control Estimator.
========

### Installation

Install using `devtools::install_github("maxkllgg/masc")`. 

### Help

See `help(masc)` for details and examples using this package.

### Details

This package implements the matching and synthetic control (masc) estimator of Kellogg, Mogstad, Pouliot, and Torgovitsky (2019), hereafter Kellogg et al., (2019).
The motivation for this estimator is the complimentary behavior of matching and synthetic control estimators: matching estimators
limit interpolation bias but are sensitive to extrapolation bias, and the synthetic control (SC) estimator limits extrapolation bias
but is sensitive to interpolation bias. The SC estimator *interpolates* by using some convex weighted average of untreated units to
create a synthetic untreated unit with pre-treatment characteristics similar to that of the treated unit. To do this, the SC estimator
may put weight on control units very dissimilar to the treated unit on pre-treatment characteristics. Interpolation bias arises,
 as shown in Kellogg et al. (2019), if the conditional mean of the outcome is non-linear in pre-treatment characteristics.
Extrapolation bias arises, on the other hand, when the counterfactual constructed for the treated unit is not
identical to the treated unit on pre-treatment characteristics. For more details, see Kellogg et al. (2019).

### Recommended packages

The default implementation of this package uses gurobi, which offers free academic licenses.
For more information, see the following link:

 https://www.gurobi.com/products/optimization-modeling-language-resources-support/r/



### Equations and Computation
See [our computational note](https://github.com/maxkllgg/masc/blob/master/computation/computation.pdf).

### Questions and Issues
Please feel free to submit an issue on the Github repository.
