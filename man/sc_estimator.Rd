% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimator.R
\name{sc_estimator}
\alias{sc_estimator}
\title{solving for standard synthetic control weights}
\usage{
sc_estimator(
  treated,
  donors,
  treated.covariates,
  donors.covariates,
  treatment,
  nogurobi = FALSE,
  psdtol = 1e+06,
  BarConvTol = 1e-08,
  BarIterLimit = 1e+05,
  sigf.ipop = 5,
  margin.ipop = 5e-04,
  ...
)
}
\arguments{
\item{donors}{See \link{masc}.}

\item{nogurobi}{A logical value. If true, uses \link[LowRankQP]{LowRankQP} to solve the synthetic control estimator,
rather than \code{gurobi}.}

\item{treated:}{See \link{masc}.}

\item{covariates.donors:}{See \link{masc}.}

\item{treated.donors:}{See \link{masc}.}

\item{treatment:}{An integer. The period T' in which forecasting begins (either the true treatment
period or the first period after a cross-validation fold).}
}
\value{
A list. The named component \code{weights.sc} contains the vector of synthetic control weights.
 Weights are ordered in the same manner as the columns in \code{Z0}.
 The \code{objval.sc} component contains the objective value (pre-period fit) of the synthetic control.
}
\description{
Solves for synthetic control weights up to the designated treatment period,
using only information in the outcome paths (no other covariates).
}
\references{
Kellogg, M., M. Mogstad, G. Pouliot, and A. Torgovitsky. Combining Matching and Synthetic Control to Trade
 off Biases from Extrapolation and Interpolation. Working Paper, 2019.
}
\seealso{
Other masc functions: 
\code{\link{cv_masc}()},
\code{\link{masc_by_phi}()},
\code{\link{masc}()},
\code{\link{solve_masc}()}
}
\concept{masc functions}
