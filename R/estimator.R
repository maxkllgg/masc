##################################
#OPTIMIZATION PROBLEM
#################################

#Nearest Neighbor weight vector
Wbar <- function(donors, treated, treatment, K) {
  dist <-
    apply((donors[1:(treatment - 1), ] - as.vector(treated)[1:(treatment -
                                                                 1)]) ^ 2, MARGIN = 2, FUN = sum)
  W <- rep(0, length(dist))
  W[which(dist %in% sort(dist)[1:K])] <-
    1 / length(which(dist %in% sort(dist)[1:K]))
  return(W)
}


#' solving for matching and synthetic control weights
#'
#' Solves separately for synthetic control weights and matching weights for a given
#' matching estimator, using information in \code{data} up to the treatment period
#' designated in that object.
#'
#' @param data A \code{list} containing three named elements: \describe{
#' \item{donors:}{A \eqn{TxN} matrix of outcome paths for untreated units, each column being a control unit.}
#' \item{treated:}{A \eqn{Tx1} matrix of outcome paths for the treated unit.}
#' \item{treatment:}{An integer. The period T' in which treatment occurs (\eqn{T'<T}).}
#' }
#'
#' @param tune_pars A \code{list} containing one element:
#' \describe{\item{m:}{an integer identifying the candidate matching (nearest neighbor) estimator.}}
#'
#' @return A list containing the weights associated with the two estimators. Weights are ordered
#' in the same manner as the columns in \code{data}. The \code{weights.sc} and \code{weights.match}
#' named components contain the vector of synthetic control and nearest neighbor weights respectively.
#' @family masc functions
#' @references Kellogg, M., M. Mogstad, G. Pouliot, and A. Torgovitsky. Combining Matching and Synthetic Control to Trade
#'  off Biases from Extrapolation and Interpolation. Working Paper, 2019.
#' @export
solve_masc <-
  function(data,
           tune_pars) {
    ####Pulling objects out of lists###
    treatment <- data$treatment
    donors <- data$donors[1:(treatment - 1), ]
    treated <- as.matrix(data$treated[1:(treatment - 1), ])
    m <- tune_pars$m
    ###########################################
    ##DEFINING FIXED GUROBI PARAMETERS:
    psdtol = 1e6
    BarConvTol = 1e-8
    BarIterLimit = 1e5
    ##########################################
    obj.synthetic <-
      function(donors, treated, treatment) {
          return(-2 * t(treated) %*% donors)
      }
    objcon.synthetic <-
      function(donors, treated, treatment) {
          return(t(treated) %*% treated)
      }
    Q.synthetic <-
      function(donors, treated, treatment) {
          return(t(donors) %*% donors)
      }



    params <- list()

    modelreg <- list(
      A = matrix(
        rep(1, length(donors[1, ]))
        ,
        nrow = 1,
        ncol = ncol(donors),
        byrow = TRUE
      ),
      sense = c('='),
      rhs = c(1),
      lb = rep(0, length(donors[1, ])),
      ub = rep(1, length(donors[1, ])),
      obj = obj.synthetic(
        donors = donors,
        treated = treated,
        treatment = treatment),
      objcon = objcon.synthetic(
        donors = donors,
        treated = treated,
        treatment = treatment),
      Q = Q.synthetic(
        donors = donors,
        treated = treated,
        treatment = treatment)
    )



    continue <- 0

    tparams <-
      gurobi::gurobi(
        modelreg,
        params = list(
          OutputFlag = 0,
          PSDTol = psdtol,
          BarConvTol = BarConvTol,
          BarIterLimit = BarIterLimit
        )
      )
    tparams <- rename(tparams, c('x' = 'pars'))
    if (tparams$status != 'OPTIMAL') {
      print('WARNING: FAILED TO SOLVE PROBLEM')
      continue
    }
    params$weights.sc <- tparams$pars
    params$weights.match <- Wbar(donors, treated, treatment, m)
    params$objval.sc <- tparams$objval
    return(params)
  }
