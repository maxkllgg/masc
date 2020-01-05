##################################
#OPTIMIZATION PROBLEM
#################################

#Nearest Neighbor weight vector
Wbar <- function(donors, treated, treatment, m) {
  dist <-
    apply((donors[1:(treatment - 1), ] - as.vector(treated)[1:(treatment -
                                                                 1)]) ^ 2, MARGIN = 2, FUN = sum)
  W <- rep(0, length(dist))
  W[which(dist %in% sort(dist)[1:m])] <-
    1 / length(which(dist %in% sort(dist)[1:m]))
  return(W)
}


#' solving for standard synthetic control weights
#'
#' Solves for synthetic control weights up to the designated treatment period.
#'
#' @param donors A \eqn{TxN} matrix of outcome paths for untreated units, each column being a control unit.
#' @param treated: A \eqn{Tx1} matrix of outcomes for the treated unit.
#' @param treatment: An integer. The period T' in which forecasting begins (either the true treatment
#' period or the first period after a cross-validation fold).
#'
#' @param nogurobi A logical value. If true, uses \link[LowRankQP]{LowRankQP} to solve the synthetic control estimator,
#' rather than \code{gurobi}.

#' @return A list. The named component \code{weights.sc} contains the vector of synthetic control weights.
#'  Weights are ordered in the same manner as the columns in \code{donors}.
#'  The \code{objval.sc} component contains the objective value (pre-period fit) of the synthetic control.
#' @family masc functions
#' @references Kellogg, M., M. Mogstad, G. Pouliot, and A. Torgovitsky. Combining Matching and Synthetic Control to Trade
#'  off Biases from Extrapolation and Interpolation. Working Paper, 2019.
#' @export
sc_estimator<-function(donors,treated,treatment,nogurobi) {
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

  if(nogurobi==FALSE){
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

    tparams <- plyr::rename(tparams, c('x' = 'pars'))
    if (tparams$status != 'OPTIMAL') {
      print('WARNING: FAILED TO SOLVE PROBLEM')
      continue
    }
    params$weights.sc <- tparams$pars
    params$objval.sc <- tparams$objval
  }
  else {
    tparams <- LowRankQP::LowRankQP(Vmat=2*modelreg$Q,
                                    dvec=modelreg$obj,
                                    Amat=modelreg$A,
                                    bvec=modelreg$rhs,
                                    uvec=modelreg$ub,
                                    method = 'LU')
    params$weights.sc <- as.vector(tparams$alpha)
    params$objval.sc<-mean((treated-donors%*%params$weights.sc)^2)
  }
  return(params)
}


#' solving for matching and synthetic control weights
#'
#' Solves separately for synthetic control weights and matching weights for a given
#' matching estimator, using data up to the treatment period.
#'
#' @param donors A \eqn{TxN} matrix of outcome paths for untreated units, each column being a control unit.
#' @param treated: A \eqn{Tx1} matrix of outcomes for the treated unit.
#' @param treatment: An integer. The period T' in which forecasting begins (either the true treatment
#' period or the first period after a cross-validation fold).
#'
#'#'@param sc_est A \code{function} which constructs weights associated with a synthetic control-type estimator. See
#'\link{sc_estimator} for input and output if you'd prefer to substitute your own estimator.
#'
#' @param tune_pars A \code{list} containing one element:
#' \describe{\item{m:}{an integer identifying the candidate matching (nearest neighbor) estimator.}}
#'
#' @return A list containing the weights associated with the two estimators. Weights are ordered
#' in the same manner as the columns in \code{donors}. The \code{weights.sc} and \code{weights.match}
#' named components contain the vector of synthetic control and nearest neighbor weights respectively.
#'  The \code{objval.sc} component contains the objective value (pre-period fit) of the synthetic control.
#' @family masc functions
#' @references Kellogg, M., M. Mogstad, G. Pouliot, and A. Torgovitsky. Combining Matching and Synthetic Control to Trade
#'  off Biases from Extrapolation and Interpolation. Working Paper, 2019.
#' @export
solve_masc <-
  function(donors,treated,treatment, sc_est,
           tune_pars,nogurobi=FALSE) {
    ####Pulling objects out of lists###
    treatment <- treatment
    donors <- donors[1:(treatment - 1), ]
    treated <- as.matrix(treated[1:(treatment - 1), ])
    m <- tune_pars$m

    params <- sc_est(donors,treated,treatment,nogurobi=nogurobi)
    params$weights.match <- Wbar(donors, treated, treatment, m)
    return(params)
  }
