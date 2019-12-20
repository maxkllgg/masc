#' Solve for optimal model averaging parameter and associated weights and errors conditional on matching estimator.
#'
#' Implements the matching and synthetic control (masc) estimator of Kellogg, Mogstad,
#'  Pouliot, and Torgovitsky (2019), \emph{conditional on a given matching estimator} characterized by m.
#'  \link{masc} loops over evaluations of this function for each candidate matching estimator,
#'  and selects the one which minimizes cross-validation error.
#'
#' @param data A \code{list} containing three named elements: \describe{
#' \item{donors:}{A \eqn{TxN} matrix of outcome paths for untreated units, each column being a control unit.}
#' \item{treated:}{A \eqn{Tx1} matrix of outcome paths for the treated unit.}
#' \item{treatment:}{An integer. The period T' in which treatment occurs (\eqn{T'<T}).}
#' }
#'
#' @param tune_pars A \code{list} containing 3 elements. You must specify the first, and you must specify
#' exctly one of the remaining two elements: \describe{
#' \item{m:}{ a vector of integers representing the set of candidate nearest neighbor estimators.}
#' \item{min_preperiods:}{an integer. The smallest number of estimation periods allowed in a fold used for cross-validation.}
#' \item{set_f:}{a \code{list} containing a single element, a vector of integers. Identifies the set of folds used
#'  for cross-validation. Each integer identifies a fold by the last time period used in estimation.}
#' }
#'
#'
#' @return returns a list containing five objects:
#' \describe{
#' \item{phi:}{selected value for the model averaging parameter (1 is pure matching, 0 pure synthetic control).}
#' \item{m:}{selected matching (nearest neighbor) estimator.}
#' \item{weights:}{The vector length N containing weights placed on each control unit.}
#' \item{pred.error:}{The vector of treatment effects implied by the masc counterfactual, for periods T' to T.}
#' \item{cv.error:}{The cross-validation error value returned by the optimal set of tuning parameter values.}
#' }
#'
#' @family masc functions
#' @references Kellogg, M., M. Mogstad, G. Pouliot, and A. Torgovitsky. Combining Matching and Synthetic Control to Trade
#'  off Biases from Extrapolation and Interpolation. Working Paper, 2019.
#' @export
#' @import data.table
#' @import plyr
cv_masc <-
  function(data,
           tune_pars, nogurobi=FALSE) {
    estimator <- solve_masc
    min_preperiods = tune_pars$min_preperiods
    set_f = tune_pars$set_f
    if ((is.na(min_preperiods) & any(is.na(set_f))) |
        (!is.na(min_preperiods) & !any(is.na(set_f))))
      stop("You must specify precisely one of min_preperiods and set_f.")
    if (!is.na(min_preperiods)) set_f <- min_preperiods:(data$treatment - 2)


    treat <- data$treatment
    set_f <- set_f[which(data$treatment - set_f - 1 >= 1)]
    ####################################
    finalweights<-estimator(data=list(treatment=treat,
                                   donors=data$donors,
                                   treated=as.matrix(data$treated)),
                         tune_pars=list(
                           m=tune_pars$m
                         ), nogurobi = nogurobi)

    Y_sc<-NULL
    Y_match<-NULL
    foldweights<-list()
    Y_treated<-as.vector(data$treated[set_f+1,])
    #SOLVE MATCHING AND SC ESTIMATORS FOR EACH FOLD:
    position<-1
    for(k in set_f){
      treat<-k+1
      foldweights[[position]]<-estimator(data=list(treatment=treat,
                                                 donors=data$donors,
                                                 treated=as.matrix(data$treated)),
                                       tune_pars=list(
                                         m=tune_pars$m
                                       ), nogurobi = nogurobi)
      Y_sc[position]<-   data$donors[k+1,]%*%foldweights[[position]]$weights.sc

      Y_match[position]<- data$donors[k+1,]%*%foldweights[[position]]$weights.match

      position=position+1
    }
    #ANALYTIC EXPRESSION FOR PHI:
    phi<-(Y_treated-Y_sc)%*%(Y_match-Y_sc)/((Y_match-Y_sc)%*%(Y_match-Y_sc))
    phi<-max(0,phi)
    phi<-min(1,phi)

    #IMPLIED COUNTERFACTUAL AND ERRORS:
    cv.error<- mean((Y_treated-phi*Y_match-(1-phi)*Y_sc)^2)
    weight<-phi*finalweights$weights.match+(1-phi)*finalweights$weights.sc

    output<-list(phi=phi,m=tune_pars$m,weights=weight)
    output$pred.error<-data$treated[data$treatment:nrow(data$donors),] - data$donors[data$treatment:nrow(data$donors),]%*%weight
    output$cv.error<-cv.error

    return(output)
  }


#' Cross-validated Estimation of the Matching and Synthetic Control Estimator.
#'
#' Implements the matching and synthetic control (masc) estimator of Kellogg, Mogstad,
#'  Pouliot, and Torgovitsky (2019).
#'
#'  The \code{masc} estimator takes a convex combination of a
#'  nearest neighbor estimator and a synthetic control estimator. That combination
#'  is parametrized by a model averaging parameter which takes a value of 1 when
#'  the estimator is equivalent to a pure matching estimator, and 0 when it is equivalent
#'  to a pure synthetic control estimator.
#'   his includes selecting the nearest neighbor estimator and model
#'  averaging parameter by a rolling-origin cross-validation procedure. Computationally,
#'  we minimize the cross-validation criterion in two steps. First, for each candidate
#'  nearest neighbor estimator, we solve for the model averaging parameter using an analytic
#'  expression (see Equation 15 of the working paper). Then, we select the candidate nearest
#'  neighbor estimator which produces the smallest cross-validation criterion value.
#'
#'  This implementation by default uses \code{gurobi} interfaced with R to solve for the synthetic control estimator.
#'  Gurobi and its associated R package are available on the gurobi website:
#'  \url{https://cran.r-project.org/web/packages/prioritizr/vignettes/gurobi_installation.html}
#'
#'  It is recommended to use this implementation of the \code{masc} estimator. However, the function may
#'  alternatively use an implementation based on \link[LowRankQP]{LowRankQP} for convenience.
#'
#' @param data A \code{list} containing three named elements: \describe{
#' \item{donors:}{A \eqn{TxN} matrix of outcome paths for untreated units, each column being a control unit.}
#' \item{treated:}{A \eqn{Tx1} matrix of outcome paths for the treated unit.}
#' \item{treatment:}{An integer. The period T' in which treatment begins (\eqn{T'<=T}).}
#' }
#'
#' Note that currently this estimator is designed so that \code{donors} and \code{treated} should only contain
#' the time series values for a single outcome for each unit. It does not allow for including other covariates.
#'
#'
#' @param tune_pars_list A \code{list} containing 3 elements. You must specify the first, and you must specify
#' exctly one of the remaining two elements: \describe{
#' \item{m:}{ a vector of integers representing the set of candidate nearest neighbor estimators.}
#' \item{min_preperiods:}{an integer. The smallest number of estimation periods allowed in a fold used for cross-validation.}
#' \item{set_f:}{a \code{list} containing a single element, a vector of integers. Identifies the set of folds used
#'  for cross-validation. Each integer identifies a fold by the last time period used in estimation.}
#' }
#'
#' @param nogurobi A logical value. If true, uses \link[LowRankQP]{LowRankQP} to solve the synthetic control estimator,
#' rather than \code{gurobi}.
#'
#' @param alloutput A logical value. If true, output includes a list \code{all.results} containing
#' full set of output associated with each candidate nearest neighbor estimator.
#'
#' @return By default, returns a list containing five objects:
#' \describe{
#' \item{phi:}{selected value for the model averaging parameter (1 is pure matching, 0 pure synthetic control).}
#' \item{m:}{selected matching (nearest neighbor) estimator.}
#' \item{weights:}{The vector length N containing weights placed on each control unit.}
#' \item{pred.error:}{The vector of treatment effects implied by the masc counterfactual, for periods T' to T.}
#' \item{cv.error:}{The cross-validation error value returned by the optimal set of tuning parameter values.}
#' }
#' Additionally, returns the following object if \code{alloutput=TRUE}:
#' \describe{
#' \item{all.results:}{A list containing the five above output components for each candidate matching estimator. Values for
#' each candidate are appended as columns.}
#' }
#'
#' @family masc functions
#' @references Kellogg, M., M. Mogstad, G. Pouliot, and A. Torgovitsky. Combining Matching and Synthetic Control to Trade
#'  off Biases from Extrapolation and Interpolation. Working Paper, 2019.
#'
#' @examples
#' ##First example: third introductory exercise of Kellogg et al. (2019),
#' ## Figures 3 and 4:
#'
#' data3 <- data.frame(t = 1:10, treated = 1:10,
#' control1 = c(15,12,14,22,21,28,29,30,29,31),
#' control2 = c(10,5,10,12,16,20,22,21,23,21),
#' control3 = c(-7,-11,-7,-3,-4,-9,-7,-3,-7,-11),
#' control4 = c(-15,-13,-14,-16,-15,-12,-13,-13,-14,-14))
#'
#' #controls are above stored as differences from treated unit. Translating into levels:
#' data3$control1 <- data3$control1 + data3$treated
#' data3$control2 <- data3$control2 + data3$treated
#' data3$control3 <- data3$control3 + data3$treated
#' data3$control4 <- data3$control4 + data3$treated
#'
#'#defining treatment period:
#'treatperiod <- 6
#'
#' data<-list(treated = as.matrix(data3$treated),
#'            donors = as.matrix(data3[,-c(1,2)]),
#'             treatment=treatperiod)
#'
#'result<-masc(data=data, tune_pars_list=list(m=1:3,
#'                                            min_preperiods=3,
#'                                            set_f=NA))
#'
#'#an equivalent specification:
#'result<-masc(data=data, tune_pars_list=list(m=1:3,
#'                                            min_preperiods=NA,
#'                                            set_f=list(3:4)))
#'
#'
#' #Weights selected, for controls 1 through 4 respectively:
#' print(round(result$weights,2))
#'
#' ##Second example: Terrorism in the Basque Region, from
#' ##Abadie and Gardeazabal (2003).
#'
#' #First, load the Synth package, which includes the dataset:
#' if (requireNamespace("Synth",quietly=TRUE) & requireNamespace("data.table",quietly=TRUE)){
#' library(Synth)
#' library(data.table)
#' data(basque)
#' basque<-as.data.table(basque)
#' basque <- basque[regionno!=1,]
#' basque[,regionname:= gsub(" (.*)","",regionname)]
#' #Grabbing region names:
#' names<- c(unique(basque[regionno==17,regionname]),unique(basque[regionno!=17,regionname]))
#' basque <- cbind(basque[regionno==17,gdpcap],
#'                                             t(reshape(basque[regionno!=17,.(regionno,year,gdpcap)],
#'                                              idvar='regionno', timevar='year',direction='wide')[,-"regionno",with=FALSE]))
#'
#' data <- list(treated = as.matrix(basque[,1]),
#'              donors = as.matrix(basque[,-1]),
#'              treatment = 16)
#'
#' result <- masc(data=data, tune_pars_list=list(m=1:10,
#'                                              min_preperiods=8,
#'                                              set_f=NA))
#' names(result$weights)<-names[-1]
#'
#' #weights on control units:
#' print(round(result$weights,3))
#'
#' #Treatment effects of terrorism on GDP per capita
#' #in thousands of 1986 US dollars, over 1970-1975:
#' #(first 6 years of treatment)
#' print(result$pred.error[1:6,])
#'}
#'
#' @export
masc <-
  function(data,
           tune_pars_list = list(m = NA,
                                 min_preperiods=NA,
                                 set_f=NA),
           nogurobi = FALSE,
           alloutput = FALSE) {

    if ((!is.na(tune_pars_list$min_preperiods) &
         !any(is.na(tune_pars_list$set_f))) |
        (is.na(tune_pars_list$min_preperiods) &
         any(is.na(tune_pars_list$set_f))))
      stop("You must specify precisely one of min_preperiods and set_f. The other must be NA.")
      tune_pars_joint <- list()
      position = 1

      for (d in 1:length(tune_pars_list$m)) {
          for (g in 1:length(tune_pars_list$min_preperiods)) {
            for (j in 1:length(tune_pars_list$set_f)) {
              tune_pars_joint[[position]] <- list(
                m = tune_pars_list$m[d],
                min_preperiods = tune_pars_list$min_preperiods[g],
                set_f = tune_pars_list$set_f[[j]]
              )
              position = position + 1

            }
          }
      }
    allresults <-
      lapply(tune_pars_joint, function(x)
        cv_masc(data = data,
          tune_pars = x, nogurobi=nogurobi))
    #allresults is a list of unnamed things, each with named list (weights, pred error, fold stuff, tune pars, cv errors)
    output <-
      allresults[[which.min(lapply(allresults, function(x)
        x$cv.error))]]
    if(alloutput==TRUE){
    output$all.results <- list()

    output$all.results$phi<-sapply(allresults, function(x)
      x$phi)
    names(output$all.results$phi)<-paste0(tune_pars_list$m,"-NN")

    output$all.results$weights<-sapply(allresults, function(x)
      x$weights)
    colnames(output$all.results$weights)<-paste0(tune_pars_list$m,"-NN")

    output$all.results$cv.error <-
      sapply(allresults, function(x)
        x$cv.error)
    names(output$all.results$cv.error)<-paste0(tune_pars_list$m,"-NN")

      output$all.results$pred.error <-
        sapply(allresults, function(x)
          x$pred.error)
      colnames(output$all.results$pred.error)<-paste0(tune_pars_list$m,"-NN")
    }
    return(output)
  }
