#' Report cross-validation error and Prediction Error, Moving between Matching and Synthetic Controls
#'
#'Reports treatment effect and cross-validation errors for estimators
#' of the form of the matching and synthetic control (masc) estimator of Kellogg, Mogstad,
#'  Pouliot, and Torgovitsky (2019). For a set of masc-type estimators defined by a synthetic control
#'  estimator, a matching estimator (m) and a weight (phi), this function returns output associated with
#'  the masc estimator constructed by placing a weight of phi on the matching estimator and (1-phi) on
#'  the synthetic control estimator.
#'
#'
#' @param treated A \eqn{Tx1} matrix of outcomes for the treated unit.
#'
#' @param donors A \eqn{TxN} matrix of outcome paths for untreated units, each column being a control unit.
#'
#' @param treatment An integer. The period T' in which forecasting begins.
#'  If \code{NULL} or T'>T, then we assume all data is pre-treatment.
#'
#'@param sc_est A \code{function} which constructs weights associated with a synthetic control-type estimator. See
#'\link{sc_estimator} for input and output if you'd prefer to substitute your own estimator.
#'
#' @param tune_pars A \code{list} containing 3 elements. You must specify the first, and you may specify
#' only one of the last two elements. The last two elements describe the folds we include in the cross-validation procedure.
#' Each fold \code{f} is denoted by the last period it uses for estimation. That is, fold \code{f} will fit estimators using
#' data from period 1 through period \code{f}, and forecast into period \code{f+1}. \describe{
#' \item{m:}{ a vector of integers. Denotes the set of nearest neighbor estimators from which we are allowed to pick.
#' E.g., \code{tune_pars_list$m=c(1,3,5)} would allow us to pick from 1-NN, 3-NN, or 5-NN.
#' Alternatively, \code{tune_pars_list$m} permits a logical vector. In this case, e.g., \code{tune_pars_list$m=c(FALSE,TRUE,TRUE)}
#' would allow us to pick from 2-NN or 3-NN.
#' If \code{NULL}, we default to allowing all possible nearest neighbor estimators.}
#' \item{min_preperiods:}{an integer. The smallest number of estimation periods allowed in a fold used for cross-validation.
#' We use all folds from fold \code{min_preperiods} up to the latest possible fold \code{treatment-2}.}
#' \item{set_f:}{a \code{list} containing a single element, a vector of integers. Identifies the set of folds used
#'  for cross-validation. As above, each integer identifies a fold by the last time period it uses in estimation.
#'  E.g., set_f=c(7,8,9) would implement cross-validation using fold 7, fold 8, and fold 9.}
#' }
#' If neither \code{min_preperiods} nor \code{set_f} are specified, then we set \code{min_preperiods} to \code{ceil(treatment/2)}.
#' In other words, we pick the first cross-validation fold so that it is estimated on the first half of the pre-period data.
#'
#'
#' @param nogurobi A logical value. If true, uses \link[LowRankQP]{LowRankQP} to solve the synthetic control estimator,
#' rather than \code{gurobi}.
#'
#' @param phivals A vector of real values between 0 and 1. Indexes a weighted
#' average of the synthetic control estimator with a matching estimator, where \code{phival}
#'indicates the weight on the matching estimator (\code{1-phival} being the
#' weight on synthetic controls).
#'
#' @param treatinterval A vector. Indicates the post-treatment periods used when reporting average
#' treatment effects in the column \code{pred} of the output. E.g., \code{treatinterval=1:5} causes the \code{pred}
#' column to report the average treatment effect over the first 5 treatment periods. If \code{NULL}, we default
#' to averaging over all post-treatment periods.
#'
#' @return Returns a \code{data.frame} with each row defined by a value of \code{m} and \code{phi} taken respectively
#' from \code{tune_pars$m} and \code{tune_pars$phivals}. The columns
#' \code{cv.error} and \code{pred} return respectively the cross-validation error and
#' a measure of prediction error (AKA treatment effect) associated with the masc estimator defined by \code{m} and \code{phi}.
#'
#' @family masc functions
#' @references Kellogg, M., M. Mogstad, G. Pouliot, and A. Torgovitsky. Combining Matching and Synthetic Control to Trade
#'  off Biases from Extrapolation and Interpolation. Working Paper, 2019.
#'
#
#'@examples
#'  ##Example: Terrorism in the Basque Region, from
#' ##Abadie and Gardeazabal (2003).
#'
#' #First, load the Synth package, which includes the dataset:
#'
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
#'
#' result <- masc(treated=basque[,1], donors=basque[,-1],treatment=16, tune_pars_list=list(m=1:10,
#'                                              min_preperiods=8))
#' names(result$weights)<-names[-1]
#'
#' #weights on control units:
#' print(round(result$weights,3))
#'
#' #Treatment effects of terrorism on GDP per capita
#' #in thousands of 1986 US dollars, over 1970-1975:
#' #(first 6 years of treatment)
#' print(result$pred.error[1:6,])
#'
#' #Selected tuning parameters?
#' print(paste0("Selected matching estimator: ",result$m_hat))
#' print(paste0("Selected weight on matching: ",result$phi_hat))
#'
#' #Now, examine the shape of A) the CV error (mean square prediciton error in pre-period) and
#' # B) average prediction error (AKA treatment effect) over the first 5 treatment years,
#' #both over values of phi, fixing the matching estimator (moving from matching to synthetic controls)
#' phis<-seq(0,1,length.out=100)
#'phi_table<-masc_by_phi(treated=basque[,1], donors=basque[,-1],treatment=16, tune_pars=list(m=result$m_hat,
#'                                              min_preperiods=8,phis=phis))
#'#Printing CV error and prediction error over values of phi. CV error is clearly lowest at intermediary values of phi,
#'#suggesting an estimator between matching and synthetic controls does best at forecasting. The average medium-run treatment
#'#effect is monotonically increasing as we move away from synthetic control and toward matching.
#'print(phi_table)
#'}
#'
#'
#' @export
#' @import data.table
#' @import plyr
masc_by_phi<-function(treated, donors,treatment=NULL,
                      tune_pars=list(min_preperiods=NULL,set_f=NULL,
                                     m=NULL,phis=seq(from=0,to=1,length.out=100)),sc_est=sc_estimator,
                      treatinterval=NULL){
  phi_table<-NULL
  if(is.logical(tune_pars$m)) tune_pars$m<-which(tune_pars$m)
  if(is.null(tune_pars$m))tune_pars$m<-1:ncol(donors)
for(mpos in 1:length(tune_pars$m)){
  mval <- tune_pars$m[mpos]
  temp_table<-data.frame(phi=tune_pars$phis)
  temp_table$m<-mval
  temp_table$cv<-NA
  temp_table$pred<-NA
  if(all(is.null(treatinterval),!is.null(treatment),treatment<=length(treated))) treatinterval<-1:(length(treated)-treatment+1)
  if(any(treatment>length(treated),is.null(treatment))) treatinterval<-NA
  for(ph in 1:length(tune_pars$phis)){
  tunevals<-tune_pars
  tunevals$m<-mval
phi_result<-cv_masc(treated=treated, donors=donors,treatment=treatment,sc_est=sc_est,tune_pars=tunevals,
        phival=tune_pars$phis[ph])
        temp_table[ph,]$cv<-phi_result$cv.error
        if(!any(is.na(treatinterval))) temp_table[ph,"pred"]<-mean(phi_result$pred.error[treatinterval])
        else temp_table[ph,"pred"]<-NA
}
  phi_table<-rbind(phi_table,temp_table)
}
return(phi_table)
}

#' Solve for optimal model averaging parameter and associated weights and errors conditional on matching estimator.
#'
#' Implements the matching and synthetic control (masc) estimator of Kellogg, Mogstad,
#'  Pouliot, and Torgovitsky (2019), \emph{conditional on a given matching estimator} characterized by m.
#'  \link{masc} loops over evaluations of this function for each candidate matching estimator,
#'  and selects the one which minimizes cross-validation error.
#' @param treated A \eqn{Tx1} matrix of outcomes for the treated unit.
#'
#' @param donors A \eqn{TxN} matrix of outcome paths for untreated units, each column being a control unit.
#'
#' @param treatment An integer. The period T' in which forecasting begins.
#'  If \code{NULL} or T'>T, then we assume all data is pre-treatment.
#'
#'@param sc_est A \code{function} which constructs weights associated with a synthetic control-type estimator. See
#'\link{sc_estimator} for input and output if you'd prefer to substitute your own estimator.
#' @param tune_pars A \code{list} containing 3 elements. You must specify the first, and you may specify
#' only one of the last two elements. The last two elements describe the folds we include in the cross-validation procedure.
#' Each fold \code{f} is denoted by the last period it uses for estimation. That is, fold \code{f} will fit estimators using
#' data from period 1 through period \code{f}, and forecast into period \code{f+1}. \describe{
#' \item{m:}{ an integer representing the  nearest neighbor estimator used.}
#' \item{min_preperiods:}{an integer. The smallest number of estimation periods allowed in a fold used for cross-validation.
#' We use all folds from fold \code{min_preperiods} up to the latest possible fold \code{treatment-2}.}
#' \item{set_f:}{a \code{list} containing a single element, a vector of integers. Identifies the set of folds used
#'  for cross-validation. As above, each integer identifies a fold by the last time period it uses in estimation.
#'  E.g., set_f=c(7,8,9) would implement cross-validation using fold 7, fold 8, and fold 9.}
#' }
#' If neither \code{min_preperiods} nor \code{set_f} are specified, then we set \code{min_preperiods} to \code{ceil(treatment/2)}.
#' In other words, we pick the first cross-validation fold so that it is estimated on the first half of the pre-period data.
#'
#'
#' @param nogurobi A logical value. If true, uses \link[LowRankQP]{LowRankQP} to solve the synthetic control estimator,
#' rather than \code{gurobi}.
#'
#' @param phival A real value between 0 and 1. If specified, hard-codes the masc estimator to take the specified weighted
#' average of matching and synthetic controls, where \code{phival} indicates the weight on matching (\code{1-phival} being the
#' weight on synthetic controls).
#'
#' @return returns a list containing five objects:
#' \describe{
#' \item{phi_hat:}{selected value for the model averaging parameter (1 is pure matching, 0 pure synthetic control).}
#' \item{m_hat:}{selected matching estimator (number of nearest neighbor).}
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
  function(treated,donors,treatment=NULL, sc_est=sc_estimator,
           tune_pars=list(min_preperiods=NULL,set_f=NULL,m=NULL), nogurobi=FALSE, phival=NULL) {
    estimator <- solve_masc
    treated<-as.matrix(treated)
    if(is.null(tune_pars$min_preperiods)) tune_pars$min_preperiods<-NA
    if(is.null(tune_pars$set_f)) tune_pars$set_f<-NA
    if(is.null(tune_pars$m)) stop("Must specify a matching estimator in tune_pars")
    min_preperiods = tune_pars$min_preperiods
    set_f = tune_pars$set_f
    if (!is.na(min_preperiods) & !any(is.na(set_f)))
      stop("One of min_preperiods and set_f must be non-missing.")

    if (is.na(min_preperiods) & any(is.na(set_f))) min_preperiods<-ceil(treatment/2)
    if (!is.na(min_preperiods)) set_f <- min_preperiods:(treatment - 2)

  if(is.null(treatment)) treatment<-nrow(treated)+1
    treat <- treatment
    set_f <- set_f[which(treatment - set_f - 1 >= 1)]
    ####################################
    finalweights<-estimator(treatment=treat,
                                   donors=donors,
                                   treated=as.matrix(treated),
                            sc_est=sc_est,
                         tune_pars=list(
                           m=tune_pars$m
                         ), nogurobi = nogurobi)

    Y_sc<-NULL
    Y_match<-NULL
    foldweights<-list()
    Y_treated<-as.vector(treated[set_f+1,])
    #SOLVE MATCHING AND SC ESTIMATORS FOR EACH FOLD:
    position<-1
    for(k in set_f){
      treat<-k+1
      foldweights[[position]]<-estimator(treatment=treat,
                                                 donors=donors,
                                                 treated=as.matrix(treated),
                                         sc_est=sc_est,
                                       tune_pars=list(
                                         m=tune_pars$m
                                       ), nogurobi = nogurobi)
      Y_sc[position]<-   donors[k+1,]%*%foldweights[[position]]$weights.sc

      Y_match[position]<- donors[k+1,]%*%foldweights[[position]]$weights.match

      position=position+1
    }
    #ANALYTIC EXPRESSION FOR PHI:
    if(is.null(phival)){
    phi<-(Y_treated-Y_sc)%*%(Y_match-Y_sc)/((Y_match-Y_sc)%*%(Y_match-Y_sc))
    phi<-max(0,phi)
    phi<-min(1,phi)
    }
    else phi = phival
    #IMPLIED COUNTERFACTUAL AND ERRORS:
    cv.error<- mean((Y_treated-phi*Y_match-(1-phi)*Y_sc)^2)
    weight<-phi*finalweights$weights.match+(1-phi)*finalweights$weights.sc

    output<-list(phi_hat=phi,m_hat=tune_pars$m,weights=weight)
    if(treatment <= nrow(donors)) output$pred.error<-treated[treatment:nrow(donors),] - donors[treatment:nrow(donors),]%*%weight
    else output$pred.error<-NA
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
#'   This function selects the nearest neighbor estimator and model
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
#' @param treated A \eqn{Tx1} matrix of outcomes for the treated unit.
#'
#' @param donors A \eqn{TxN} matrix of outcome paths for untreated units, each column being a control unit.
#'
#' Note that currently this estimator is designed so that \code{donors} and \code{treated} should only contain
#' the time series values for a single outcome for each unit. It does not allow for including other covariates.

#' @param treatment An integer. The period T' in which forecasting begins. If \code{NULL} or T'>T, then we
#'  assume all data is pre-treatment.
#'
#'
#'@param sc_est A \code{function} which constructs weights associated with a synthetic control-type estimator. See
#'\link{sc_estimator} for input and output if you'd prefer to substitute your own estimator.
#'
#' @param tune_pars_list A \code{list} containing 3 elements. You may specify
#' only one of \code{min_preperiods} and \code{set_f}. Those elements describe the folds we include in the cross-validation procedure.
#' Each fold \code{f} is indexed by the last period it uses for estimation. That is, fold \code{f} will fit estimators using
#' data from period 1 through period \code{f}, and forecast into period \code{f+1}. \describe{
#' \item{m:}{ a vector of integers. Denotes the set of nearest neighbor estimators from which we are allowed to pick.
#' E.g., \code{tune_pars_list$m=c(1,3,5)} would allow us to pick from 1-NN, 3-NN, or 5-NN.
#' Alternatively, \code{tune_pars_list$m} permits a logical vector. In this case, e.g., \code{tune_pars_list$m=c(FALSE,TRUE,TRUE)}
#' would allow us to pick from 2-NN or 3-NN.
#' If \code{NULL}, we default to allowing all possible nearest neighbor estimators.}
#' \item{min_preperiods:}{an integer. The smallest number of estimation periods allowed in a fold used for cross-validation.
#' We use all folds from fold \code{min_preperiods} up to the latest possible fold \code{treatment-2}.}
#' \item{set_f:}{a \code{list} containing a single element, a vector of integers. Identifies the set of folds used
#'  for cross-validation. As above, each integer identifies a fold by the last time period it uses in estimation.
#'  E.g., set_f=c(7,8,9) would implement cross-validation using fold 7, fold 8, and fold 9.}
#' }
#' If neither \code{min_preperiods} nor \code{set_f} are specified, then we set \code{min_preperiods} to \code{ceil(treatment/2)}.
#' In other words, we pick the first cross-validation fold so that it is estimated on the first half of the pre-period data.
#'
#' @param nogurobi A logical value. If true, uses \link[LowRankQP]{LowRankQP} to solve the synthetic control estimator,
#' rather than \code{gurobi}.
#'
#' @param alloutput A logical value. If true, output includes a list \code{all.results} containing
#' full set of output associated with each candidate nearest neighbor estimator.
#'
#' @return By default, returns a list containing five objects:
#' \describe{
#' \item{phi_hat:}{selected value for the model averaging parameter (1 is pure matching, 0 pure synthetic control).}
#' \item{m_hat:}{selected matching (nearest neighbor) estimator.}
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
#'result<-masc(treated=data3$treated,donors=data3[,-c(1,2)],treatment=treatperiod,
#'             tune_pars_list=list(m=1:3,min_preperiods=3))
#'
#'#an equivalent specification:
#'result<-masc(treated=data3$treated,donors=data3[,-c(1,2)],treatment=treatperiod,
#'             tune_pars_list=list(m=1:3,set_f=list(3:4)))
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
#'
#' result <- masc(treated=basque[,1], donors=basque[,-1],treatment=16, tune_pars_list=list(m=1:10,
#'                                              min_preperiods=8))
#' names(result$weights)<-names[-1]
#'
#' #weights on control units:
#' print(round(result$weights,3))
#'
#' #Treatment effects of terrorism on GDP per capita
#' #in thousands of 1986 US dollars, over 1970-1975:
#' #(first 6 years of treatment)
#' print(result$pred.error[1:6,])
#'
#' #Selected tuning parameters?
#' print(paste0("Selected matching estimator: ",result$m_hat))
#' print(paste0("Selected weight on matching: ",result$phi_hat))
#'
#' #Now, examine the shape of A) the CV error (mean square prediciton error in pre-period) and
#' # B) average prediction error (AKA treatment effect) over the first 5 treatment years,
#' #both over values of phi, fixing the matching estimator (moving from matching to synthetic controls)
#' phis<-seq(0,1,length.out=100)
#'phi_table<-masc_by_phi(treated=basque[,1], donors=basque[,-1],treatment=16, tune_pars=list(m=result$m_hat,
#'                                              min_preperiods=8),phivals=phis)
#'#Printing CV error and prediction error over values of phi. CV error is clearly lowest at intermediary values of phi,
#'#suggesting an estimator between matching and synthetic controls does best at forecasting. The average medium-run treatment
#'#effect is monotonically increasing as we move away from synthetic control and toward matching.
#'print(phi_table)
#'}
#'
#'
#' @export
masc <-
  function(treated,donors,treatment=NULL,
           sc_est=sc_estimator,
           tune_pars_list = list(),
           nogurobi = FALSE,
           alloutput = FALSE) {

    treated<-as.matrix(treated)
    donors<-as.matrix(donors)
    if(is.null(treatment)) treatment<-nrow(treated)+1
    if ((!is.null(tune_pars_list$min_preperiods) &
         !any(is.null(tune_pars_list$set_f))))
      stop("You may specify precisely one of min_preperiods and set_f in tune_pars_list. The other must be empty (NULL)")

    if(is.null(tune_pars_list$min_preperiods) &
        any(is.null(tune_pars_list$set_f))) tune_pars_list$min_preperiods<-ceil(treatment/2)

    if(is.null(tune_pars_list$set_f)) tune_pars_list$set_f<-NA
    if(is.null(tune_pars_list$min_preperiods)) tune_pars_list$min_preperiods<-NA

    if(is.null(tune_pars_list$m)) tune_pars_list$m<-1:ncol(donors)
    if(is.logical(tune_pars_list$m)) tune_pars_list$m<-which(tune_pars_list$m)
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
        cv_masc(treated = treated,donors=donors,treatment=treatment, sc_est=sc_est,
          tune_pars = x, nogurobi=nogurobi))
    #allresults is a list of unnamed things, each with named list (weights, pred error, fold stuff, tune pars, cv errors)
    output <-
      allresults[[which.min(lapply(allresults, function(x)
        x$cv.error))]]
    if(alloutput==TRUE){
    output$all.results <- list()

    output$all.results$phi<-sapply(allresults, function(x)
      x$phi_hat)
    names(output$all.results$phi_hat)<-paste0(tune_pars_list$m,"-NN")

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
