#' explain.barp
#'
#' This function is implemented to calculate the contribution of each variable
#' in the BARP(Bayesian Additive Regression Tree with post-stratification) model
#' using the permutation method.
#'
#' @param object A BARP model (Bayesian Additive Regression Tree) estimated
#' using the barps function, a modified version of the barp function from the BARP library with a fixed seed.
#' @param feature_names The name of the variable for which you want to check the contribution.
#' The default value is set to NULL, which means the contribution of all variables in X will be calculated.
#' @param X The dataset containing all independent variables used as input when estimating the BART model.
#' @param nsim The number of MCMC simulation iterations, which is fixed at 1 by default
#' in the case of the BARP model.
#' @param pred_wrapper A function used to estimate the predicted values of the model.
#' @param newdata New data containing the variables included in the model.
#' This is used when checking the contribution of newly input data using the model.
#' The default value is set to NULL, meaning that the input X data,
#' i.e., the data used for model estimation, will be used by default.
#' @param parallel The default value is set to FALSE,
#' but it can be changed to TRUE for parallel computation.
#' @return phis A list containing the Shapley values for each variable.
#' @return newdata The data used to check the contribution of variables.
#' If a variable has two categories, it is dummy-coded, and if it has three or more categories,
#'  categorical variables are one-hot encoded.
#' @return fnull The expected value of the model's predictions.
#' @return fx The prediction value for each observation.
#' @return factor_names The name of the categorical variable.
#' If the data contains only continuous or dummy variables, it is set to NULL.
#' @references Štrumbelj, Erik, and Igor Kononenko. "Explaining prediction models and individual
#' predictions with feature contributions." Knowledge and information systems 41 (2014): 647-665.
#' Bisbee, James. "Barp: Improving mister p using bayesian additive regression trees." American Political Science Review 113.4 (2019): 1060-1065.
#' @export

explain.barp <- function(object, feature_names = NULL, X = NULL, nsim = 1, pred_wrapper = NULL,
                         census = NULL, geo.unit = NULL,   parallel = FALSE,   ...) {


   if (is.null(census)) {
    stop("Please enter your census data.", call. = FALSE)
  }

  # Only nsim = 1
  if (nsim > 1) stop ("It stops because nsim > 1.", call. = FALSE)

  # prediction function check
  if (is.null(pred_wrapper)) {
    stop("Prediction function required for approximate Shapley values.", call. = FALSE)
  }

  if (is.null(X)) {
    X <- object $ trees $ X
  }
  # Convert to data frame
  if (inherits(X, what = "tbl_df")) {
    X <- as.data.frame(X)
  }

    # expacted value
    fx <-  object$barp.dat$SL.bartMachine_1_All

    # predicted value
    fnull <-  object $ pred.opn [,1:2]
    names (  fnull ) [2] = "f_null"

    # variable names
    feature_names = names(object $ trees $ X)

    # Variable sorting
    censusdata <- as.data.frame(census[,  feature_names])
    X <- X [,  feature_names]

    # Categorical variable processing
    newdata <- as.data.frame(bartMachine:::pre_process_new_data(censusdata,object $ trees))

  # Set up the 'foreach' "do" operator
  `%.do%` <- if (isTRUE(parallel)) `%dopar%` else `%do%`

  set.seed(object $ setSeed)
  # Compute approximate Shapley values
  phis_temp <-  foreach(i = feature_names  ) %.do% {
    fastshap :::explain_column ( object$trees, X = X, column = i,
                              pred_wrapper = pred_wrapper,newdata = censusdata)
  }
  names( phis_temp) <- feature_names


  # Deal with other NULL arguments
  geo_check = stringr::str_detect(colnames(newdata),paste0("^", geo.unit))
  featurenames <- colnames(newdata) [-which(geo_check)]


  phis  <-  foreach(i = 1:length(featurenames) ) %.do% {
    ind = featurenames[i]
    phi_idx = which( str_detect( ind , feature_names ))

    temp = matrix (0, nrow = dim (newdata)[1],ncol= object$trees$num_iterations_after_burn_in )
    temp [which(newdata[,ind]==1),]  = phis_temp[[phi_idx]] [which(newdata[,ind]==1),]
    temp
  }
  names( phis) =   featurenames


  out = list (phis = phis, newdata = newdata,   fnull =  fnull, fx = fx ,factor_names = featurenames)

  class( out ) = "explainbarp"

  return( out )
}


