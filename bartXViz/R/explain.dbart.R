#' explain.bart
#'
#' This function is used to calculate the contribution of each variable
#' in the Bayesian Additive Regression Trees (BART) model using permutation.
#'
#' @param object A BART model (Bayesian Additive Regression Tree) estimated
#' using the bart function from the dbarts library.
#' @param feature_names The name of the variable for which you want to check the contribution.
#' The default value is set to NULL, which means the contribution of all variables in X will be calculated.
#' @param X The dataset containing all independent variables used as input when estimating the BART model.
#' @param nsim The number of MCMC simulation iterations, which is fixed at 1 by default
#' in the case of the BART model.
#' @param pred_wrapper A function used to estimate the predicted values of the model.
#' @param newdata New data containing the variables included in the model.
#' This is used when checking the contribution of newly input data using the model.
#' The default value is set to NULL, meaning that the input X data,
#' i.e., the data used for model estimation, will be used by default.
#' @param parallel The default value is set to FALSE,
#' but it can be changed to TRUE for parallel computation.
#' @return phis A list containing the Shapley values for each variable.
#' @return newdata The data used to check the contribution of variables.
#' If a variable has categories, categorical variables are one-hot encoded.
#' @return fnull The expected value of the model's predictions.
#' @return fx The prediction value for each observation.
#' @return factor_names The name of the categorical variable.
#' If the data contains only continuous or dummy variables, it is set to NULL.
#' @references Štrumbelj, Erik, and Igor Kononenko. "Explaining prediction models and individual
#' predictions with feature contributions." Knowledge and information systems 41 (2014): 647-665.
#' @export

explain.bart <- function(object, feature_names = NULL, X = NULL, nsim = 1,  pred_wrapper = NULL,
                         newdata = NULL,   parallel = FALSE, ...) {

  # Only nsim = 1
  if (nsim > 1) stop ("It stops because nsim > 1.",
                      "Because the BART model uses posterior samples,",
                      "it is used by setting nsim=1.", call. = FALSE)


  # Compute baseline/average training prediction (fnull) and predictions
  # associated with each explanation (fx); if `adjust = FALSE`, then the
  # baseline is not needed and defaults to zero.
  if (is.null(X)) {
    stop("Training features required for approximate Shapley values.", call. = FALSE)
  }

  if (inherits(X, what = "tbl_df")) {
    X <- as.data.frame(X)
  }

  if (is.null(pred_wrapper)) {
    stop("Prediction function required for approximate Shapley values. Please ",
         "supply one via the `pred_wrapper` argument", call. = FALSE)
  }

  if (is.null(newdata)) {
       newdata <-  X
  }

  fx  <-   colMeans ( predict(object, newdata =  newdata ))
  fnull <-  mean( colMeans (object$yhat.train)) # baseline value (i.e., avg training prediction)

  # Deal with other NULL arguments
  if (is.null(feature_names)) {
    feature_names = colnames(X)
  }



  # Set up the 'foreach' "do" operator
  `%.do%` <- if (isTRUE(parallel)) `%dopar%` else `%do%`

  # Compute approximate Shapley values # ,  ...
  phis_temp <- foreach(i = feature_names ) %.do% {
     t(fastshap:::explain_column(object, X = X, column = i, pred_wrapper = pred_wrapper,
                                 newdata = newdata))
    # number of post by obs = n matrix , list = number of variable
  }
  names( phis_temp ) <- feature_names


  featurenames = names(as.data.frame(dbarts:::makeModelMatrixFromDataFrame (newdata)))


  if (( sum(sapply(newdata, is.factor)) > 0 | sum(sapply(newdata, is.character)) > 0) &
      length(which(str_detect( featurenames, '[0-9]$' ))) >= 1 ) {

    factor_names = names(X) [which( ( names(X) %in% featurenames  ) ==FALSE)]

    phis  <-  foreach(i = 1:length(featurenames) ) %.do% {
      ind = featurenames[i]
      phi_idx = which( str_detect( ind , feature_names ))
      value_factor =  tail(unlist(strsplit(ind, split = "[.]")),1)

       temp = matrix (0, nrow = dim (newdata)[1],ncol= dim(object $ yhat.train) [1] )
       if(length(which(str_detect(ind , feature_names))) > 2) {

         temp [which(newdata[,phi_idx ]==value_factor),]  = phis_temp[[phi_idx]] [which(newdata[,phi_idx]==  value_factor),]
         temp
       } else {
         temp  = phis_temp[[i]]
         temp
       }

    }
    names( phis) =   featurenames
    newdata = as.data.frame(dbarts:::makeModelMatrixFromDataFrame (newdata))

  } else if ( sum(sapply(newdata, is.factor)) == 0 ){
    phis = phis_temp

    if( sum(featurenames%in%feature_names) != length(featurenames) ){
      factor_names = names(X) [which( ( names(X) %in% featurenames  ) ==FALSE)]
    }else{
      factor_names = NULL
    }

  }


  out = list (phis = phis, newdata = newdata,   fnull =  fnull,  fx = fx, factor_names = factor_names)

  class(out) = "explainBART"
  return (out)

}
