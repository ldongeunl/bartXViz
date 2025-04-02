#' explain.default
#'
#' This function is used to calculate the contribution of each variable
#' in model using permutation.
#'
#' @param object model
#' @param feature_names The name of the variable for which you want to check the contribution.
#' The default value is set to NULL, which means the contribution of all variables in X will be calculated.
#' @param X The dataset containing all independent variables used as input when estimating the BART model.
#' @param nsim The number of MCMC simulation iterations, which is fixed at 1 by default.
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
#' @export
explain.default  <- function(object, feature_names = NULL, X = NULL,
                             nsim = 1, pred_wrapper = NULL,
                             newdata = NULL, parallel = FALSE,...) {

  # Compute baseline/average training prediction (fnull) and predictions
  # associated with each explanation (fx); if `adjust = FALSE`, then the
  # baseline is not needed and defaults to zero.
  if (is.null(X)) {
    stop("Training features required for approximate Shapley values. Please ",
         "supply them via the `X` argument; see `?fastshap::explain` for ",
         "details.", call. = FALSE)
  }
  if (inherits(X, what = "tbl_df")) {
    X <- as.data.frame(X)
  }
  if (is.null(newdata)) {
    newdata <- as.data.frame( X )  # explain all rows of background data set
  }


  fx <- pred_wrapper(object, newdata = newdata)
  fnull <-  mean(pred_wrapper(object, newdata = X))


  # Deal with other NULL arguments
  if (is.null(feature_names)) {
    feature_names = colnames(X)
  }
  if (is.null(pred_wrapper)) {
    stop("Prediction function required for approximate Shapley values. Please ",
         "supply one via the `pred_wrapper` argument; see ",
         "`?fastshap::explain` for details.", call. = FALSE)
  }

  # Set up the 'foreach' "do" operator
  `%.do%` <- if (isTRUE(parallel)) `%dopar%` else `%do%`

   # decoded data ----------------- # mltools: one_hot
   temp_new = label_data  ( newdata) $  decoded_data
   temp_X = label_data  ( X )$  decoded_data

   temp_feature = names( temp_X )

  # Compute approximate Shapley values # ,  ...
   phis_temp <- foreach(i = temp_feature ) %.do% {
    replicate(nsim, {  # replace later with vapply()
      explain_column_default (object, X = temp_X , column = i,
                              pred_wrapper = pred_wrapper, newdata =  temp_new)
    })
    # number of post by obs = n matrix , list = number of variable
  }
  names( phis_temp) <- temp_feature


  tmp_var = label_data  ( newdata) $ factor_check

  idx = NULL
  if ( sum (tmp_var $ n >= 2) > 1 ){

  for (j in  unique(  tmp_var$var_tmp [  tmp_var $ uniq_len ==2 &   tmp_var $ n >=2]) ){
    idx = c(idx ,which(str_detect( feature_names, paste0("^", j) )))
  }
 }

  phis  <-  foreach(i = 1:length( feature_names ) ) %.do% {

    temp = matrix (0, nrow = dim (newdata)[1], ncol = nsim )

    if ( i %in% idx ) {
      temp [which(newdata[, i]==1),]  =
        phis_temp[[tmp_var$var_tmp[i]]] [which(newdata[,i]==  1),]
      temp

    } else {
      temp =  phis_temp[[feature_names[i]]]

    }

  }
  names( phis) =   feature_names

  factor_names = NULL
  factor_names = names(X) [which((names(X) %in% names( temp_new)) ==FALSE)]

  out = list (phis = phis, newdata = newdata,  fnull =  fnull, fx = fx, factor_names = factor_names )

  class(out) = "explain"
  return (out)

}

#' @rdname explain.default
#'
#' @export
explain.lm <- function(object, feature_names = NULL, X, nsim = 1, pred_wrapper,
                       newdata = NULL, exact = FALSE, baseline = NULL,
                       parallel = FALSE, ...) {
  if (isTRUE(exact)) {  # use Linear SHAP
    phis <- if (is.null(newdata)) {
      stats::predict(object, type = "terms", ...)
    } else {
      stats::predict(object, newdata = newdata, type = "terms", ...)
    }
    baseline <- attr(phis, which = "constant")  # mean response for all training data

    explain (object, feature_names = feature_names, X = X, nsim = nsim,
              pred_wrapper = pred_wrapper, newdata = newdata, parallel = parallel, ...)
  }
}

#' @rdname explain.default
#'
#' @export

explain.xgb.Booster <- function(object, feature_names = NULL, X = NULL, nsim = 1,
                                pred_wrapper, newdata = NULL, exact = FALSE,   parallel = FALSE, ...) {
  if (isTRUE(exact)) {  # use Tree SHAP
    if (is.null(X) && is.null(newdata)) {
      stop("Must supply `X` or `newdata` argument (but not both).",
           call. = FALSE)
    }
    X <- if (is.null(X)) newdata else X
    phis <- stats::predict(object, newdata = X, predcontrib = TRUE,
                           approxcontrib = FALSE, ...)

    pred_y <- predict(object, X)

      return(structure(list(
        "phis" = phis,
        "newdata" = newdata[, feature_names, drop = FALSE],
        "fnull" =  unique(phis[, "BIAS"]),
        "fx" =  pred_y,
        "factor_names" =NULL
      ), class = "explain"))

  } else {
    explain (object, feature_names = feature_names, X = X, nsim = nsim,
                    pred_wrapper = pred_wrapper, newdata = newdata,
                    parallel = parallel, ...)
  }
}

#' @rdname explain.default
#'
#' @export


explain.lgb.Booster <- function(object, feature_names = NULL, X = NULL, nsim = 1,
                                pred_wrapper, newdata = NULL,
                                exact = FALSE,  parallel = FALSE, ...) {
  if (isTRUE(exact)) {  # use Tree SHAP
    if (is.null(X) && is.null(newdata)) {
      stop("Must supply `X` or `newdata` argument (but not both).",
           call. = FALSE)
    }
    X <- if (is.null(X)) newdata else X

    # Adapt LightGBM predict() interface
    if (utils::packageVersion("lightgbm") > package_version("3.3.2")) {
      phis <- stats::predict(object, X, type = "contrib", ...)
    } else {
      phis <- stats::predict(object, X, predcontrib = TRUE, ...)
    }

    colnames(phis) <- c(colnames(X), "BIAS")

    pred_y = predict(object, X)

      return(structure(list(
        "phis" = phis,
        "newdata" = newdata[, feature_names, drop = FALSE],
        "fnull" =  unique(phis[, "BIAS"]),
        "fx" =  pred_y,
        "factor_names" = NULL
      ), class = "explain"))

  } else {
    explain (object, feature_names = feature_names, X = X, nsim = nsim,
                    pred_wrapper = pred_wrapper, newdata = newdata,
                    parallel = parallel, ...)
  }
}
