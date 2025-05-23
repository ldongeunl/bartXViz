#' explain.bartMachine
#'
#' This function is used to calculate the contribution of each variable
#' in the Bayesian Additive Regression Trees (BART) model using permutation.
#'
#' @param object A BART model (Bayesian Additive Regression Tree) estimated
#' using the bartMachine function from the bartMachine library.
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

explain.bartMachine <- function(object, feature_names = NULL,  X = NULL,
                                nsim = 1, pred_wrapper = NULL,
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


  if (is.null(newdata)) {  # explain all rows of background data set
    newdata <- X
  }

  if (is.null(pred_wrapper)) {
    stop("Prediction function required for approximate Shapley values. Please ",
         "supply one via the `pred_wrapper` argument", call. = FALSE)
  }


  fx <-   predict(object, new_data = newdata)
  # baseline value (i.e., avg training prediction)
  fnull <- mean(as.numeric(object $ y_hat_train))

  # Deal with other NULL arguments
  if (is.null(feature_names)) {
    feature_names = colnames(X)
  }



  # Set up the 'foreach' "do" operator
  `%.do%` <- if (isTRUE(parallel)) `%dopar%` else `%do%`
  # Compute approximate Shapley values # ,  ...
  phis_temp <- foreach(i = feature_names ) %.do% {
    fastshap:::explain_column(object, X = X, column = i, pred_wrapper = pred_wrapper,
                              newdata = newdata)
    # number of post by obs = n matrix , list = number of variable
  }
  names(phis_temp) <- feature_names

  temp_new =  as.data.frame(bartMachine:::pre_process_new_data(newdata,object))
  featurenames = names(  temp_new )

  tmp_var = str_replace ( featurenames, paste0("_",'[:digit:]'),"")
  tmp_var = data.frame (  tmp_var )%>% group_by(tmp_var)%>% mutate (n=n())%>%distinct(.keep_all = T)


  idx = NULL
  for (j in tmp_var$tmp_var[tmp_var$n == 2]){
    idx = c(idx,min(which(str_detect( featurenames, paste0("^", j) ))))
  }

  if(is.null (idx) == F){
    temp_new =   temp_new [,-idx]
    featurenames = names(  temp_new )
  }



  if ( (sum(sapply(newdata, is.factor)) > 0  | sum(sapply(newdata, is.character)) > 0) & (
       length(which(str_detect( featurenames, '[0-9]$'   ))) >= 2 )) {


    phis  <-  foreach(i = 1:length(  featurenames ) ) %.do% {

      ind = names(temp_new)[i]
      phi_idx = which( str_detect( ind , feature_names ))


      if( tmp_var$n[tmp_var$tmp_var ==   feature_names[i]] >  2 ){
        value_factor =  tail(unlist(strsplit( ind , split = "_")),1)
        temp_var = str_replace (  ind , paste0("_",value_factor),"")
      }


      temp = matrix (0, nrow = dim (newdata)[1],
                     ncol= object$num_iterations_after_burn_in )

      if(length(which(str_detect( featurenames, paste0("^",temp_var)  ))) > 2) {

        temp [which(newdata[,phi_idx ]==value_factor),]  =
          phis_temp[[phi_idx]] [which(newdata[,phi_idx]==  value_factor),]
        temp

      } else  if(length(which(str_detect( featurenames, paste0("^",temp_var)  ))) == 2) {

        temp [which(newdata[,phi_idx ]==1),]  =
          phis_temp[[phi_idx]] [which(newdata[,phi_idx]==  1),]
        temp

      } else {
        temp  = phis_temp[[i]]
        temp
      }
    }
    names( phis) =   featurenames


    if (isTRUE(is.na(tmp_var$tmp_var[tmp_var$n == 2]))  | is.null(tmp_var$tmp_var[tmp_var$n == 2])){
      names( phis) = names(temp_new)

       idx = which( str_detect(  names ( temp_new), paste0("^",temp_var) ))[1]
      newdata =  temp_new[,-idx]
      factor_names = names(X) [which((names(X) %in% names(newdata)) ==FALSE)]

    } else if ( length( tmp_var$tmp_var[tmp_var$n == 2] ) > 1  ){

      names( phis) = names(X)

      factor_names = names(X) [which((names(X) %in% names( temp_new)) ==FALSE)]
      names(temp_new) = str_replace ( names( temp_new), paste0("_",'[:digit:]',"$"),"")

    }
    }else if ( sum(sapply(newdata, is.factor)) == 0 ){
    phis = phis_temp
    factor_names = NULL
  }


  out = list (phis = phis, newdata = newdata,   fnull =  fnull,  fx = fx, factor_names = factor_names)

  class(out) = "explainbartMachine"
  return (out)


}
