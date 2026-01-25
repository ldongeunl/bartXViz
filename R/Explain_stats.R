#' Numerical summary of Shapley values from an Explain object
#'
#' @description
#' Computes global numerical summaries of Shapley values using two averaging criteria:
#' observation-based and posterior-sample-based.
#'
#' @param x An object belonging to the \code{Explain} class or its subclasses, containing the Shapley values.
#' @param probs Enter the probability for the quantile interval. Default is \code{0.95}.
#' @param abs Logical. If \code{TRUE}, summarizes absolute Shapley values (importance-style).
#' @param na.rm Logical. Remove NA values in summaries. Default is \code{TRUE}.
#' @param geo.unit (Explainbarp only) Name of the stratification variable used in post-stratification.
#' @param geo.id (Explainbarp only) One value of interest among the values of the stratification variable.
#'
#' @return
#' A named list with two elements:
#' \item{obs}{A data.frame containing observation-based numerical summaries of Shapley values for each variable.}
#' \item{post}{A data.frame containing posterior-sample-based numerical summaries of Shapley values for each variable.}
#' @export
#' @examples
#' \dontrun{
#' ## Friedman data
#' set.seed(2025)
#' n <- 200
#' p <- 5
#' X <- data.frame(matrix(runif(n * p), ncol = p))
#' y <- 10 * sin(pi* X[ ,1] * X[,2]) +20 * (X[,3] -.5)^2 + 10 * X[ ,4] + 5 * X[,5] + rnorm(n)
#' 
#' ## Using the dbarts library
#' model <- dbarts::bart(X,y,keeptrees = TRUE ,  ndpost = 200)
#' 
#' ## prediction wrapper function
#' pfun <- function(object, newdata) {
#'        predict(object , newdata)
#'        }
#'        
#'## Calculate shapley values
#'model_exp <-  Explain  ( model, X = X,  pred_wrapper =  pfun )
#'
#'# Numerical summaries of summarised Shapley values
#'Explain_stats ( model_exp,  probs = 0.95)
#'}
Explain_stats <- function(x, probs = 0.95, abs = TRUE, na.rm = TRUE,
                          geo.unit = NULL, geo.id = NULL) {
  UseMethod("Explain_stats")
}

# ---- internal helper shared by both methods ----
.Explain_stats_core <- function(object, feature_names, newdata_for_mask,
                                probs = 0.95, abs = TRUE, na.rm = TRUE) {
  
  alpha <- (1 - probs) / 2
  q_lo <- alpha
  q_hi <- 1 - alpha
  
  factor_names <- object$factor_names
  length_cate <- character(0)
  
  if (!is.null(factor_names) && length(factor_names) > 0) {
    nd_names <- colnames(newdata_for_mask)
    for (fn in factor_names) {
      hits <- nd_names[grepl(fn, nd_names)]
      if (length(hits) > 1) length_cate <- c(length_cate, hits)
    }
    length_cate <- unique(length_cate)
  }
  
  mask_one_hot <- function(mat, varname) {
    if (!is.null(length_cate) &&
        varname %in% length_cate &&
        varname %in% colnames(newdata_for_mask)) {
      
      temp_i <- which(newdata_for_mask[, varname] != 1)
      if (length(temp_i) > 0) mat[temp_i, ] <- NA
    }
    mat
  }
  
  obs_list  <- vector("list", length(feature_names))
  post_list <- vector("list", length(feature_names))
  names(obs_list)  <- feature_names
  names(post_list) <- feature_names
  
  for (v in feature_names) {
    
    phi <- object$phis[[v]]
    if (is.null(phi)) next
    if (!is.matrix(phi)) {
      stop(sprintf("object$phis[['%s']] must be a matrix (#obs x #post).", v),
           call. = FALSE)
    }
    
    phi <- mask_one_hot(phi, v)
    
    # --------------------
    # obs-based 
    # --------------------
    phi_obs <- if (isTRUE(abs)) base::abs(phi) else phi
    post_means <- apply(phi_obs, 2, mean, na.rm = na.rm)
    
    obs_list[[v]] <- data.frame(
      variable = v,
      Average  = "obs",
      mean     = mean(post_means, na.rm = na.rm),  
      q_lower  = as.numeric(quantile(post_means, q_lo, na.rm = na.rm, names = FALSE)),
      q_upper  = as.numeric(quantile(post_means, q_hi, na.rm = na.rm, names = FALSE)),
      stringsAsFactors = FALSE
    )
    
    # --------------------
    # post-based 
    # --------------------
    obs_use <- apply(phi, 1, function(x) {
      m <- mean(x, na.rm = na.rm)
      if (isTRUE(abs)) base::abs(m) else m
    })
    
    post_list[[v]] <- data.frame(
      variable = v,
      Average  = "post",
      mean     = mean(obs_use, na.rm = na.rm),   
      q_lower  = as.numeric(quantile(obs_use, q_lo, na.rm = na.rm, names = FALSE)),
      q_upper  = as.numeric(quantile(obs_use, q_hi, na.rm = na.rm, names = FALSE)),
      stringsAsFactors = FALSE
    )
  }
  
  obs_keep  <- Filter(Negate(is.null), obs_list)
  post_keep <- Filter(Negate(is.null), post_list)
  
  if (length(obs_keep) == 0 || length(post_keep) == 0) {
    stop("No valid Shapley summaries could be computed.", call. = FALSE)
  }
  
  obs_df  <- do.call(rbind, obs_keep)
  post_df <- do.call(rbind, post_keep)
  
  res <- list(obs = obs_df, post = post_df)
  
  # ranking
  res <- lapply(res, function(df) {
    df <- df[order(-df$mean), , drop = FALSE]
    df$rank <- rank(-df$mean, ties.method = "first")
    rownames(df) <- NULL
    df
  })
  
  res
}

#' @rdname Explain_stats
#' @export
#'

Explain_stats.default <- function(x, probs = 0.95, abs = TRUE, na.rm = TRUE,
                                  geo.unit = NULL, geo.id = NULL) {
  
  if (missing(x)) {
    stop("The object containing Shapley values obtained from 'Explain()' is missing.",
         call. = FALSE)
  }
  if (!is.null(geo.unit) || !is.null(geo.id)) {
    stop("Arguments 'geo.unit' and 'geo.id' are only used for 'Explainbarp' objects.",
         call. = FALSE)
  }
  
  .allowed <- c("Explain", "ExplainBART", "ExplainbartMachine")
  
  if (!any(vapply(.allowed, function(cl) inherits(x, cl), logical(1)))) {
    stop(
      "The input object must be an Explain-class object (e.g., 'Explain', 'ExplainBART', or 'ExplainbartMachine').",
      call. = FALSE
    )
  }
  
  if (is.null(x$phis) || !is.list(x$phis) || length(x$phis) == 0) {
    stop("No Shapley values found: 'x$phis' is missing or empty.", call. = FALSE)
  }
  if (is.null(x$newdata)) {
    stop("No 'newdata' found in the Explain object.", call. = FALSE)
  }
  
  object <- x
  feature_names <- names(object$phis)
  if (is.null(feature_names) || any(feature_names == "")) {
    feature_names <- names(object$newdata)
  }
  
  .Explain_stats_core(
    object = object,
    feature_names = feature_names,
    newdata_for_mask = object$newdata,
    probs = probs, abs = abs, na.rm = na.rm
  )
}

#' @rdname Explain_stats
#' @export
Explain_stats.Explainbarp <- function(x, probs = 0.95, abs = TRUE, na.rm = TRUE,
                                      geo.unit = NULL, geo.id = NULL) {
  
  if (missing(x)) {
    stop("The object containing Shapley values obtained from 'Explain()' is missing.",
         call. = FALSE)
  }
  if (is.null(x$phis) || !is.list(x$phis) || length(x$phis) == 0) {
    stop("No Shapley values found: 'x$phis' is missing or empty.", call. = FALSE)
  }
  if (is.null(x$newdata)) {
    stop("No 'newdata' found in the Explainbarp object.", call. = FALSE)
  }
  if (is.null(geo.unit) || is.null(geo.id)) {
    stop("For 'Explainbarp' objects, both 'geo.unit' and 'geo.id' must be provided.",
         call. = FALSE)
  }
  
  object <- x
  
  # --- identify stratum index  ---
  prefix <- paste0(geo.unit, "_")
  nd_names <- colnames(object$newdata)
  geo_check <- ifelse(startsWith(nd_names, prefix),
                      substring(nd_names, nchar(prefix) + 1L),
                      "")
  
  geo_cols <- which(geo_check == as.character(geo.id))
  if (length(geo_cols) == 0) {
    stop("Could not find a poststratification dummy column matching geo.id in 'newdata'.",
         call. = FALSE)
  }
  geo_cols <- geo_cols[1L]
  
  index <- which(object$newdata[, geo_cols] == 1)
  if (length(index) == 0) {
    stop("No observations found for the specified stratum (geo.unit / geo.id).",
         call. = FALSE)
  }
  
  # feature columns only (exclude geo dummies)
  feature_names <- nd_names[geo_check == ""]
  if (length(feature_names) == 0) {
    stop("No feature columns found after excluding poststratification dummy columns.",
         call. = FALSE)
  }
  
  new_data <- object$newdata[index, geo_check == "", drop = FALSE]
  
  # interval endpoints for stats output
  alpha <- (1 - probs) / 2
  q_lo <- alpha
  q_hi <- 1 - alpha
  
  # --- make temp_phis exactly like your original code ---
  temp_phis <- lapply(feature_names, function(v) {
    phi <- object$phis[[v]]
    if (is.null(phi)) return(NULL)
    phi[index, , drop = FALSE]
  })
  names(temp_phis) <- feature_names
  temp_phis <- Filter(Negate(is.null), temp_phis)
  
  if (length(temp_phis) == 0) {
    stop("No Shapley values available for the selected stratum.", call. = FALSE)
  }
  
  # apply the SAME masking rule as your original code
  for (v in names(temp_phis)) {
    if (!v %in% colnames(new_data)) next
    temp_i <- which(new_data[, v] != 1)
    if (length(temp_i) > 0) temp_phis[[v]][temp_i, ] <- NA
  }
  
  # obs-based
  obs_df <- do.call(rbind, lapply(names(temp_phis), function(v) {
    phi <- temp_phis[[v]]
    phi_use <- if (isTRUE(abs)) base::abs(phi) else phi
    post_means <- apply(phi_use, 2, mean, na.rm = na.rm)
    
    data.frame(
      variable = v,
      Average = "obs",
      mean = mean(post_means, na.rm = na.rm),
      q_lower = as.numeric(quantile(post_means, q_lo, na.rm = na.rm, names = FALSE)),
      q_upper = as.numeric(quantile(post_means, q_hi, na.rm = na.rm, names = FALSE)),
      stringsAsFactors = FALSE
    )
  }))
  
  obs_df <- obs_df[is.finite(obs_df$mean), , drop = FALSE]
  
  # Posterior sample-based
  post_df <- do.call(rbind, lapply(names(temp_phis), function(v) {
    phi <- temp_phis[[v]]
    obs_means <- apply(phi, 1, mean, na.rm = na.rm)
    obs_use <- if (isTRUE(abs)) base::abs(obs_means) else obs_means
    
    data.frame(
      variable = v,
      Average = "post",
      mean = mean(obs_use, na.rm = na.rm),
      q_lower = as.numeric(quantile(obs_use, q_lo, na.rm = na.rm, names = FALSE)),
      q_upper = as.numeric(quantile(obs_use, q_hi, na.rm = na.rm, names = FALSE)),
      stringsAsFactors = FALSE
    )
  }))
  
  post_df <- post_df[is.finite(post_df$mean), , drop = FALSE]
  
  res <- list(obs = obs_df, post = post_df)
  
  # ranking within each table
  res <- lapply(res, function(df) {
    df <- df[order(-df$mean), , drop = FALSE]
    df$rank <- rank(-df$mean, ties.method = "first")
    rownames(df) <- NULL
    df
  })
  
  res
}
