#'*************************************************************************************
#' Implementations for quantile prediction, 
#' main function predict.Cvinereg for class Cvinereg
#' functions are implemented based on vinereg package 
#' major change qcvine function
#' Aug 2019 version by J.Zhou
#'*************************************************************************************
#' @export
#'
#' @importFrom kde1d pkde1d qkde1d
#' @importFrom stats predict
predict.Cvinereg <- function(object, newdata, alpha = 0.5, uscale = FALSE, ...) {
  if (missing(newdata))
    return(fitted.Cvinereg(object, alpha = alpha, uscale = uscale))
  stopifnot(length(alpha) > 0)
  if (any(is.na(alpha)) & inherits(object$model_frame[[1]], "ordered"))
    stop("cannot predict mean for ordered response.")
  
  # check if margins were estimated
  if (!is.null(object$margins[[1]]$u) & (!uscale)) {
    warning("no margins have been estimated, setting uscale = TRUE")
    uscale <- TRUE
  }
  
  # check if all variables in the model are in newdata
  if (is.matrix(newdata))
    newdata <- as.data.frame(newdata)
  missing_vars <- setdiff(colnames(object$model_frame)[-1], colnames(newdata))
  if (length(missing_vars) > 0)
    stop("'newdata' is missing variables ",
         paste(missing_vars, collapse = ", "))
  
  # predict the mean if alpha contains NA
  if (any(is.na(alpha))) {
    alpha <- alpha[!is.na(alpha)]
    preds <- predict_mean_Cvine(object, newdata, uscale)
  } else {
    preds <- NULL
  }
  
  if (length(alpha) > 0) {
    stopifnot(is.numeric(alpha))
    stopifnot(all(alpha > 0) & all(alpha < 1))
    
    # expand factors and make ordered variables numeric
    x <- cctools::expand_as_numeric(newdata[colnames(object$model_frame)[-1]])
    
    # remove unused variables
    selected_vars <- match(object$order, colnames(x))
    u <- x <- x[, selected_vars, drop = FALSE]
    
    # transform to uniform scale
    if (!uscale) {
      for (j in seq_len(ncol(x)))
        u[, j] <- pkde1d(x[, j], object$margins[[selected_vars[j] + 1]])
    }
    
    # calculate quantile on uniform scale
    q_hat <- qcvine(u, alpha, vine = object$vine)
    
    # transform to original scale
    if (!uscale) {
      q_hat <- lapply(q_hat, qkde1d, obj = object$margins[[1]])
      if (inherits(object$model_frame[[1]], "ordered")) {
        # when response is discrete, we need to adjust the quantiles
        lvls <- levels(object$model_frame[[1]])
        q_hat <- lapply(q_hat, with_levels, lvls = lvls)
      }
    }
    
    ## always return as data frame
    q_hat <- as.data.frame(q_hat)
    names(q_hat) <- alpha
    if (!is.null(preds)) {
      preds <- cbind(preds, q_hat)
      alpha <- c("mean", alpha)
    } else {
      preds <- q_hat
    }
  } else {
    alpha <- "mean"
  }
  
  preds <- as.data.frame(preds)
  names(preds) <- alpha
  preds
}

#' @rdname predict.Cvinereg
#' @importFrom stats fitted
#' @export
fitted.Cvinereg <- function(object, alpha = 0.5, ...) {
  predict.Cvinereg(object, newdata = object$model_frame, alpha = alpha, ...)
}


predict_mean_Cvine <- function(object, newdata, uscale) {
  preds <- predict.Cvinereg(object, newdata, alpha = 1:10/ 11, uscale)
  data.frame(mean = rowMeans(preds))
}

with_levels <- function(q, lvls) {
  q <- ceiling(q)
  q <- pmax(q, 1)
  q <- pmin(q, length(lvls))
  q <- ordered(lvls[q], levels = lvls)
  q
}

# quantile estimation
qcvine <- function(u, alpha, vine) {
  d <- dim(vine)[1]
  if (ncol(u) != d - 1)
    stop("Dimensions of u and vine are not compatible")
  
  vine$structure <- as_rvine_structure(gen_cvine_mat(d))

  ## obtain diagonal entries in V matrix
  n <- nrow(u)
  V <- array(NA, dim = c(d, d, n))
  V[d, 2:d, ] <- t(u[, (d - 1):1])
  if (d > 2) {
    for (j in (d - 1):2) {
      for (k in 2:j) {
        tmp <- cbind(V[j + 1, j + 1, ], V[j + 1, k, ])
        V[j, k, ]  <- hbicop(tmp, 1, vine$pair_copulas[[d - j]][[k]])
      }
    }
    tmp <- t(apply(V, 3, diag)[-1, ]) 
  } else {
    tmp <- V[d, 2, ]
  }

  tmp <- pmin(pmax(tmp, 1e-10), 1 - 1e-10)
  
  # return as list (will be processed further)
  lapply(alpha,
         function(a) {
           inv.tmp <- inverted.tmp <- cbind(a, tmp)
           for (i in 1:(d-1)) {
             inverted.tmp[, i + 1] <- hbicop(cbind(inv.tmp[, i + 1], inverted.tmp[, i]), 1, 
                                     vine$pair_copulas[[d - i]][[1]],
                               inverse = TRUE)
           }
           return(inverted.tmp[, d])
         })
}




#qcvine <- function(u, alpha, vine) {
#  d <- dim(vine)[1]
#  if (ncol(u) != d - 1)
#    stop("Dimensions of u and vine are not compatible")
  
#  cpits <- rosenblatt(cbind(0.5, u[, (d - 1):1]), vine)[, -1, drop = FALSE]
 # print(cpits)
#  tmp <- pmin(pmax(cpits, 1e-10), 1 - 1e-10)
#  print(cor(tmp))
  
 # q_hat <- lapply(
 #   alpha,
#    function(a) inverse_rosenblatt(cbind(a, cpits), vine)[, 1]
#  )
  
#  q_hat <- as.data.frame(q_hat)
#  names(q_hat) <- alpha
#  q_hat
#}
