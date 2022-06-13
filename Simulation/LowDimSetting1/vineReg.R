
#' 
#' @noRd
#' @importFrom stats pchisq
#' @importFrom rvinecopulib as_rvine_structure
finalize_Cvinereg_object <- function(formula, model_frame, margins, vine,
                                    status, var_nms) {
  ## adjust model matrix and names
  reorder <- status$selected_vars
  reorder[order(reorder)] <- seq_along(status$selected_vars)
  vine$structure <- as_rvine_structure(
    gen_cvine_mat(elements = c(1, rev(reorder + 1)))
  )
  covariate_names <- names(sort(status$selected_vars))
  vine$names <- c(var_nms[1], covariate_names)
  
  ## compute fit statistics
  nobs <- nrow(model_frame)
  var_edf <- status$edf
  var_cll <- status$cll
  var_caic <- -2 * var_cll + 2 * var_edf
  var_cbic <- -2 * var_cll + log(nobs) * var_edf
  var_p_value <- pchisq(2 * var_cll, var_edf, lower.tail = FALSE)
  var_p_value[1] <- NA
  cll <- sum(var_cll)
  edf <- sum(var_edf)
  caic <- sum(var_caic)
  cbic <- sum(var_cbic)
  
  stats <- list(
    nobs = nobs,
    edf = edf,
    cll = cll,
    caic = caic,
    cbic = cbic,
    var_edf = var_edf,
    var_cll = var_cll,
    var_caic = var_caic,
    var_cbic = var_cbic,
    var_p_value = var_p_value
  )
  
  ## return results as S3 object
  out <- list(
    formula = formula,
    selcrit = status$selcrit,
    model_frame = model_frame,
    margins = margins,
    vine = vine,
    stats = stats,
    order = var_nms[status$selected_vars + 1],
    selected_vars = status$selected_vars
  )
  class(out) <- "Cvinereg"
  out
}


check_order <- function(order, var_nms) {
  if (!all(order %in% var_nms))
    stop("unknown variable name in 'order'; ",
         "allowed values: '", paste(var_nms[-1], collapse = "', '"), "'.")
  if (any(order == var_nms[1]))
    stop("response variable '", var_nms[1],
         "' must not appear in 'order'.")
}


fit_margins <- function(x, par_1d, cores, uscale) {
  d <- ncol(x)
  par_1d <- process_par_1d(par_1d, d)
  if (uscale) {
    # data are uniform, no need to estimate margins
    margs <- lapply(
      seq_len(d),
      function(i) list(u = x[, i], loglik = 0, edf = 0)
    )
  } else {
    fit_margin <- function(k) {
      arg_lst <- list(
        x = x[, k],
        xmin = par_1d$xmin[k],
        xmax = par_1d$xmax[k],
        bw   = par_1d$bw[k],
        mult = par_1d$mult[k],
        deg  = par_1d$deg[k]
      )
      arg_lst[sapply(arg_lst, is.null)] <- NULL
      m <- do.call(kde1d, arg_lst)
      m$x_cc <- x[, k]
      m
    }
    margs <- furrr::future_map(seq_len(d), fit_margin)
  }
  
  names(margs) <- colnames(x)
  margs
}

get_pits <- function(margin_models, cores) {
  if (!is.null(margin_models[[1]]$u)) {
    # data are uniform, no need for PIT
    u <- sapply(margin_models, function(m) m$u)
  } else {
    get_pit <- function(m) pkde1d(m$x_cc, m)
    u <- furrr::future_map(margin_models, ~ pkde1d(.$x_cc, .))
    u <- do.call(cbind, u)
  }
  
  colnames(u) <- names(margin_models)
  u
}


process_par_1d <- function(pars, d) {
  if (!is.null(pars$xmin)) {
    if (length(pars$xmin) != d)
      stop("'xmin'  must be a vector with one value for each variable")
  }
  if (!is.null(pars$xmax)) {
    if (length(pars$xmax) != d)
      stop("'xmin'  must be a vector with one value for each variable")
  }
  if (!is.null(pars$xmax)) {
    if (length(pars$xmax) != d)
      stop("'xmin' must be a vector with one value for each variable")
  }
  if (length(pars$bw) != d && !is.null(pars$bw))
    stop("'bw' must be a vector with one value for each variable")
  
  if (is.null(pars$mult))
    pars$mult <- 1
  if (length(pars$mult) == 1)
    pars$mult <- rep(pars$mult, d)
  if (length(pars$mult) != d)
    stop("mult has to be of length 1 or the number of variables")
  
  if (is.null(pars$deg))
    pars$deg <- 2
  if (length(pars$deg) == 1)
    pars$deg <- rep(pars$deg, d)
  if (length(pars$deg) != d)
    stop("deg has to be of length 1 or the number of variables")
  
  pars
}

# u_k <- pkde1d(x[, k], fit)
# list(fit = fit, u = pkde1d(x[, k], fit))
#
initialize_fit <- function(u) {
  list(
    # 1-dimensional (= empty) vine
    vine = list(pair_copulas = list(list()), structure = as.matrix(1)),
    # array for storing pseudo-observations
    psobs = list(
      direct = array(u[, 1], dim = c(1, 1, nrow(u))),
      indirect = array(NA, dim = c(1, 1, nrow(u)))
    ),
    # conditional log-likelihood of the model
    cll = 0
  )
}

initialize_status <- function(margin_fits, selcrit) {
  list(
    # remaining variable indices to select from
    remaining_vars = seq_len(length(margin_fits) - 1),
    # variables indices included in the model
    selected_vars = NULL,
    # selection criterion
    selcrit = selcrit,
    # conditional logliklihood (only unconditional margin so for)
    clls = margin_fits[[1]]$loglik,
    # number of parameters in current model
    edf = margin_fits[[1]]$edf,
    # TRUE when no improvement is possible
    optimum_found = FALSE
  )
}

update_status <- function(status, new_vines) {
  clls <- sapply(new_vines, function(vine) vine$cll)
  edf <- sapply(new_vines, function(vine) vine$edf)
  n <- nrow(new_vines[[1]]$psobs[[1]])
  crits <- calculate_crits(clls, edf, n, status$selcrit)
  
  if (max(crits) < 0) {
    # optimum found, keep old fit
    status$optimum_found <- TRUE
  } else {
    status$best_ind <- which.max(crits)
    status$selected_vars <- c(
      status$selected_vars,
      status$remaining_vars[status$best_ind]
    )
    status$remaining_vars <- setdiff(
      status$remaining_vars,
      status$selected_vars
    )
    status$clls = c(status$clls, clls[status$best_ind])
    status$edf = c(status$edf, edf[status$best_ind])
  }
  
  status
}

calculate_crits <- function(clls, edf, n, selcrit) {
  clls - switch(
    selcrit,
    "loglik" = 0,
    "aic" = edf,
    "bic" = edf * log(n) / 2,
    "keep_all" = -Inf
  )
}

