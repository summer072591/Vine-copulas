#'******************************************************************************************************
#' the following several functions are implemented for qunatile prediction 
#' some are based on the implementations in the "vinereg" and "rvinecopulib" package
#' main fitting function Cvinereg modified from vinereg function in vinereg package
#' used for fitting a c-vine with a speficied root node order
#' Aug 2019 version by J.Zhou
#'******************************************************************************************************
source("cpit.R")
source("eval_utils.R")
source("vinereg.R")
source("generic C_vine.R")
source("parallel_script_update_rvinecopulib.R")
#source("parallel_script_update.R")
source("predict.Cvinereg.R")
source("final_cvine_fit.R")
library(rvinecopulib)
library(vinereg)
library(kde1d)

#'******************************************************************************************************

#' Generates a C-Vine Matrix of dimension p. An order can be specified with the elements argument
#' this function is based on the gen_dvine_mat function implemented by D.Klaus
#' using elements argument can specify a root order

#'******************************************************************************************************


gen_cvine_mat <- function(p = length(elements), elements = seq_len(p)) {
  
  if (!setequal(seq_len(p), elements))
    stop("'elements' must consist of the numbers from 1 to p.")
  
  mat <- diag(rev(elements))
  for (i in seq_len(p - 1) + 1) mat[i - 1, (i : p)] <- rev(elements)[i - 1]
  
  ## reverse for vinecopulib notation
  mat[, p : 1]
}

#'******************************************************************************************************

#'******************************************************************************************************
Cvinereg <-
  function (formula, data, family_set = "nonparametric", selcrit = "loglik", 
          order = NA, par_1d = list(), cores = 1, uscale = FALSE,...) {
  
  if (!missing(data)) {
    mf <- model.frame(formula, data)
  }
  else {
    mf <- model.frame(formula, parent.frame())
  }
  if (!(is.ordered(mf[[1]]) | is.numeric(mf[[1]]))) 
    stop("response must be numeric or ordered")
  x <- cctools::cont_conv(mf)
  d <- ncol(x)
  cores <- min(cores, future::availableCores())
  if (cores > 1) {
    future::plan(future::multiprocess, workers = cores)
    on.exit(future::plan(), add = TRUE)
  }
  else {
    future::plan(future::sequential)
  }
  margin_models <- fit_margins(x, par_1d, cores, uscale)
  u <- get_pits(margin_models, cores)  ## original order
  

  ## root node order selection 
  #if(is.na(order)) {
    
  #  if(missing(candidate.number)) candidate.number <- dim(u)[2] - 1
  #  if(missing(selection.step)) selection.step <- dim(u)[2] - 1
    
   # selected.vars <- Cvine.selection(udata = u, candidate.number = candidate.number, 
  #                  selection.step = selection.step, method = method,
  #                  k.redu.method = k.redu.method)
  
  #  print(selected.vars$selection.var)
  #  order <- colnames(u)[selected.vars$selection.var + 1]
   # print(order)
    
  #}
  
  
  var_nms <- colnames(u)
  current_fit <- initialize_fit(u)
  status <- initialize_status(margin_models, selcrit)
  

  check_order(order, var_nms)
  status$selcrit <- "keep_all"
  current_fit <- final_fit_cvine(u, order, current_fit, family_set, ## this part changes to cvine fit
                               selcrit, ...)
  status <- update_status(status, list(current_fit))
      
      
  status$selected_vars <- sapply(order, function(nm) which(var_nms == nm) - 1, simplify = TRUE)
    
  status$remaining_vars <- which(!((1 : (d - 1)) %in% status$selected_vars))

  finalize_Cvinereg_object(formula = formula, model_frame = mf, 
                          margins = margin_models, vine = current_fit$vine, status = status, 
                          var_nms = var_nms)
}

#'******************************************************************************************************

# simulate data
#x <- matrix(rnorm(200), 100, 2)
#y <- x %*% c(1, -2)
#dat <- data.frame(y = y, x = x, z = rbinom(100, 1, prob = 0.5))
#udat <- apply(dat, 2, function(x) rank(x) / (length(x)+ 1))[, -1]

# fit vine regression model
#fit <- Cvinereg(y ~ ., dat, order = c("x.2", "x.1", "z"), family_set = "nonparametric")
#summary(fit)#fit


#fit2 <- vinereg(y ~ ., dat,order = c("x.2", "x.1", "z"),  family_set = "nonparametric")

#plot_effects(fit)
#plot_effects(fit2)

# model predictions
#mu_hat  <- predict(fit, newdata = dat, alpha = NA)          # mean
#med_hat <- predict(fit, newdata = dat, alpha = 0.5)         # median
#plot(cbind(y, mu_hat))

#mu_hat2  <- predict(fit2, newdata = dat, alpha = NA)          # mean
#med_hat2 <- predict(fit2, newdata = dat, alpha = 0.5)         # median
#plot(cbind(y, mu_hat2))

