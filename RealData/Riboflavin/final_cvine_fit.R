##**************************************************************************************
##**************************************************************************************

## this part is to fit the final cvine tree of based on the order of selected variables
## the function is based on xtcd_vine and gen_dvine_mat in vinereg package
## Aug 2019 version by J.Zhou

##**************************************************************************************

final_fit_cvine <- function(udata, order, old_fit, family_set, selcrit, ...) {
  
  d <- dim(udata)[2] 
  n <- dim(udata)[1]
  
  psobs <- list(
    direct = array(NA, dim = c(d, d, n)),
    indirect = array(NA, dim = c(d, d, n)))
  
  psobs$direct[d, 1, ] <- udata[, 1]
  psobs$direct[d, 2:d, ] <- t(udata[, rev(order)])
  
  old_fit$vine$pair_copulas <- list()
  
  edf <- 0
  for (i in ((d - 1) : 1)) {
    
    old_fit$vine$pair_copulas[[d - i]] <- list()
    
    for (var in seq_len(i)) {
      # get data for current edge
       u_e <- matrix(NA, n, 2)
       u_e[, 2] <- psobs$direct[(i + 1), var, ]     ## leaf
       u_e[, 1] <- psobs$direct[(i + 1), (i + 1), ]   ## root node
       
      # correct bandwidth for regression context
      # (optimal rate is n^(-1/5) instead of n^(-1/6))
       mult <- ifelse(is.null(list(...)$mult), 1, list(...)$mult)
       dots <- modifyList(list(mult = n^(1/6 - 1/5) * mult), list(...))
       args <- modifyList(
               list(
               data = u_e,
               family_set = family_set,
               selcrit = selcrit
               ),
               dots
               )
       args$threshold <- NULL
    
      # bivariate fit 
      if (!is.null(dots$threshold)) {
           if (abs(cor(args$data, method = "kendall")[1, 2]) < dots$threshold) {
           pc_fit <- bicop_dist()
           class(pc_fit) <- c("bicop", "bicop_dist")
           pc_fit$data <- args$data
           } else {
           pc_fit <- do.call(bicop, args)
           }
      } else {
           pc_fit <- do.call(bicop, args)
      }
       
      # store edges
      old_fit$vine$pair_copulas[[d - i]][[var]] <- pc_fit
  
      edf <- edf + pc_fit$npars
    
      # pseudo observations for next tree
      psobs$direct[i, var, ] <- hbicop(u_e, 1, pc_fit)
    } ## end of loop of tree i
    
  }## end of loop of all trees
  
  list(
    vine = vinecop_dist(old_fit$vine$pair_copulas, gen_cvine_mat(d)),
    psobs = psobs,
    cll = logLik(pc_fit),
    edf = edf
  )
}
