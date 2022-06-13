##**************************************************************************************
##**************************************************************************************

## this part is to fit the final cvine tree of based on the order of selected variables
## the function is based on xtcd_vine and gen_dvine_mat in vinereg
## we use the xtnd_vine with the selected order with a different vine matrix 
## combine with a new cvine matrix 

##**************************************************************************************

xtnc_vine <- function(new_var, old_fit, family_set, selcrit, ...) {
  d <- dim(old_fit$psobs$direct)[1] + 1
  n <- length(new_var)
  
  psobs <- list(
    direct = array(NA, dim = c(d, d, n)),
    indirect = array(NA, dim = c(d, d, n))
  )
  psobs$direct[-1, -d, ] <- old_fit$psobs$direct
  psobs$indirect[-1, -d, ] <- old_fit$psobs$indirect
  psobs$direct[d, d, ] <- new_var
  
  print(psobs)
  
  old_fit$vine$pair_copulas[[d - 1]] <- list()
  edf <- 0
  for (i in rev(seq_len(d - 1))) {
    # get data for current edge
    u_e <- matrix(NA, n, 2)
    u_e[, 1] <- psobs$direct[i + 1, d, ]
    u_e[, 2] <- if (i == d - 1) {
      psobs$direct[i + 1, i + 1, ]
    } else {
      psobs$indirect[i + 1, i + 1, ]
    }
    
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
    
    
    # fit 
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
    
    old_fit$vine$pair_copulas[[d - i]][[i]] <- pc_fit
    
    edf <- edf + pc_fit$npars
    
    
    # pseudo observations for next tree
    psobs$direct[i, i, ] <- hbicop(u_e, 2, pc_fit)
    psobs$indirect[i, i, ] <- hbicop(u_e, 1, pc_fit)
  }
  
  list(
    vine = vinecop_dist(old_fit$vine$pair_copulas, gen_cvine_mat(d)),
    psobs = psobs,
    cll = logLik(pc_fit),
    edf = edf
  )
}

##******************************************************************************
## fit the final C-vine
## the last variable should be the response 
## this part should be organized later as 
##******************************************************************************

final_tree_fit <- function(udata, family_set, order...) {
  
  ## load initialize_fit function
  ## order specifies the order of selected predictors
  ##       should contain only none zero values and the last one should be repsonse
  
  current_fit <- initialize_fit(udata)
  for(i in order) current_fit <- xtnc_vine(udata[, i], current_fit, family_set)
  
  current_fit
}

#order <- selected.vars$selection.var
#order <- c((order[order != 0] + 1), 1)
#tree_fit <- final_tree_fit(udata, family_set = "nonparametric", order = order)

## the one with the largest index is the response variable
#plot(tree_fit$vine, tree = 1 : length(order))
###

## the prediction function etc would be the same as those implemented in vinereg
## difference would be the vinereg function 
## searching procedure is different, final fit to cvine is different 