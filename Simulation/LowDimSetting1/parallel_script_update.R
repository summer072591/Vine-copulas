##******************************************************************************************************
##               Candidate selection based on correlations
##
##******************************************************************************************************
candidate.selection.index <- function(correlations, candidate.number, root.index) {
  
  ## perform the candidate selection here 
  ## the inputs are the partial correlations between the candidates and the root calculated before
  ## the candidate number is the number of candidates to be selected
  ## the root.index is a vector recording the varibles which have been selected out already
  ## return the index of the candidate variables
  abs.correlation <- abs(correlations)
  
  if(missing(root.index)) root.index <- rep(0, length(correlations))
  ## ensure that the previous roots do not get selected
  abs.correlation[root.index] <- -Inf
  
  ## check how many variables are available for selection
  selectable.number <- length(correlations) - length(which(root.index!= 0))

  if(candidate.number <= selectable.number) 
    top.correlations.selection <- sort(abs.correlation, decreasing = TRUE)[1 : candidate.number]
  else top.correlations.selection <- sort(abs.correlation, decreasing = TRUE)[1 : selectable.number]
  
  return(which(abs.correlation %in% top.correlations.selection))
  
}

##******************************************************************************************************
##                 Partial selection calculation
## 
##******************************************************************************************************
## correlation selection 
## input is the udata matrix where the first column is related to the response
## udata should be either marginal udata or conditional udata
## ****************************************************************************************************
partial.cor.cal <- function(udata, selection.pointer, 
                            selection.var, method = c("pearson", "spearman", "kendall")) {
  
  if(selection.pointer == 1) return(correlation.cal(candidate = udata[, -1], 
                                                    root = udata[, 1], method = method))
  else return(apply(udata[, -1], 2, 
                    function(x) partial.correlation(pair = cbind(x, udata[, 1]), 
                                                    conditioning = udata[, selection.var[selection.pointer - 1] + 1], method = method)))
  
  
}

## ****************************************************************************************************
## ****************************************************************************************************

## used to calculate the paired partial correlation
## support u data
partial.correlation <- function(pair, conditioning, method) {
  ## input pair is a n \times 2 matrix, used to calculate the partial correlation 
  cor.prep <- c(cor(pair[, 1], pair[, 2], method = method), 
                correlation.cal(candidate = pair, root = conditioning, method = method))
  
  return((cor.prep[1] - cor.prep[2] * cor.prep[3]) / sqrt((1 - (cor.prep[2])^2) * (1 - (cor.prep[3])^2)))
  
}

## ****************************************************************************************************
##        Correlation calculation
##******************************************************************************************************
correlation.cal <- function(candidate, root, method = c("spearson", "kendal", "spearman")) {
  
  ## this function is used to calcualte the ordinary correlation between the candidate variables 
  ## with the root
  ## for selecting the top candidate variables with the highest correlation with the root
  ## in order to reduce the computation complexity
  
  ## candidate should be a matrix type, each column corresponds to a variable
  
  return(apply(candidate, 2, function(x) cor(x, root, method = method)))
  
}


##******************************************************************************************************
##        kernel density estimation -- marginal transformation to approximate PIT
##        alternatively one can perform discrete PIT
##******************************************************************************************************
kernel.marginal.transform <- function(data) {
  
  ## ***************************************************************************************************
  ## ***************************************************************************************************   
  ## data - a matrix with the first column being the response 
  ## requires R package ks 
  if("ks" %in% rownames(installed.packages()) == FALSE) install.packages("ks")
  library(ks)
  
  udata <- apply(data, 2, function(x) {
    kernel.est.tmp <- kcde(x)
    return(predict(kernel.est.tmp, x = x))
  }
  )
  
  return(udata)
}


##******************************************************************************************************
##******************************************************************************************************
##      the following several functions are implemented for qunatile prediction 
##      some are based on the implementations in the "vinereg" and "rvinecopulib" package
##******************************************************************************************************


##******************************************************************************************************
##                                   Main selection implementation
##******************************************************************************************************



Cvine.selection <- function(udata, candidate.number, selection.step, 
                           method = c("pearson", "spearman", "kendall"),
                           k.redu.method = c("all.partial.cor", "mix.partial.cor", "none"),
                           k.redu.percent = c(.1, .1)) {

  ## ***************************************************************************************************
  ## ***************************************************************************************************         
  ## udata - a matrix with the first column being the repsonse on the u scale
  ## method - the correlation calculation measure
  ## candidate number - the number of variables selected based on partial correlation
  ## selection step - the number of variables to be selected
  ## k.redu.method - method of reducing variables used to calculate the partial log-likelihood,
  ##                 "all.partial.cor" denotes selection based on only partial correlation 
  ##                 "mix.partial.cor" denotes mix of partial correlation selection
  ##                  and random selection
  ##                 "none" denotes no reduction
  ## k.redu.percent - percentage of the variables left, only used when k.redu.method != "none"
  ##                  the first entry denotes proportion selected by partial correlation
  ##                  the second entry denotes random selection; != 0 when k.redu.method == "mix.partial.cor"  
  ## requires R packages "kdecopula" and "doParallel"
  ## ***************************************************************************************************
  ## ***************************************************************************************************
  start.time <- Sys.time()
  p <- dim(udata)[2] - 1
  n <- dim(udata)[1]
    
  loglik.large.old <- loglik.large.new <- 
            loglik.compare <- list(var.index = -Inf, loglik = rep(-Inf, p), conditional.dist = NA)

  selection.var <- recycle.loglik <- numeric(p)
  
  
  if("kdecopula" %in% rownames(installed.packages()) == FALSE) install.packages("kdecopula")
  if("doParallel" %in% rownames(installed.packages()) == FALSE) install.packages("doParallel")
  
  
  library(kdecopula)
  library(doParallel)
  cores <- detectCores()
  
  if(cores[1] > 1) {
    cl <- makeCluster(cores[1] - 1) #not to overload your computer
   
  } else {
    cl <- 1
  }
  registerDoParallel(cl)
  ## ***************************************************************************************************
  ## ***************************************************************************************************
 
  for(selection.pointer in 1:selection.step) {
    
    print(paste("Selecting variable", selection.pointer))
  
     ## partial correlation or marginal correlation calculation 
  
     if(selection.pointer %in% c(1, 2)) partial.cor.selection <- 
         partial.cor.cal(udata = udata, selection.pointer = selection.pointer, 
                      selection.var = selection.var, method = method)
     else partial.cor.selection <- 
         partial.cor.cal(loglik.large.old$conditional.dist, selection.pointer = selection.pointer, 
                      selection.var = selection.var, method = method)
  
     ## the candidate variables in each selection step
     cand.index <- candidate.selection.index(correlations = partial.cor.selection, 
                              candidate.number = candidate.number, root.index = selection.var)
  

     
  ## ***************************************************************************************************
  ## ***************************************************************************************************
  
     if(selection.pointer == 1) loglik.large.new$conditional.dist <- udata
  
     loglik.large.old <- loglik.large.new
     loglik.large.new <- loglik.compare <- list()
  
     loglik.large.new.tmp <- 
       foreach(cand.pointer = cand.index, .packages = c("kdecopula", "doParallel"), 
               .export = c("partial.correlation", "correlation.cal"), .inorder = TRUE) %dopar% {
  #   for(cand.pointer in cand.index) {
      
           loglik.compare$conditional.dist <- matrix(rep(0, n * (p + 1)), c(n, (p + 1)))
        

           loglik.compare$var.index <- cand.pointer
           loglik.compare$loglik <- rep(-Inf, p)
      
           ## first component of the log-likelihood 
      
           if(selection.pointer == 1)  
             loglik.component1 <- logLik(kdecop(cbind(udata[, (cand.pointer + 1)], udata[, 1])))
           else loglik.component1 <- recycle.loglik[cand.pointer]
      
      
           ## second components of the log-likelihood
           cond.distribution <- matrix(rep(0, n*2), c(n, 2))
           
      
           ## selection
           k.index <-  which(!(1 : p %in% selection.var | 1 : p == cand.pointer))
           ## partial correlation selection
           ## reduce calculation of order O(p)
           if(k.redu.method == "none") k.index.reduced <- k.index
           else if(k.redu.method %in% c("all.partial.cor", "mix.partial.cor")) {
             
             if(selection.pointer == p) k.index.reduced <- k.index
             else {
              k.index.cor.record <- abs(sapply(k.index, 
                                 function(x) 
                                   partial.correlation(pair = cbind(loglik.large.old$conditional.dist[, (x + 1)], 
                                                                    loglik.large.old$conditional.dist[, 1]), 
                                                       conditioning = loglik.large.old$conditional.dist[, (cand.pointer + 1)], 
                                                       method = method)))
              redu.length <- floor(k.redu.percent * length(k.index))
              
             if(k.redu.method == "all.partial.cor") {
               
               ## only the first argument is used 
               if(redu.length[1] < 1) { 
                 k.index.reduced <- k.index
               } else {
                 k.index.reduced <- k.index[order(k.index.cor.record, decreasing = TRUE) %in% (1:redu.length[1])]
                 }
               
              
             }else if(k.redu.method == "mix.partial.cor") {
               
               ## first select based on partial correlation 
               if(redu.length[1] < 1) { 
                 k.index.reduced <- k.index
               } else {
               k.index.reduced <- k.index[order(k.index.cor.record, decreasing = TRUE) %in% (1:redu.length[1])]
               }
               ## the rest are chosen randomly 
               if(redu.length[2] >= 1) {
               k.index.reduced <- c(k.index.reduced, 
                                    sample(k.index[!(order(k.index.cor.record, decreasing = TRUE) %in% (1:redu.length[2]))],
                                           size = redu.length[2],
                                           replace = FALSE))
               } ## else don't do anything 
               
             }
             }
           }


      
           dist.fit.tmp.1 <- kdecop(cbind(loglik.large.old$conditional.dist[, (cand.pointer + 1)], 
                                          loglik.large.old$conditional.dist[, 1]))
      
       #    for(sample.pointer in 1:n) {
        
       #           hevaluate <- hkdecop(
       #             c(loglik.large.old$conditional.dist[sample.pointer, (cand.pointer + 1)],
       #               loglik.large.old$conditional.dist[sample.pointer, 1]), dist.fit.tmp.1, cond.var = 1)
        
        #          cond.distribution[sample.pointer, 1] <- hevaluate
           
        
         #  }
           
           
      hevaluate <- apply(cbind(loglik.large.old$conditional.dist[, (cand.pointer + 1)],
                   loglik.large.old$conditional.dist[, 1]), 
                   1, function(x) hkdecop(x, dist.fit.tmp.1, cond.var = 1))

      cond.distribution[, 1] <- hevaluate
      
      loglik.compare$conditional.dist[, 1] <- cond.distribution[, 1] 

      
      
      # for(k.pointer in k.index) {
      loglik.compare.tmp <- 
        foreach(k.pointer = k.index.reduced, .packages = "kdecopula", 
                .combine = function(...) mapply('cbind', ..., SIMPLIFY=FALSE),
                .multicombine=TRUE, .inorder = TRUE) %dopar% {
                  
                  
                  dist.fit.tmp.2 <- 
                    kdecop(cbind(loglik.large.old$conditional.dist[, (cand.pointer + 1)], 
                                 loglik.large.old$conditional.dist[, (k.pointer + 1)]))
                                      
                 # for(sample.pointer in 1:n) {
                 #   hevaluate <- 
                 #     hkdecop(c(loglik.large.old$conditional.dist[sample.pointer, (cand.pointer + 1)],
                 #               loglik.large.old$conditional.dist[sample.pointer, (k.pointer + 1)]),
                 #             dist.fit.tmp.2, cond.var = 1)
                    
                 #   cond.distribution[sample.pointer, 2] <- hevaluate
                    
              
                 #   }  
                                      
                  hevaluate <- apply(cbind(loglik.large.old$conditional.dist[, (cand.pointer + 1)],
                                       loglik.large.old$conditional.dist[, (k.pointer + 1)]), 1,
                                     function(x) hkdecop(x, dist.fit.tmp.2, cond.var = 1))
                  
                  cond.distribution[, 2] <- hevaluate                     
                  
                  loglik.component2 <- logLik(kdecop(cond.distribution))
                                      
                                      
                  loglik <- loglik.component1 + loglik.component2
                                      
                  conditional.dist <- cond.distribution[, 2]

                  list(loglik = loglik, conditional.dist = conditional.dist)
                  
                  
                }## the end of k.index loop
      
      loglik.compare$loglik[k.index.reduced] <- loglik.compare.tmp$loglik
      loglik.compare$conditional.dist[, (1 + k.index.reduced)] <- loglik.compare.tmp$conditional.dist

 
      loglik.compare
      
      } ## the end of candidate loop
  
  index <- which.max(sapply(loglik.large.new.tmp, function(x) max(x$loglik)))
  loglik.large.new <- loglik.large.new.tmp[[index]]
  
  
  

  selection.var[selection.pointer] <- loglik.large.new$var.index
  
  

  ## due to the reduced selection, the conditional distribution matrix is not complete
  ## need to fill in the columns corresponding to the complementary of k.index.reduced
  k.index <-  which(!(1 : p %in% selection.var | loglik.large.new$loglik != -Inf))

  cand.pointer <- selection.var[selection.pointer]
  if(selection.pointer == 1)  
    loglik.component1 <- logLik(kdecop(cbind(udata[, (cand.pointer + 1)], udata[, 1])))
  else loglik.component1 <- recycle.loglik[cand.pointer]
  
  ## second components of the log-likelihood
  cond.distribution <- cbind(loglik.large.new$conditional.dist[, 1], rep(0, n))

  
  loglik.compare.tmp <- 
    foreach(k.pointer = k.index, .packages = "kdecopula", 
              .combine = function(...) mapply('cbind', ..., SIMPLIFY=FALSE),
              .multicombine=TRUE, .inorder = TRUE) %dopar% {
                
              
                
                dist.fit.tmp.2 <- 
                  kdecop(cbind(loglik.large.old$conditional.dist[, (cand.pointer + 1)], 
                               loglik.large.old$conditional.dist[, (k.pointer + 1)]))

                #for(sample.pointer in 1:n) {
                #  hevaluate <- 
                #    hkdecop(c(loglik.large.old$conditional.dist[sample.pointer, (cand.pointer + 1)],
                #              loglik.large.old$conditional.dist[sample.pointer, (k.pointer + 1)]),
                #            dist.fit.tmp.2, cond.var = 1)
                  
               #   cond.distribution[sample.pointer, 2] <- hevaluate
                  
                #}  
                hevaluate <- apply(cbind(loglik.large.old$conditional.dist[, (cand.pointer + 1)],
                                     loglik.large.old$conditional.dist[, (k.pointer + 1)]), 1,
                                   function(x) hkdecop(x,  dist.fit.tmp.2, cond.var = 1))
                  
                cond.distribution[, 2] <- hevaluate
                

                loglik.component2 <- logLik(kdecop(cond.distribution))
                
                
                loglik <- loglik.component1 + loglik.component2
                
                conditional.dist <- cond.distribution[, 2]

                list(loglik = loglik, conditional.dist = conditional.dist)
                
                
              }## the end of k.index loop for completing the conditional dist matrix 
  
  ## fill in the conditional dist matrix
  loglik.large.new$loglik[k.index] <- loglik.compare.tmp$loglik
  loglik.large.new$conditional.dist[, (1 + k.index)] <- loglik.compare.tmp$conditional.dist

  ## recycle the loglikelihood 
  ## will be used for calculating the loglikelihood for the succeeding trees 
  recycle.loglik <- recycle.loglik + loglik.large.new$loglik
  

  } ## the end of variable search step loop 
  
  stopCluster(cl)
  
  end.time <- Sys.time()
  return(list(selection.var = selection.var, cond.dist = loglik.large.new$conditional.dist,
              run.time = end.time - start.time))
  
  
}



##******************************************************************************************************
##******************************************************************************************************
##  Next is an simulation example
##
##******************************************************************************************************
##******************************************************************************************************

#p <- 6
#n <- 100
#beta <- c(3, 1.5, 0, 0, 2, rep(0, p - 5))
#beta <- c(1:7, rep(0, p - 7))
# setting 1
#beta <- c(5, 5, 5, rep(0, p - 3))

## setting 2
#set.seed(1234)
#beta <- c(5, 5, 5, rep(0, 72), sample(c(2^(-2), 2^(-1), 1, 2), size = 75, replace = TRUE))



#p <- 150
#sigma.X <- function(p) {
  
  ## argument is p, covariance matrix is p \times p
#  sigma.tmp <- matrix(rep(0, p * p), c(p, p))
#  values <- 0.5^(0 : (p - 1))
  
  
#  for(i in 1 : p)
#    for (j in 1 : i) sigma.tmp[i, j] <- values[abs(i - j) + 1]
  
  ## flip the lower triangle to the upper triangle
#  diag.sigma.tmp <- diag(sigma.tmp)
#  sigma.tmp <- t(sigma.tmp) + sigma.tmp
#  diag(sigma.tmp) <- diag.sigma.tmp
  
#  return(sigma.tmp)
  
#}

#set.seed(123)
#X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = sigma.X(p))
#noise <- rnorm(n)

#y <- X %*% beta + noise 

#u.X <- apply(X, 2, function(x) rank(x) / (length(x) + 1))
#u.y <- rank(y) / (length(y) + 1)
#udata <- cbind(u.y, u.X)
#udata <- kernel.marginal.transform(cbind(y, X))

#selected.vars <- 
#  Cvine.selection(udata = udata, candidate.number = 5, selection.step = 5, method = "spearman",
#                  k.redu.method = "all.partial.cor", k.redu.percent = c(.2, .1))
#selected.vars





