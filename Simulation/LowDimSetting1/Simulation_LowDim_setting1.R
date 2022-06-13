#' **************************************************************************************************
#' this part corresponds to M5 in Kraus and Czado (2017)
#' 
#' **************************************************************************************************
#' 
#' p <- 70
#beta <- c(3, 1.5, 0, 0, 2, rep(0, p - 5))

library(np)
p <- 4


interval_score <- function(fit, lb, ub, alpha = 0.05) {
  # https://www.stat.washington.edu/raftery/Research/PDF/Gneiting2007jasa.pdf
  # Equation 43
  S <- (ub - lb) + 2 / alpha * (lb - fit) * as.numeric(fit < lb) + 
    2 / alpha * (fit - ub) * as.numeric(fit > ub)
  return(mean(S))
}

## quantile loss function 
quantile_loss <- function(x, tau) {
  
  return(x * (tau - I(x<= 0)))
}


sigma.X <- function(p) {
  
  ## argument is p, covariance matrix is p \times p
  sigma.tmp <- matrix(rep(0, p * p), c(p, p))
  values <- 0.5^(0 : (p - 1))
  
  
  for(i in 1 : p)
    for (j in 1 : i) sigma.tmp[i, j] <- values[abs(i - j) + 1]
  
  ## flip the lower triangle to the upper triangle
  diag.sigma.tmp <- diag(sigma.tmp)
  sigma.tmp <- t(sigma.tmp) + sigma.tmp
  diag(sigma.tmp) <- diag.sigma.tmp
  
  return(sigma.tmp)
  
}


## training and test sample in total
n <- 450



lowDim_single <- function(s) {
  
   set.seed(s)
  
   ## adjust the standard deviation accordingly
   sd.noise <- .1
   noise <- rnorm(n, sd = sd.noise)
   X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = sigma.X(p))
  
   y <- sqrt(abs(2 * X[, 1] - X[, 2] + 0.5)) + (-0.5 * X[, 3] + 1) * (0.1 * X[, 4]^3) + noise

   y.quant05.theory <- sqrt(abs(2 * X[, 1] - X[, 2] + 0.5)) + (-0.5 * X[, 3] + 1) * (0.1 * X[, 4]^3) + sd.noise * qnorm(0.5)
   y.quant095.theory <- sqrt(abs(2 * X[, 1] - X[, 2] + 0.5)) + (-0.5 * X[, 3] + 1) * (0.1 * X[, 4]^3) + sd.noise * qnorm(0.95)



   ## sample data
   ## adjust accordingly
   n.train <- 300
   data.train <- data.frame(cbind(y, X)[1:n.train, ])
  
   n.test <- 150
   data.test <- data.frame(cbind(y, X)[(n.train + 1) : n, ])

   #y.quant05.theory <- y.quant05.theory[(n.train + 1) : n]
   #y.quant095.theory <- y.quant095.theory[(n.train + 1) : n]

 
 
   colnames(data.train) <- colnames(data.test) <- c("y", "x1", "x2", "x3", "x4")
  
  
   source("quantile prediction.R")
   udata.train <- kernel.marginal.transform(data.train)
  
   selected.vars <- 
     Cvine.selection(udata = udata.train, candidate.number = 4, selection.step = 4, method = "spearman",
                    k.redu.method = "none")
  
   order <- colnames(data.train)[selected.vars$selection.var + 1]
   fit_cvine <- Cvinereg(y ~ ., data = data.train, order = order, family_set = "tll")
   yhat.cvine <- predict(fit_cvine,  newdata = data.test[, -1], alpha = NA)
   yhat.cvine.quant <- predict(fit_cvine,  newdata = data.test[, -1], alpha = c(.025, .05, 0.5, 0.95, .975))


  
   fit_dvine <- vinereg(y ~ ., data.train, family_set = "nonparametric")
   yhat.dvine <- predict(fit_dvine,  newdata = data.test[, -1], alpha = NA)
   yhat.dvine.quant <- predict(fit_dvine,  newdata = data.test[, -1], alpha = c(.025, .05, 0.5, 0.95, .975))

  
   
   ols <- lm(y ~ .- 1, data = data.train)
   beta <- ols$coefficients
   
   yhat.ols <- as.matrix(data.test[, -1]) %*% beta
   
   bw <- npregbw(xdat = data.train[, -1], ydat = data.train[, 1])
   yhat.np <- fitted(npreg(bws = bw, newdata = data.test[, -1]))
   
   write.table(t(c(s, mean((yhat.np - data.test[, 1])^2))), file = "mean_np_lowdim_setting1.txt",
               append = TRUE, col.names = FALSE, row.names = FALSE)

   
   write.table(t(c(s, mean((yhat.ols - data.test[, 1])^2))), file = "mean_ols_lowdim_setting1.txt",
               append = TRUE, col.names = FALSE, row.names = FALSE)
   
   
   write.table(t(c(s, order)), file = "selected_vars_lowdim_setting1.txt",
               append = TRUE, col.names = FALSE, row.names = FALSE)
   
   
   write.table(t(c(s, mean((yhat.cvine$mean - data.test[, 1])^2),
                   mean((yhat.dvine$mean - data.test[, 1])^2))), file = "mean_compare_lowdim_setting1.txt",
               append = TRUE, col.names = FALSE, row.names = FALSE)
   
   
   ## cvine quantile prediction record
   write.table(t(c(s, yhat.cvine.quant$`0.025`)), file = "quant0025_cvine_lowdim_setting1.txt",
               append = TRUE, col.names = FALSE, row.names = FALSE)
   
   write.table(t(c(s, yhat.cvine.quant$`0.05`)), file = "quant005_cvine_lowdim_setting1.txt",
               append = TRUE, col.names = FALSE, row.names = FALSE)
   
   
   write.table(t(c(s, yhat.cvine.quant$`0.5`)), file = "quant05_cvine_lowdim_setting1.txt",
               append = TRUE, col.names = FALSE, row.names = FALSE)
   
   write.table(t(c(s, yhat.cvine.quant$`0.95`)), file = "quant095_cvine_lowdim_setting1.txt",
               append = TRUE, col.names = FALSE, row.names = FALSE)
   
   write.table(t(c(s, yhat.cvine.quant$`0.975`)), file = "quant0975_cvine_lowdim_setting1.txt",
               append = TRUE, col.names = FALSE, row.names = FALSE)
   
   
   
   ## dvine quantile prediction record
   write.table(t(c(s, yhat.dvine.quant$`0.025`)), file = "quant0025_dvine_lowdim_setting1.txt",
               append = TRUE, col.names = FALSE, row.names = FALSE)
   
   write.table(t(c(s, yhat.dvine.quant$`0.05`)), file = "quant005_dvine_lowdim_setting1.txt",
               append = TRUE, col.names = FALSE, row.names = FALSE)
   
   write.table(t(c(s, yhat.dvine.quant$`0.5`)), file = "quant05_dvine_lowdim_setting1.txt",
               append = TRUE, col.names = FALSE, row.names = FALSE)
   
   write.table(t(c(s, yhat.dvine.quant$`0.95`)), file = "quant095_dvine_lowdim_setting1.txt",
               append = TRUE, col.names = FALSE, row.names = FALSE)
   
   write.table(t(c(s, yhat.dvine.quant$`0.975`)), file = "quant0975_dvine_lowdim_setting1.txt",
               append = TRUE, col.names = FALSE, row.names = FALSE)
   
   
   
   write.table(t(c(s, yhat.cvine$mean)), file = "mean_cvine_lowdim_setting1.txt",
               append = TRUE, col.names = FALSE, row.names = FALSE)
   
   write.table(t(c(s, mean(sapply(yhat.cvine.quant$`0.05` - data.test[, 1], quantile_loss, tau = 0.05)),
                   mean(sapply(yhat.dvine.quant$`0.05` - data.test[, 1], quantile_loss, tau = 0.05)))),
               file = "quantile_loss_tau005_compare_lowdim_setting1.txt",
               append = TRUE, col.names = FALSE, row.names = FALSE)
   
   write.table(t(c(s, mean(sapply(yhat.cvine.quant$`0.5` - data.test[, 1], quantile_loss, tau = 0.5)),
                   mean(sapply(yhat.dvine.quant$`0.5` - data.test[, 1], quantile_loss, tau = 0.5)))),
               file = "quantile_loss_tau05_compare_lowdim_setting1.txt",
               append = TRUE, col.names = FALSE, row.names = FALSE)
   
   write.table(t(c(s, mean(sapply(yhat.cvine.quant$`0.95` - data.test[, 1], quantile_loss, tau = 0.95)),
                   mean(sapply(yhat.dvine.quant$`0.95` - data.test[, 1], quantile_loss, tau = 0.95)))),
               file = "quantile_loss_tau095_compare_lowdim_setting1.txt",
               append = TRUE, col.names = FALSE, row.names = FALSE)
   
  
   
   ## Interval prediction
   write.table(t(c(s, interval_score(fit = data.test[, 1], lb = yhat.cvine.quant$`0.975`,
                                     ub = yhat.cvine.quant$`0.025`, alpha = 0.05),
                   interval_score(fit = data.test[, 1], lb = yhat.dvine.quant$`0.975`,
                                  ub = yhat.dvine.quant$`0.025`, alpha = 0.05))),
               file = "interval_score_compare_lowdim_setting1.txt",
               append = TRUE, col.names = FALSE, row.names = FALSE)
   

}

