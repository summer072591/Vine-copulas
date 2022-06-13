
# Interval score
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



## leave-one-out cross validation
ribo  <- function(s) {
  set.seed(s)
  ribo.data <- read.table("reducedRiboData.txt", header = TRUE, sep = " ")
  p <- dim(ribo.data)[2] - 1
  n <- dim(ribo.data)[1]
  
  index <- sample.int(n = n, size = 10, replace = FALSE)
  
  data.train <- ribo.data[-index, ]
  #pairs(data.train[, 1:5])
  
  data.test <- ribo.data[index, ]
  
  
  source("quantile prediction.R")
  #dat <- data.frame(y = y, X = X.1)
  udata.train <- kernel.marginal.transform(data.train)
  
  selected.vars <- 
    Cvine.selection(udata = udata.train, candidate.number = 10, selection.step = p, method = "spearman",
                    k.redu.method = "mix.partial.cor", k.redu.percent = c(.2, .1))
  
  order <- colnames(data.train)[selected.vars$selection.var + 1]
  
  fit_cvine <- Cvinereg(ribo.y  ~ ., data = data.train, order = order, family_set = "nonparametric")
  yhat.cvine <- predict(fit_cvine,  newdata = data.test[, -1], alpha = NA)
  yhat.cvine.quant <- predict(fit_cvine,  newdata = data.test[, -1], alpha = c(0.025, 0.05, 0.5, 0.95, 0.975))
  
  
  fit_dvine <- vinereg(ribo.y ~ ., data.train, family_set = "nonparametric")
  
  # summary(fit)
  yhat.dvine <- predict(fit_dvine,  newdata = data.test[, -1], alpha = NA)
  yhat.dvine.quant <- predict(fit_dvine,  newdata = data.test[, -1], alpha = c(0.025, 0.05, 0.5, 0.95, 0.975))
  
  
  write.table(t(c(s, order)), file = "selected_vars_highdim.txt",
              append = TRUE, col.names = FALSE, row.names = FALSE)
  
  
  write.table(t(c(s, mean((yhat.cvine$mean - data.test[, 1])^2), 
                  mean((yhat.dvine$mean - data.test[, 1])^2))), file = "mean_compare_highdim.txt",
              append = TRUE, col.names = FALSE, row.names = FALSE)
  
  ## cvine quantile prediction record
  write.table(t(c(s, yhat.cvine.quant$`0.025`)), file = "quant0025_cvine_highdim.txt",
              append = TRUE, col.names = FALSE, row.names = FALSE)
  
  write.table(t(c(s, yhat.cvine.quant$`0.05`)), file = "quant005_cvine_highdim.txt",
              append = TRUE, col.names = FALSE, row.names = FALSE)
  
  
  write.table(t(c(s, yhat.cvine.quant$`0.5`)), file = "quant05_cvine_highdim.txt",
              append = TRUE, col.names = FALSE, row.names = FALSE)
  
  write.table(t(c(s, yhat.cvine.quant$`0.95`)), file = "quant095_cvine_highdim.txt",
              append = TRUE, col.names = FALSE, row.names = FALSE)
  
  write.table(t(c(s, yhat.cvine.quant$`0.975`)), file = "quant0975_cvine_highdim.txt",
              append = TRUE, col.names = FALSE, row.names = FALSE)
  
  
  
  ## dvine quantile prediction record
  write.table(t(c(s, yhat.dvine.quant$`0.025`)), file = "quant0025_dvine_highdim.txt",
              append = TRUE, col.names = FALSE, row.names = FALSE)
  
  write.table(t(c(s, yhat.dvine.quant$`0.05`)), file = "quant005_dvine_highdim.txt",
              append = TRUE, col.names = FALSE, row.names = FALSE)
  
  write.table(t(c(s, yhat.dvine.quant$`0.5`)), file = "quant05_dvine_highdim.txt",
              append = TRUE, col.names = FALSE, row.names = FALSE)
  
  write.table(t(c(s, yhat.dvine.quant$`0.95`)), file = "quant095_dvine_highdim.txt",
              append = TRUE, col.names = FALSE, row.names = FALSE)
  
  write.table(t(c(s, yhat.dvine.quant$`0.975`)), file = "quant0975_dvine_highdim.txt",
              append = TRUE, col.names = FALSE, row.names = FALSE)
  
  
  
  write.table(t(c(s, yhat.cvine$mean)), file = "mean_cvine_highdim.txt",
              append = TRUE, col.names = FALSE, row.names = FALSE)
  
  write.table(t(c(s, mean(sapply(yhat.cvine.quant$`0.05` - data.test[, 1], quantile_loss, tau = 0.05)), 
                  mean(sapply(yhat.dvine.quant$`0.05` - data.test[, 1], quantile_loss, tau = 0.05)))), 
              file = "quantile_loss_tau005_compare_highdim.txt",
              append = TRUE, col.names = FALSE, row.names = FALSE)
  
  write.table(t(c(s, mean(sapply(yhat.cvine.quant$`0.5` - data.test[, 1], quantile_loss, tau = 0.5)), 
                  mean(sapply(yhat.dvine.quant$`0.5` - data.test[, 1], quantile_loss, tau = 0.5)))), 
              file = "quantile_loss_tau05_compare_highdim.txt",
              append = TRUE, col.names = FALSE, row.names = FALSE)
  
  write.table(t(c(s, mean(sapply(yhat.cvine.quant$`0.95` - data.test[, 1], quantile_loss, tau = 0.95)), 
                  mean(sapply(yhat.dvine.quant$`0.95` - data.test[, 1], quantile_loss, tau = 0.95)))), 
              file = "quantile_loss_tau095_compare_highdim.txt",
              append = TRUE, col.names = FALSE, row.names = FALSE)
  
 
  
  ## Interval prediction
  write.table(t(c(s, interval_score(fit = data.test[, 1], lb = yhat.cvine.quant$`0.975`, 
                                    ub = yhat.cvine.quant$`0.025`, alpha = 0.05), 
                  interval_score(fit = data.test[, 1], lb = yhat.dvine.quant$`0.975`, 
                                 ub = yhat.dvine.quant$`0.025`, alpha = 0.05))), 
              file = "interval_score_compare_highdim.txt",
              append = TRUE, col.names = FALSE, row.names = FALSE)
  
  
  
  
  
}


