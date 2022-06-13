library(hdi)
data(riboflavin)

## 71 rows and 4088 columns 
ribo <- riboflavin

## y--log-tranformed riboflavin production rate
## x--measuring the logarithm of the expression level of 4088 genes
## n = 71, p = 4088


ribo.y <- ribo$y
ribo.x <- as.matrix(ribo$x)

ribo.sd <- apply(ribo.x, 2, sd)
ribo.rank <- order(ribo.sd, decreasing = TRUE)

## the top 500 gene expressions with the highest variance
## can adjust by users
n.selection <- 100
dim(ribo.X <- ribo.x[, ribo.rank %in% 1 : 100])

ribo.redu <- as.matrix(cbind(ribo.y, ribo.X))

## write selected data into .txt file
write.table(ribo.redu, file = "reducedRiboData.txt", col.names = TRUE, row.names = FALSE)
