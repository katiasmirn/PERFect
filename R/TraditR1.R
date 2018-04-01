##############################################
#Traditional Filtering Rule 1
##############################################
TraditR1 <- function(X, rel = FALSE, thresh=5){
  #Select without setting filtering criteria
  TaxaAll <- colnames(X)
  NNzero <- apply(X, 2, nnzero)
  if(rel == TRUE){n <- dim(X)[1]
                  NNzero <- NNzero/n}
  filtX <- X[,NNzero >= thresh]
  return(filtX)
}