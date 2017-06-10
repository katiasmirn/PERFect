##############################################
#Traditional Filtering Rule 1
##############################################
TraditR1 <- function(X, rel = FALSE, thresh=5){
  #Select without setting filtering criteria
  TaxaAll <- names(Counts)
  NNzero <- apply(Counts, 2, nnzero)
  if(rel == TRUE){n <- dim(X)[1]
                  NNzero <- NNzero/n}
  filtX <- X[,NNzero >= thresh]
  return(filtX)
}