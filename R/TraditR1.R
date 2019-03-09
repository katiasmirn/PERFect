#' @export
##############################################
#Traditional Filtering Rule 1
##############################################
TraditR1 <- function(X, rel = FALSE, thresh=5){
  # Check the format of X
  if(!(class(X) %in% c("matrix"))){X <- as.matrix(X)}
  #   stop('X must be a data frame or a matrix')
  # if(!(class(X) == "matrix")){X <- as.matrix(X)}

  # Check the format of rel
  if(class(rel) != "logical")
    stop('rel argument must be a logical value')

  # Check the format of thresh
  if(!is.numeric(thresh)) stop('thresh argument must be a numerical value')

  #Select without setting filtering criteria
  TaxaAll <- colnames(X)
  NNzero <- apply(X, 2, nnzero)
  if(rel == TRUE){n <- dim(X)[1]
                  NNzero <- NNzero/n}
  filtX <- X[,NNzero >= thresh]
  return(filtX)
}

