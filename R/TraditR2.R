#' @export
##############################################
#Traditional Filtering Rule 2
##############################################
TraditR2 <- function(X,  Ab_min = 0.001){

  # Check the format of X
  if(!(class(X) %in% c("matrix"))){X <- as.matrix(X)}
  #   stop('X must be a data frame or a matrix')
  # if(!(class(X) == "matrix")){X <- as.matrix(X)}

  # Check the format of Ab_min
  if(!is.numeric(Ab_min)) stop('Ab_min argument must be a numerical value')

  #check if X is a relative abundance matrix
  if(!(all(apply(X, 1, sum) ==1))) {X <- sweep(X, MARGIN=1, apply(X,1,sum), '/')}
  n <- dim(X)[1]
  #Min abundance level is always 0, select max expression of a taxon
  Abund <- apply(X, 2, max)
  Nsamp <- apply(X, 2, nnzero)
  Psamp <- Nsamp/n
  #Select only taxa with abundance level > 0.001
  if(is.null(names(Abund))){names(Abund) <- 1:length(Abund)}
  Taxa <- names(Abund[Abund > Ab_min])
  Abund <- Abund[Taxa]
  Nsamp <- Nsamp[Taxa]
  Psamp <- Psamp[Taxa]
  #Check additional conditions
  C1 <- names(Abund[Abund > 0.01 & Nsamp >=1])
  C2 <- names(Abund[Abund > 0.001 & Psamp >=0.02])
  C3 <- names(Abund[Psamp >=0.05])
  Taxa <- unique(c(C1, C2, C3))
  filtX <- X[,Taxa]
  return(filtX)
}
