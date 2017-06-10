##################################
#Calculate FL for taxa left in the set
#################################
FL_J <- function(X, J, leave = TRUE){
  Ind <- which(names(X) %in%  J)
  if(leave == TRUE){X_R <- X[,Ind]}
  else{X_R <- X[,-Ind]}
  #calculate corresponding norm
  Netw <- t(as.matrix(X))%*%as.matrix(X)
  Netw_R <- t(as.matrix(X_R))%*%as.matrix(X_R)
  FL <- 1 - tr(t(Netw_R)%*%Netw_R)/tr(t(Netw)%*%Netw)
  
  return(FL)
}
