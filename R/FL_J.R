##################################
#Calculate FL for taxa left in the set
#################################
FL_J <- function(X, J){
  Ind <- which(colnames(X) %in%  J)
  X_R <- X[,-Ind]
  #calculate corresponding norm
  Netw <- t(as.matrix(X))%*%as.matrix(X)
  Netw_R <- t(as.matrix(X_R))%*%as.matrix(X_R)
  FL <-  1 - (tr(t(Netw_R)%*%Netw_R)/tr(t(Netw)%*%Netw))
  
  return(FL)
}
