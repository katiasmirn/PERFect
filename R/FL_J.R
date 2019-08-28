##################################
#Calculate FL for taxa left in the set
#################################
#' Filtering Loss for a set of filtered taxa J
#'
#' @description This function calculates filtering loss due to removing a group of J taxa.
#'
#' @usage FL_J(X, J)
#'
#' @param X OTU table, where taxa are columns and samples are rows of the table.
#' It should be a in data frame format with columns corresponding to taxa names.
#' @param J A vector of J taxa to be removed. It must be subset of column names of X.
#'
#' @return FL Filtering loss value
#'
#' @references Smirnova, E., Huzurbazar, H., Jafari, F. ``PERFect: permutation  filtration of microbiome data", to be submitted.
#'
#' @author Ekaterina Smirnova
#'
#' @seealso \code{\link{FiltLoss}}
#'
#' @examples
#' data(mock2)
#'
#' # Proportion data matrix
#' Prop <- mock2$Prop
#'
#' # Counts data matrix
#' Counts <- mock2$Counts
#'
#' #arrange counts in order of increasing number of samples taxa are present in
#' NP <- NP_Order(Counts)
#' Counts <- Counts[,NP]
#'
#' # Extract the taxa names to be removed
#' J <- colnames(Counts)[1:30]
#'
#' #Calculate filtering loss due to removing these taxa
#' FL_J(Counts,J)
#'
#' @export
#'

FL_J <- function(X, J){

  # Check the format of X
  if(!(class(X) %in% c("matrix"))){X <- as.matrix(X)}
  #   stop('X must be a data frame or a matrix')
  # if(!(class(X) == "matrix")){X <- as.matrix(X)}

  # Check the format of J
  if(class(J) != "character")
    stop('J argument must be a character vector containing names of taxa to be removed')

  Ind <- which(colnames(X) %in%  J)
  X_R <- X[,-Ind]
  #calculate corresponding norm
  Netw <- t(as.matrix(X))%*%as.matrix(X)
  Netw_R <- t(as.matrix(X_R))%*%as.matrix(X_R)
  #FL <-  1 - (psych::tr(t(Netw_R)%*%Netw_R)/psych::tr(t(Netw)%*%Netw))
  FL <-  1 - (sum(Netw_R*Netw_R)/sum(Netw*Netw))
  return(FL)
}
