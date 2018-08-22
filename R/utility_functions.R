########################################
#Function to calculate j^th DFL loss
#' @export
DiffFiltLoss_j <- function(Perm_Order,Netw, j){
  J_jp1 <-  Perm_Order[1:(j+1)]
  #calculate corresponding norm ratios this is the value of the matric ||
  DFL <-  (2*t(Netw[J_jp1[max(NROW(J_jp1))],-J_jp1])%*%Netw[J_jp1[max(NROW(J_jp1))],-J_jp1]+
             Netw[J_jp1[max(NROW(J_jp1))], J_jp1[max(NROW(J_jp1))]]^2)

  return(DFL)
}
##################################
#Randomly permute labels n times for a fixed j
#################################
#' @export
Perm_j_s <- function(j, Netw, k,p, p2 = NULL){
  #Netw = X'X - data matrix with taxa in columns and samples in rows
  #k - number of permutations used
  #create a list of k*p arraments of orders
  if(is.null(p2)){p2 <- p}
  labl <- sapply(1:k,function(x) NULL)
  labl <- lapply(labl,function(x)  sample(1:p,p2))
  FL_j <- sapply(labl,DiffFiltLoss_j, Netw = Netw, j=j)
  return(FL_j)
}

###################################################
#Obtain sampling distribution using permutations of DFL
###################################################
#' @export
sampl_distr <- function(X, k){
p <- dim(X)[2]
Netw <- t(X)%*%X
full_norm <- tr(t(Netw)%*%Netw)#this is full norm value
#For each taxon j, create a distribution of its DFL's by permuting the labels
res_all <- lapply(1:(p-1),function(x) x)
FL_j <- lapply(res_all, function(x) Perm_j_s(j = x, Netw =Netw, k=k, p =p, p2 = x+1))
#divide by the full matrix norm values
res_pres <- lapply(FL_j, function(x) {x/full_norm})
return(res_pres)
}

#######################################
#Compare results of filtering rules
####################################
#' @export
ResultsComparison <-function(X, Counts,Ab_min = 0.001,rel =FALSE, thresh=5,prop =FALSE,
                             res_sim, res_perm){

  #Traditional filtering
  Filt_R1 <- TraditR1(X=Counts,  rel =rel, thresh=thresh)#if rel = TRUE then use propo of samples > thresh (say 0.05 = 5%)
  Filt_R2 <- TraditR2(X=Counts,  Ab_min = Ab_min)#prop =FALSE if counts matrix is used
  #taxa left after traditional rules 1 and 2 applied
  Trad_R1 <- names(Filt_R1)
  Trad_R2 <- names(Filt_R2)

  #taxa left after traditional rules 1 and 2 applied
  Loss_Trad_R1 <- FL_J(X = X, J = Trad_R1)
  Loss_Trad_R2 <- FL_J(X = X, J = Trad_R2)

  #PERFect loss
  PERFect_sim_taxa <- names(res_sim$filtX)
  PERFect_perm_taxa <- names(res_perm$filtX)
  #corresponding filtering loss
  Loss_PERFect_sim <- FL_J(X = X, J = PERFect_sim_taxa)
  Loss_PERFect_perm <- FL_J(X = X, J = PERFect_perm_taxa)

  #table for the paper
  Res <- cbind(rbind(length(Trad_R1), length(Trad_R2), length(PERFect_sim_taxa), length(PERFect_perm_taxa)),
               rbind(Loss_Trad_R1, Loss_Trad_R2, Loss_PERFect_sim, Loss_PERFect_perm))
  colnames(Res) <- c("N_taxa", "FL")
  rownames(Res) <- c("Trad_R1", "Trad_R2", "PERFect_sim", "PERFect_perm")
  return(list(Res = Res,  Filt_R1 = Filt_R1, Filt_R2 = Filt_R2, res_sim = res_sim, res_perm = res_perm))
}
###########################################################################################
#Order functions
##########################################################################################
#' Taxa importance ordering by the number of occurrences of the taxa in the n samples
#'
#' @usage NP_Order(Counts)
#'
#' @param Counts OTU COUNTS table, where taxa are columns and samples are rows of the table.
#' It should be a in data frame format with columns corresponding to taxa names.
#'
#' @return NP Taxa names in increasing  order of the number of samples taxa are present in.
#'
#' @references Smirnova, E., Huzurbazar, H., Jafari, F. ``PERFect: permutation  filtration of microbiome data", to be submitted.
#'
#' @author Ekaterina Smirnova
#'
#' @examples
#' data(mock2)
#' # Proportion data matrix
#' Prop <- mock2$Prop
#'
#' # Counts data matrix
#' Counts <- mock2$Counts
#'
#' #arrange counts in order of increasing number of samples taxa are present in
#' NP <- NP_Order(Counts)
#' @export

NP_Order <- function(Counts){
  #arrange counts in order of increasing number of samples taxa are present in
  NP <- names(sort(apply(Counts, 2, nnzero)))
  return(NP)
}
#############################################################
#' Taxa importance ordering by PERFect p-values
#'
#' @description This function orders taxa by increasing significance of simultaneous PERFect p-values.
#'
#' @usage pvals_Order(Counts, res_sim)
#'
#' @param Counts OTU COUNTS table, where taxa are columns and samples are rows of the table.
#' It should be a in data frame format with columns corresponding to taxa names.
#'
#' @param res_sim Output of \code{PERFect_sim()} function.
#'
#' @return Order_pvals Taxa names in increasing order of p-values significance.
#'
#' @references Smirnova, E., Huzurbazar, H., Jafari, F. ``PERFect: permutation  filtration of microbiome data", to be submitted.
#'
#' @author Ekaterina Smirnova
#'
#' @seealso \code{\link{PERFect_sim}}, \code{\link{PERFect_perm}}
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
#' # Perform simultaenous filtering of the data
#' res_sim <- PERFect_sim(X=Counts)
#'
#' #order according to p-values
#' pvals_sim <- pvals_Order(Counts, res_sim)
#'
#' @export

pvals_Order <- function(Counts, res_sim){
  vec <- names(res_sim$pvals)
  Order <- c(setdiff(colnames(Counts), vec), vec)
  Order_pvals <- Order[c(1,sort.int(res_sim$pvals, index.return = TRUE, decreasing = TRUE)$ix +1)]
  return(Order_pvals)
}

#######################################################################
#' Taxa importance ordering by the number ofconnected taxa
#'
#' @usage NC_Order(Counts)
#'
#' @param Counts OTU COUNTS table, where taxa are columns and samples are rows of the table.
#' It should be a in data frame format with columns corresponding to taxa names.
#'
#' @return NC Taxa names in increasing  order of the number of connected taxa
#'
#' @references Smirnova, E., Huzurbazar, H., Jafari, F. ``PERFect: permutation  filtration of microbiome data", to be submitted.
#'
#' @author Ekaterina Smirnova
#'
#' @examples
#' data(mock2)
#' # Proportion data matrix
#' Prop <- mock2$Prop
#'
#' # Counts data matrix
#' Counts <- mock2$Counts
#'
#' #arrange counts in order of increasing number of samples taxa are present in
#' NC <- NC_Order(Counts)
#'
#' @export

NC_Order <- function(Counts){

  Netw <- t(Counts)%*%Counts
  Netw2 <- Netw
  diag(Netw2) <- 0
  NC_val <- sort(apply(Netw2, 2, nnzero), decreasing = FALSE)
  NC <- names(NC_val)
  return(NC)
}

###############################################################################
#' Taxa importance ordering by the weighted number of connected taxa
#'
#' @usage NCw_Order(Counts)
#'
#' @param Counts OTU COUNTS table, where taxa are columns and samples are rows of the table.
#' It should be a in data frame format with columns corresponding to taxa names.
#'
#' @return NCW Taxa names in increasing  order of the weighted number of connected taxa.
#'
#' @references Smirnova, E., Huzurbazar, H., Jafari, F. ``PERFect: permutation  filtration of microbiome data", to be submitted.
#'
#' @author Ekaterina Smirnova
#'
#' @examples
#' data(mock2)
#' # Proportion data matrix
#' Prop <- mock2$Prop
#'
#' # Counts data matrix
#' Counts <- mock2$Counts
#'
#' #arrange counts in order of increasing number of samples taxa are present in
#' NCw <- NCw_Order(Counts)
#'
#' @export

NCw_Order <- function(Counts){
  #first calculate NC values
  Netw <- t(Counts)%*%Counts
  Netw2 <- Netw
  diag(Netw2) <- 0
  NC_val <- sort(apply(Netw2, 2, nnzero), decreasing = FALSE)
  NC <- names(NC_val)
  #then weight by the number of samples taxon is present in
  n <- dim(Counts)[1]
  n_Tilde <- apply(Counts[,NC], 2, nnzero)
  NCW_val <- sort((n_Tilde*NC_val)/n, decreasing = FALSE)
  NCW <- names(NCW_val)
}

##############################################################
##############################################################
#' @export
filt_pval <- function(X, pvals, alpha, Order = "NP",   Order.user = NULL, pvals_sim = NULL){

  #Order columns by importance
  if(Order == "NP") {Order.vec <- NP_Order(X)}
  if(Order == "pvals") {Order.vec <- pvals_Order(X, pvals_sim)}
  if(Order == "NC"){Order.vec <- NC_Order(X)}
  if(Order == "NCw"){Order.vec <- NCw_Order(X)}
  else if (!is.null(Order.user)) {Order.vec = Order.user} #user-specified ordering of columns of X
  X <- X[,Order.vec]#properly order columns of X
  #remove all-zero OTU columns
  nzero.otu <- apply(X, 2, nnzero) != 0
  X <- X[, nzero.otu]
  p <- dim(X)[2]
  Order.vec <- Order.vec[nzero.otu]

  #select taxa that are kept in the data set at significance level alpha
  Ind <- which(pvals <=alpha)
  if (length(Ind !=0)) {Ind <- min(Ind)}
  else{Ind <- dim(X)[2]-1
  warning("no taxa are significant at a specified alpha level")}
  #if jth DFL is significant, then throw away all taxa 1:j
  filtX <- X[,-(1:Ind)]

  return(filtX)
}

