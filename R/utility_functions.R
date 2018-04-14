########################################
#Function to calculate j^th DFL loss
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
ResultsComparison <-function(X, Counts,Ab_min = 0.001,rel =FALSE, thresh=5,prop =FALSE,
                             res_sim, res_perm){
  
  #Traditional filtering
  Filt_R1 <- TraditR1(X=Counts,  rel =rel, thresh=thresh)#if rel = TRUE then use propo of samples > thresh (say 0.05 = 5%)
  Filt_R2 <- TraditR2(X=Counts,  Ab_min = Ab_min, prop =prop)#prop =FALSE if counts matrix is used
  #taxa left after traditional rules 1 and 2 applied
  Trad_R1 <- names(Filt_R1)
  Trad_R2 <- names(Filt_R2)
  
  #taxa left after traditional rules 1 and 2 applied
  Loss_Trad_R1 <- FL_J(X = X, J = Trad_R1, leave = TRUE)  
  Loss_Trad_R2 <- FL_J(X = X, J = Trad_R2, leave = TRUE)
  
  #PERFect loss
  PERFect_sim_taxa <- names(res_sim$filtX)
  PERFect_perm_taxa <- names(res_perm$filtX)
  #corresponding filtering loss
  Loss_PERFect_sim <- FL_J(X = X, J = PERFect_sim_taxa, leave = TRUE)  
  Loss_PERFect_perm <- FL_J(X = X, J = PERFect_perm_taxa, leave = TRUE)
  
  #table for the paper
  Res <- cbind(rbind(length(Trad_R1), length(Trad_R2), length(PERFect_sim_taxa), length(PERFect_perm_taxa)), 
               rbind(Loss_Trad_R1, Loss_Trad_R2, Loss_PERFect_sim, Loss_PERFect_perm))
  colnames(Res) <- c("N_taxa", "FL")
  rownames(Res) <- c("Trad_R1", "Trad_R2", "PERFect_sim", "PERFect_perm")
  return(list(Res = Res,  Filt_R1 = Filt_R1, Filt_R2 = Filt_R2, res_sim = res_sim, res_perm = res_perm))
}
#####################
#Order functions
####################
NP_Order <- function(Counts){
  #arrange counts in order of increasing number of samples taxa are present in 
  NP <- names(sort(apply(Counts, 2, nnzero)))
  return(NP)
}

pvals_Order <- function(Counts, res_sim){
  vec <- names(res_sim$pvals)
  Order <- c(setdiff(colnames(Counts), vec), vec)
  Order_pvals <- Order[c(1,sort.int(res_sim$pvals, index.return = TRUE, decreasing = TRUE)$ix +1)]
  return(Order_pvals)
}

NC_Order <- function(Counts){
  
  Netw <- t(Counts)%*%Counts
  Netw2 <- Netw
  diag(Netw2) <- 0
  NC_val <- sort(apply(Netw2, 2, nnzero), decreasing = FALSE)
  NC <- names(NC_val)
  return(NC)
}

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

