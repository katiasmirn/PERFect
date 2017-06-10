PERFect_perm_reorder <- function(X,  Order_alt,  res_perm, alpha = 0.05, distr = "sn"){
  
  X2 <- X
  X2 <- X2[,Order_alt]
  Order_Ind <- rep(1:length(Order))#convert to numeric indicator values
  DFL <- DiffFiltLoss(X = X2, Order_Ind, Plot = TRUE, Taxa_Names = Order_alt) 
  #re-evaluate p-values
  pvals <- rep(0, length(DFL$DFL))
  names(pvals) <- names(DFL$DFL)
  for (i in 1:length(DFL$DFL)){
    if(distr=="sn"){
      pvals[i] <- 1- psn(x=log(DFL$DFL[i]), 
                     xi = res_perm$est[[i]][1], omega = res_perm$est[[i]][2], 
                     alpha = res_perm$est[[i]][3])
    }
    if(distr == "norm"){
      #calculate p-values
      pvals[i] <- pnorm(q=log(DFL$DFL[i]), mean = res_perm$est[[i]][1], sd = res_perm$est[[i]][2], 
                        lower.tail = FALSE, log.p = FALSE)
    }
    
    } 
  
  #re-calculate filtered X
  #select taxa that are kept in the data set at significance level alpha
  Ind <- which(pvals <=alpha)
  if (length(Ind !=0)) {Ind <- min(Ind)}
  else{Ind <- dim(X2)[2]-1
     warning("no taxa are significant at a specified alpha level")}
  #if jth DFL is significant, then throw away all taxa 1:j 
  res_perm$filtX <- X2[,-(1:Ind)]
  res_perm$pvals <- pvals #end if !is.null(res_perm)

  return(res_perm)
}