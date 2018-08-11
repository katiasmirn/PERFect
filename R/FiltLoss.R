############################
#filtering loss 
############################
FiltLoss <- function(X, Order = "NP", Order.user = NULL, type = "Cumu", Plot = TRUE){
  #X - data matrix
  #type = Ind - only jth taxon removed to calculate FL
  #type = Cm - 1:j tax are removed to calculate FL
  p <- dim(X)[2]#total #of taxa
  Norm_Ratio <- rep(1, p)
  X <- as.matrix(X)
  
  #Order columns by importance
  if(Order == "NP") {Order.vec <- NP_Order(X)}
  if(Order == "pvals") {Order.vec <- pvals_Order(X, pvals_sim)}
  if(Order == "NC"){Order.vec <- NC_Order(X)}
  if(Order == "NCw"){Order.vec <- NCw_Order(X)}
  else if (!is.null(Order.user)) {Order.vec = Order.user} #user-specified ordering of columns of X
  X <- X[,Order.vec]#properly order columns of X
  
  Order_Ind <- 1:length(Order.vec) 
  Netw <- t(X)%*%X
  
  #Taxa at the top of the list have smallest number of connected nodes
  for (i in 1:p){
    if (type == "Cumu") {Ind <- Order_Ind[-(1:i)]}
    else{Ind <- Order_Ind[-i]}
    
    #define matrix X_{-J}'X_{-J} for the choice of cumulative or individual filtering loss
    Netw_R <- Netw[Ind, Ind]
    #calculate the corresponding norm
    Norm_Ratio[i] <-  tr(t(Netw_R)%*%Netw_R)

  }#end for
  
  FL <- 1 - Norm_Ratio/tr(t(Netw)%*%Netw) #result divided by the full matrix norm 
  if(Plot == TRUE){
    
    #Plot Full Norm reduction 
    df <- data.frame(Order.vec, rep(1:length(FL)), FL)
    names(df)[2] <- "x"
    Lab <- 1:length(Order.vec)
    df <- cbind(Lab, df)
    
    #Plots      
    p_FL <- ggplot(df) + geom_line(aes(x = Lab, y = FL, group =1), 
                                   colour = "dodgerblue3")+
      theme(panel.background = element_rect(fill = "white"), 
            panel.grid.major = element_line(colour = "grey90"),
            axis.text.x  = element_text( size=10,colour="black", angle = 90, hjust = 1))+
      ggtitle("") + ylab("Filtering Loss") + xlab("Taxa")+xlim(0, max(df$Lab))
    
    
  }#end if Plot = TRUE
  names(FL) <-  colnames(X)
  
  return(list(FL = FL, p_FL= p_FL))
}

