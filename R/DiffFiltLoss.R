DiffFiltLoss <- function(X,  Order_Ind, Plot = TRUE, Taxa_Names = NULL){
  #X - data matrix with taxa in columns and samples in rows
  p <- dim(X)[2]#total #of taxa
  DFL <- rep(1,p-1)
  X <- as.matrix(X)
  Netw <- t(X)%*%X

  for(j in 1:(p-1)){DFL[j] <- DiffFiltLoss_j(Order_Ind,Netw, j)}
  DFL <- DFL/tr(t(Netw)%*%Netw)
  
  if(Plot == TRUE){
  
  #Plot Full Norm reduction 
  df <- data.frame(Order_Ind[-1], rep(1:length(DFL)), DFL)
  names(df)[2] <- "x"
  Lab <- 1:length(Order_Ind[-1])
  df <- cbind(Lab, df)
  
  #Plots      
  p_FL <- ggplot(df) + geom_line(aes(x = Lab, y = DFL, group =1), 
                                 colour = "dodgerblue3")+
    theme(panel.background = element_rect(fill = "white"), 
          panel.grid.major = element_line(colour = "grey90"),
          axis.text.x  = element_text( size=10,colour="black", angle = 90, hjust = 1))+
    ggtitle("") + ylab("Differences in Filtering Loss") + xlab("Taxa")+xlim(0, max(df$Lab))
  
  
}#end if Plot = TRUE

if(!is.null(Taxa_Names)) {names(DFL) <- Taxa_Names[-1]}


return(list(DFL = DFL, p_FL= p_FL))
}