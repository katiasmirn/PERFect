#' Difference in filtering loss
#'
#' @description This function calculates differences in filtering loss due to removing a set of J taxa sequentially.
#'
#' @usage DiffFiltLoss(X, Order_Ind, Plot = TRUE, Taxa_Names = NULL)
#'
#' @param X OTU table, where taxa are columns and samples are rows of the table. It should be a in dataframe format
#' with columns corresponding to taxa names.
#'
#' @param Order_Ind Numeric column order corresponding to taxa importance arrangement.
#'
#' @param Plot A binary TRUE/FALSE value. If TRUE, the function returns plot of sequential differences in filtering loss.
#'
#' @param Taxa_Names Optional taxa labels corresponding to the columns ordering given by Order_Ind.
#'
#' @details This function calculates and plots (if Plot = TRUE) differences in filtering loss sequentially
#' for removing the first j taxa as DFL(j+1) = FL(J_{j+1}) - FL(J_j) for taxa j=1, ..., p.
#'
#' @return
#' \item{DFL}{Differences in filtering loss values}
#' \item{p_FL}{Plot of the differences in filtering loss}
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
#'
#' #obtain numeric column order corresponding to taxa importance arrangment
#' Order_Ind <- match(NP, names(Prop))
#' DFL <- DiffFiltLoss(X=Prop, Order_Ind = Order_Ind, Plot = TRUE, Taxa_Names = NP)
#'
#' #Differences in filtering loss values
#' DFL$DFL
#'
#' #Plot of the differences in filtering loss
#' DFL$p_FL
#' @export

DiffFiltLoss <- function(X,  Order_Ind, Plot = TRUE, Taxa_Names = NULL){

  # Check the format of X
  if(!(class(X) %in% c("matrix"))){X <- as.matrix(X)}
  #   stop('X must be a data frame or a matrix')
  # if(!(class(X) == "matrix")){X <- as.matrix(X)}

  # Check the format of Order_Ind
  if(class(Order_Ind) != "integer")
    stop('Order_Ind argument must be a vector of integer values')

  # Check the format of Plot
  if(class(Plot) != "logical")
    stop('Plot argument must be a logical value')

  #X - data matrix with taxa in columns and samples in rows
  p <- dim(X)[2]#total #of taxa
  DFL <- rep(1,p-1)
  X <- as.matrix(X)
  Netw <- t(X)%*%X

  for(j in 1:(p-1)){DFL[j] <- DiffFiltLoss_j(Order_Ind,Netw, j)}
  DFL <- DFL/sum(Netw*Netw)

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
