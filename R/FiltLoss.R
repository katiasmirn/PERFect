############################
#filtering loss
############################

#' Filtering Loss
#'
#' @description Sequential filtering loss calculation for removing a set of J_j taxa for J= 1, ..., p.
#'
#' @usage FiltLoss(X, Order = "NP", Order.user = NULL, type = "Cumu", Plot = TRUE)
#'
#' @param X OTU table, where taxa are columns and samples are rows of the table.
#' It should be a in data frame format with columns corresponding to taxa names.
#'
#' @param Order Taxa ordering. The default ordering is the number of occurrences (NP) of the taxa in all samples.
#'  Other types of order are number of connected taxa and weighted number of connected taxa,
#'  denoted as \code{"NC"}, \code{"NCw"} respectively. More details about taxa ordering are described in Smirnova et al.
#'  User can also specify their preference order with Order.user.
#'
#' @param Order.user User's taxa ordering. This argument takes a character vector of ordered taxa names.
#'
#' @param type Type of filtering loss calculation.
#' \describe{
#' \item{\code{"Ind"}}{Individual taxon's filtering loss FL_u(j)}
#' \item{\code{"Cumu"}}{Cumulative filtering loss FL(J) due to removing a set of taxa J}
#' }
#'
#' @param Plot Binary TRUE/FALSE value. If TRUE, the function returns plot of sequential differences in filtering loss.
#'
#' @details
#'
#' The individual filtering loss due to removing one taxon j is defined as:
#'
#' FL_u(j)= 1- (||X^T_{-j} X_{-j}||_F^2/||X^TX||_F^2),
#'
#' where X_{-j} is the matrix X without column corresponding to jth taxon and ||Z||_F is the Frobenious norm of a matrix Z.
#'
#' The cumulative filtering loss due to removing a set of taxa is defined as:
#'
#' FL(J)= 1- (||X^T_{-J}  X_{-J}||_F^2\||X^TX||_F^2),
#'
#' where X_{-J} is the n x (p-|J|) dimensional matrix obtained by removing the columns indexed by the set J from the data matrix X.
#'
#' The cumulative filtering loss is calculated sequentially for each set of taxa J_j, j=1, ..., p.
#'
#' @return
#' \item{FL}{Filtering loss values}
#' \item{p_FL}{Plot of filtering loss values}
#'
#' @references Smirnova, E., Huzurbazar, H., Jafari, F. ``PERFect: permutation  filtration of microbiome data", to be submitted.
#'
#' @author Ekaterina Smirnova
#'
#' @seealso \code{\link{DiffFiltLoss}}
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
#' #Calculate cumulative filtering loss
#' FL <- FiltLoss(X=Prop, Order = "NP", type = "Cumu", Plot = TRUE)
#'
#' #Differences in filtering loss values
#' FL$FL
#'
#' #Plot of the differences in filtering loss
#' FL$p_FL
#' @export
FiltLoss <- function(X, Order = "NP", Order.user = NULL, type = "Cumu", Plot = TRUE){
  #X - data matrix
  #type = Ind - only jth taxon removed to calculate FL
  #type = Cm - 1:j tax are removed to calculate FL
  p <- dim(X)[2]#total #of taxa
  Norm_Ratio <- rep(1, p)

  # Check the format of X
  if(!(class(X) %in% c("matrix"))){X <- as.matrix(X)}
  #   stop('X must be a data frame or a matrix')
  # if(!(class(X) == "matrix")){X <- as.matrix(X)}

  # Check the format of Order
  if(!(Order %in% c("NP","NC","NCw")))
    stop('Order argument can only be "NP", "NC", or "NCw" ')

  # Check the format of type
  if(!(type %in% c("Cumu","Ind")))
    stop('type argument can only be "Cumu" or "Ind" ')

  # Check the format of Plot
  if(class(Plot) != "logical")
    stop('Plot argument must be a logical value')

  #Order columns by importance
  if(is.null(Order.user)){
    if(Order == "NP") {Order.vec <- NP_Order(X)}
    #if(Order == "pvals") {Order.vec <- pvals_Order(X, pvals_sim)}
    if(Order == "NC"){Order.vec <- NC_Order(X)}
    if(Order == "NCw"){Order.vec <- NCw_Order(X)}
  } else {
    Order.vec <- Order.user #user-specified ordering of columns of X
  }
  X <- X[,Order.vec]#properly order columns of X

  Order_Ind <- seq_len(length(Order.vec))
  Netw <- t(X)%*%X

  #Taxa at the top of the list have smallest number of connected nodes
  for (i in seq_len(p)){
    if (type == "Cumu") {Ind <- Order_Ind[-seq_len(i)]}
    else{Ind <- Order_Ind[-i]}

    #define matrix X_{-J}'X_{-J} for the choice of cumulative or individual filtering loss
    Netw_R <- Netw[Ind, Ind]
    #calculate the corresponding norm
    #Norm_Ratio[i] <-  psych::tr(t(Netw_R)%*%Netw_R)
    Norm_Ratio[i] <-  sum(Netw_R*Netw_R)
  }#end for

  #FL <- 1 - Norm_Ratio/psych::tr(t(Netw)%*%Netw) #result divided by the full matrix norm
  FL <- 1 - Norm_Ratio/sum(Netw*Netw)
  if(Plot == TRUE){

    #Plot Full Norm reduction
    df <- data.frame(Order.vec, rep(seq_len(length(FL))), FL)
    names(df)[2] <- "x"
    Lab <- seq_len(length(Order.vec))
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

