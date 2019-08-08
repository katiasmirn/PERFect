#' Plots of PERFect p-values
#'
#' @description Graphical representation of p-values obtained by running \code{PERFect_sim()} or
#' \code{PERFect_perm()} for jth taxon colored by quantile values of individual filtering loss.
#'
#' @usage pvals_Plots(PERFect, X, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha = 0.1)
#'
#' @param PERFect Output of \code{PERFect_sim()} or \code{PERFect_perm()} function.
#'
#' @param X OTU table, where taxa are columns and samples are rows of the table.
#' It should be a in data frame format with columns corresponding to taxa names.
#' @param quantiles Quantile values for coloring, these are set to 25\%, 50\%, 80\% and
#' 90\% percentiles of the individual filtering loss values.
#' @param alpha Alpha level of the test, set to 0.1 by default.
#'
#' @return p_vals Plot of p-values
#'
#' @author Ekaterina Smirnova
#'
#' @seealso \code{\link{PERFect_sim}}, \code{\link{PERFect_perm}}
#'
#' @examples
#' data(mock2)
#' # Proportion data matrix
#' Prop <- mock2$Prop
#'
#' # Counts data matrix
#' Counts <- mock2$Counts
#' dim(Counts) # 240x46
#'
#' # Perform simultaenous filtering of the data
#' res_sim <- PERFect_sim(X=Counts)
#' dim(res_sim$filtX)      # 240x10, removing 36 taxa
#' colnames(res_sim$filtX) # signal taxa
#'
#' #permutation perfect colored by FLu values
#' pvals_Plots(PERFect = res_sim, X = Counts, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.05)
#' @export
pvals_Plots <- function(PERFect, X, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.1){

  # Check if PERFect object is input correctly
  if(length(PERFect$pvals) == 0){
    stop('PERFect object must be a result from PERFect_sim() or PERFect_perm()')
  }

  # Check the format of X
  if(!(class(X) %in% c("matrix"))){X <- as.matrix(X)}
  #   stop('X must be a data frame or a matrix')
  # if(!(class(X) == "matrix")){X <- as.matrix(X)}

  # Check the format of quantiles
  if(!is.vector(quantiles)) stop('quantiles argument must be a vector')

  # Check the format of alpha
  if(!is.numeric(alpha)) stop('alpha argument must be a numerical value')

  if(!(0 %in% quantiles)){quantiles <- c(0, quantiles)}
  if(!(1 %in% quantiles)){quantiles <- c(1, quantiles)}
  quantiles <- sort(quantiles)

  pvals <- PERFect$pvals
  Order_pvals <- names(pvals)
  #taxa that separate the filtered out data set
  taxa_filt <- Order_pvals[!(Order_pvals %in% colnames(PERFect$filtX))]
  taxa <- max(which(taxa_filt %in% Order_pvals))
  #calculate FLu values
  res_FLu <- FiltLoss(X = X[,Order_pvals], Order.user = Order_pvals, type =  "Ind", Plot = TRUE)$FL
  #create a color scale for p-values plot
  FLu_vals <- res_FLu
  breaks <- quantile(res_FLu, quantiles)
  seq1 <- round(100*quantiles)[-length(quantiles)]
  seq2 <- round(100*quantiles)[-1]
  labs <- paste(seq1, seq2, sep = "-")
  labs <- paste("[",labs, ")%", sep = "")
  FLu_vals <- cut(res_FLu, breaks, right=FALSE , labels = labs, include.lowest = TRUE)

  #overlay individual filtering loss information using points size
  #plot p-values
  Ind <- which(names(res_FLu) %in% names(pvals))
  df <- data.frame(seq(1:length(pvals)),
                   pvals,
                   res_FLu[Ind],
                   FLu_vals[Ind])
  names(df) <- c("Taxa", "p_value", "FLu", "Quantiles")
  p_pvals <- ggplot(df) + geom_point( aes(x = Taxa, y = p_value, color = Quantiles)) +
  ggtitle("Permutation PERFect p-values") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey90"),
        axis.text.x  = element_text( size=10,colour="black", angle = 90, hjust = 1))+
  guides(color=guide_legend(title="FLu Quantiles"))
  #add alpha level horizontal line info

  p_pvals <- p_pvals + geom_hline(yintercept=alpha, color="red", linetype="dashed")
  #add taxa cutoff at alpha level vertical  line info
  p_pvals <- p_pvals + ggtitle("") + geom_vline(xintercept=taxa, color="purple", linetype="dashed")
  return(list(data = df, plot = p_pvals))
}
