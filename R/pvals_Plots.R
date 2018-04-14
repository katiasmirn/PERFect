pvals_Plots <- function(PERFect, X, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.1){
  
  if(!(0 %in% quantiles)){quantiles <- c(0, quantiles)}
  if(!(1 %in% quantiles)){quantiles <- c(1, quantiles)}
  quantiles <- sort(quantiles)
  
  pvals <- PERFect$pvals
  Order_pvals <- names(pvals)
  #taxa that separate the filtered out data set
  taxa_filt <- Order_pvals[!(Order_pvals %in% colnames(PERFect$filtX))]
  taxa <- max(which(taxa_filt %in% Order_pvals))
  #calculate FLu values
  res_FLu <- FiltLoss(X = X[,Order_pvals], Order = Order_pvals, type =  "Ind", Plot = TRUE)$FL
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