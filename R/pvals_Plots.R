pvals_Plots <- function(PERFect, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.05){
  pvals <- PERFect$pvals
  Order_pvals <- names(pvals)
  #taxa that separate the filtered out data set
  taxa_filt <- Order_pvals[!(Order_pvals %in% names(PERFect$filtX))]
  taxa <- max(which(taxa_filt %in% Order_pvals))
  #calculate FLu values
  res_FLu <- FiltLoss(X = Prop[,Order_pvals], Order = Order_pvals, type =  "Ind", Plot = TRUE)$FL
  #create a color scale for p-values plot
  FLu_vals <- res_FLu
  breaks <- quantile(res_FLu, quantiles)
  FLu_vals[FLu_vals <=breaks[1]] <- 1
  FLu_vals[FLu_vals > breaks[1] & FLu_vals <=breaks[2]] <- 2
  FLu_vals[FLu_vals > breaks[2] & FLu_vals <=breaks[3]] <- 3
  FLu_vals[FLu_vals > breaks[3] & FLu_vals <=breaks[4]] <- 4
  FLu_vals[!(FLu_vals %in% c(1,2,3,4))] <- 5

  FLu_vals[FLu_vals == 1] <- "0-25%"
  FLu_vals[FLu_vals == 2] <- "25 - 50%"
  FLu_vals[FLu_vals == 3] <- "50 - 80%"
  FLu_vals[FLu_vals ==4] <- "80 - 90%"
  FLu_vals[FLu_vals ==5] <- "90 - 100%"
  #overlay individual filtering loss information using points size
  #plot p-values
  df <- data.frame(seq(1:length(pvals)), pvals, res_FLu[which(names(res_FLu) %in% names(pvals))],
                 FLu_vals[which(names(FLu_vals) %in% names(pvals))])
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
  return(p_pvals)
}