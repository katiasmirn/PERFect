##########################################
#Use permutations to build the distribution 
#of DFL
###########################################
PERFect_perm <- function(X,  Order,   quant = c(0.1,0.25, 0.5), distr = "sn", alpha = 0.05,k = 10000,
                          nbins =30, col = "red", fill = "green", hist_fill = 0.2, linecol = "blue"){
  #X- OTU table with taxa in columns and samples in rows
  #Order - ordering of taxa
  #quant - quantiles used to fit distribution
  #distr - distribution to fit DFL to
    #1. sn - skew normal
    #2. normal - normal distribution
  #alpha - alpha level of the PERFect test
  #k- number of permutations
  #the rest of parameters are graphical parameters to construct the histogram
  #nbins - number of bins used in histogram
  #fill - color of histogram bars
  #hist_fill - color intensity of histogram bars
  #linecol - color of the line for fitted density
  p <- dim(X)[2]
  pvals <- rep(0,p-1)
  hist_list <- lapply(1:(p-1),function(x) NULL)
  est_list <- list()
  fit_list <- list()
  X <- X[,Order]
  #calculate loss
  #calculate DFL values
  Order_Ind <- rep(1:length(Order))#convert to numeric indicator values
  DFL <- DiffFiltLoss(X = X, Order_Ind, Plot = TRUE, Taxa_Names = Order) 
  #name p-values
  names(pvals) <- names(DFL$DFL)
  #For each taxon j, create a distribution of its DFL's by permuting the labels 
  dfl_distr <- sampl_distr(X = X, k=k)
  #build histograms for each taxon j
  for(i in 1:(p-1)) {
    
    #add a histogram 
    lfl <- data.frame(log(dfl_distr[[i]][!dfl_distr[[i]]==0]))
    if(length(dfl_distr[[i]][dfl_distr[[i]]==0])>0){
      print(paste("taxon", i, "number of zeroes = ", 
                  length(dfl_distr[[i]][dfl_distr[[i]]==0])))}
    names(lfl) <- c("DFL")
    #plot histogram
    x = "DFL"
    ord_map = aes_string(x = x)
    #hist <- ggplot(lfl, ord_map) + geom_histogram()
    hist <- ggplot(lfl, ord_map) + geom_histogram(bins = nbins, aes(y=..density..), 
                                                  col = col, fill = fill, alpha =hist_fill)+
      theme(panel.background = element_rect(fill = "white"), 
            panel.grid.major = element_line(colour = "grey90"),
            axis.text.x  = element_text( size=10))+
      ggtitle("") + xlab("Filtering Loss") + ylab("Density")
    
    #estimate using normal
    if(distr == "norm"){
      fit <- qmedist(lfl$DFL, distr, probs=quant)  
      est <- fit$estimate
      #add density line to the plot
      hist <- hist + stat_function(fun = dnorm, args = list(mean = est[1], sd = est[2]), colour=linecol)
      hist_list[[i]] <- hist
      #calculate p-values
      pvals[i] <- pnorm(q=log(DFL$DFL[i]), mean = est[1], sd = est[2], lower.tail = FALSE, log.p = FALSE)
    }
    
    #estimate using skew normal 
    if(distr == "sn"){
      lp <- list(xi = mean(lfl$DFL), omega = sd(lfl$DFL), alpha = 1.5)
      fit <- qmedist(lfl$DFL, distr, probs=quant, start=lp)
      est <- fit$estimate
      hist <- hist + stat_function(fun = dsn, args = list(xi = est[1], omega = est[2], alpha = est[3]), colour=linecol)
      hist_list[[i]] <- hist
      #calculate p-values
      pvals[i] <- 1- psn(x=log(DFL$DFL[i]), xi = est[1], omega = est[2], alpha = est[3])
    }
    #save estimate results
    est_list[[i]] <- est
    fit_list[[i]] <- fit
  }
  
  #select taxa that are kept in the data set at significance level alpha
  Ind <- which(pvals <=alpha)
  if (length(Ind !=0)) {Ind <- min(Ind)}
  else{Ind <- dim(X)[2]-1
       warning("no taxa are significant at a specified alpha level")}
  #if jth DFL is significant, then throw away all taxa 1:j 
  filtX <- X[,-(1:Ind)]
  return(list(filtX = filtX, pvals = pvals, fit = fit_list, hist = hist_list, 
              est =est_list,   dfl_distr=dfl_distr ))
}

