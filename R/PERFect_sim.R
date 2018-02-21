##############################
#define skew normal distribution to be used in fitristplus
#############################
require(sn);
require(fitdistrplus);
dsmnorm<<-function(x, mean = 0, sd = 1, xi = 1.5, log = FALSE){return(dsnorm(x, mean = 0, sd = 1, xi = 1.5, log = FALSE))}
psmnorm<<-function(q, mean = 0, sd = 1, xi = 1.5){return(psnorm(q, mean = 0, sd = 1, xi = 1.5))}
qsmnorm<<-function(p, mean = 0, sd = 1, xi = 1.5){return(qsnorm(p, mean = 0, sd = 1, xi = 1.5))}
rsmnorm<<-function(n, mean = 0, sd = 1, xi = 1.5){return(rsnorm(n, mean = 0, sd = 1, xi = 1.5))}

#available distributions: 
#normal - "norm"
#skew-normal - "sn"
#non-central t - "t"
#non-central beta - "cauchy"
######################################
#Fit quantiles to the log data and get p-values
######################################
PERFect_sim <- function(X,  Order = "NP",   Order.user = NULL, 
                    normalize = "counts", center = FALSE,
                    quant = c(0.10, 0.25, 0.5),  distr ="sn",
                    alpha = 0.10, lag = 3, direction ="left",
                    pvals_sim = NULL,
                    nbins =30, 
                    col = "red", fill = "green", hist_fill = 0.2, linecol = "blue"){

  pDFL <- NULL
  phist <- NULL
  
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
  
  #save non-centered, non-normalized X
  X.orig <- X
  
  #normalize the data
  if(normalize == "prop"){X <- X/apply(X, 1, sum)}
  else if (normalize == "pres"){X[X!=0]<-1}
  
  #center if true
  if(center){X <- apply(X, 2, function(x) {x-mean(x)})}
  
  #calculate DFL values
  Order_Ind <- rep(1:length(Order.vec))#convert to numeric indicator values
  DFL <- DiffFiltLoss(X = X, Order_Ind, Plot = TRUE, Taxa_Names = Order.vec) 
  #alternative calculation of filtering loss using presise formula
  #Function to calculate j^th DFL loss 
  
  Taxa <- Order.vec[-length(Order.vec)]
  lfl <- data.frame(Taxa, log(DFL$DFL))
  names(lfl) <- c("Taxa", "DFL")
  #plot histogram
  hist <- ggplot(data = lfl, aes(lfl$DFL)) + geom_histogram(bins = nbins, aes(y=..density..), 
                                                           col = col, fill = fill, alpha =hist_fill)+
    theme(panel.background = element_rect(fill = "white"), 
          panel.grid.major = element_line(colour = "grey90"),
          axis.text.x  = element_text( size=10))+
    ggtitle("") + xlab("log differences in filtering loss") + ylab("Density")
  #estimate using normal 
  if(distr == "norm"){
    if(length(quant) > 2){quant <- quant[(length(quant) - 1):length(quant)]
                          print("Warning: more than 2 quantile values are given. \nLargest 2 quantiles are used.")}
    if(length(quant) < 2){stop("At least two quantile values must be specified.")}
    fit <- qmedist(lfl$DFL, distr, probs=quant)
    est <- fit$estimate
    #add density line to the plot
    hist <- hist + stat_function(fun = dnorm, args = list(mean = est[1], sd = est[2]), colour=linecol)
    #calculate p-values
    pvals <- pnorm(q=lfl$DFL, mean = est[1], sd = est[2], lower.tail = FALSE, log.p = FALSE)
  }
  #estimate using t-distribution 
  if(distr == "t"){
      if(length(quant) > 2){quant <- quant[(length(quant) - 1):length(quant)]
                            print("Warning: more than 2 quantile values are given. \nLargest 2  quantiles are used.")}
      if(length(quant) < 2){stop("At least 2 quantile value must be specified.")}
    fit <- qmedist(lfl$DFL, distr, probs=quant, start=list(df=2, ncp = mean(lfl$DFL)))
    est <- fit$estimate
    #add density line to the plot
    hist <- hist + stat_function(fun = dt, args = list(df =est[1],  ncp = est[2]), colour=linecol)
    #calculate p-values
    pvals <- pt(q=lfl$DFL,  df =est[1],  ncp = est[2],  lower.tail = FALSE, log.p = FALSE)
  }
    #estimate using cauchy distribution 
    if(distr == "cauchy"){
      if(length(quant) > 2){quant <- quant[(length(quant) - 1):length(quant)]
                        print("Warning: more than 2 quantile values are given. \nLargest 2 quantiles are used.")}
      if(length(quant) < 2){stop("At least 2 quantile value must be specified.")}
    fit <- qmedist(lfl$DFL, distr, probs=quant)
    est <- fit$estimate
    #add density line to the plot
    hist <- hist + stat_function(fun = dcauchy, args = list(location = est[1],  scale= est[2]), colour=linecol)
    #calculate p-values
    pvals <- pcauchy(q=lfl$DFL,  location =est[1],  scale = est[2],  lower.tail = FALSE, log.p = FALSE)
  }
  #estimate using skew normal 
  if(distr == "sn"){
    if(length(quant) > 3){quant <- quant[(length(quant) - 2):length(quant)]
                          print("Warning: more than 3 quantile values are given. \nLargest 3 quantiles are used.")}
    if(length(quant) < 3){stop("At least 3 quantile values must be specified.")}
    lp <- list(xi = mean(lfl$DFL), omega = sd(lfl$DFL), alpha = 1.5)
    fit <- qmedist(lfl$DFL, distr, probs=quant, start=lp)
    est <- fit$estimate
    hist <- hist + stat_function(fun = dsn, args = list(xi = est[1], omega = est[2], alpha = est[3]), colour=linecol)
    #calculate p-values
    pvals <- 1- psn(x=lfl$DFL, xi = est[1], omega = est[2], alpha = est[3])
  }
  
  #select taxa that are kept in the data set at significance level alpha
  names(pvals) <- names(DFL$DFL)
  
  #smooth p-values
  pvals_avg <- rollmean(pvals, k=lag, align=direction,  fill=NA )
  #replace na's with original values
  pvals_avg[is.na(pvals_avg)] <- pvals[is.na(pvals_avg)]
  
  Ind <- which(pvals_avg <=alpha)
  if (length(Ind !=0)) {Ind <- min(Ind)}
  else{Ind <- dim(X)[2]-1
       warning("no taxa are significant at a specified alpha level")}
  #if jth DFL is significant, then throw away all taxa 1:j 
  filtX <- X.orig[,-(1:Ind)]
  
  return(list(filtX = filtX, pvals = round(pvals_avg,5), DFL = DFL$DFL, fit=fit, hist = hist, est = est,   
              pDFL = DFL$p + ylab("Difference in Filtering Loss"))) 
}

