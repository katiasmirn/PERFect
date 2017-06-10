##############################
#define skew normal distribution to be used in fitristplus
#############################
require(sn);
require(fitdistrplus);
dsmnorm<<-function(x, mean = 0, sd = 1, xi = 1.5, log = FALSE){return(dsnorm(x, mean = 0, sd = 1, xi = 1.5, log = FALSE))}
psmnorm<<-function(q, mean = 0, sd = 1, xi = 1.5){return(psnorm(q, mean = 0, sd = 1, xi = 1.5))}
qsmnorm<<-function(p, mean = 0, sd = 1, xi = 1.5){return(qsnorm(p, mean = 0, sd = 1, xi = 1.5))}
rsmnorm<<-function(n, mean = 0, sd = 1, xi = 1.5){return(rsnorm(n, mean = 0, sd = 1, xi = 1.5))}

######################################
#Fit quantiles to the log data and get p-values
######################################
PERFect_sim <- function(X,  Order,  nbins =30, quant = c(0.25, 0.5), distr =c("norm", "sn"),
                    alpha = 0.05,
                    col = "red", fill = "green", hist_fill = 0.2, linecol = "blue"){

  pDFL <- NULL
  phist <- NULL
  X <- X[,Order]#properly order columns of X
  p <- dim(X)[2]
  #calculate DFL values
  Order_Ind <- rep(1:length(Order))#convert to numeric indicator values
  DFL <- DiffFiltLoss(X = X, Order_Ind, Plot = TRUE, Taxa_Names = Order) 
  #alternative calculation of filtering loss using presise formula
  #Function to calculate j^th DFL loss 
  
  Taxa <- Order[-length(Order)]
  lfl <- data.frame(Taxa, log(DFL$DFL))
  names(lfl) <- c("Taxa", "DFL")
  #plot histogram
  hist <- ggplot(data = lfl, aes(lfl$DFL)) + geom_histogram(bins = nbins, aes(y=..density..), 
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
    #calculate p-values
    pvals <- pnorm(q=lfl$DFL, mean = est[1], sd = est[2], lower.tail = FALSE, log.p = FALSE)
  }
  #estimate using skew normal 
  if(distr == "sn"){
    lp <- list(xi = mean(lfl$DFL), omega = sd(lfl$DFL), alpha = 1.5)
    fit <- qmedist(lfl$DFL, distr, probs=quant, start=lp)
    est <- fit$estimate
    hist <- hist + stat_function(fun = dsn, args = list(xi = est[1], omega = est[2], alpha = est[3]), colour=linecol)
    #calculate p-values
    pvals <- 1- psn(x=lfl$DFL, xi = est[1], omega = est[2], alpha = est[3])
  }
  
  #select taxa that are kept in the data set at significance level alpha
  names(pvals) <- names(DFL$DFL)
  Ind <- which(pvals <=alpha)
  if (length(Ind !=0)) {Ind <- min(Ind)}
  else{Ind <- dim(X)[2]-1
       warning("no taxa are significant at a specified alpha level")}
  #if jth DFL is significant, then throw away all taxa 1:j 
  filtX <- X[,-(1:Ind)]
  
  return(list(hist = hist, est = est, pvals = round(pvals,5), filtX = filtX, DFL = DFL$DFL, 
              pDFL = DFL$p + ylab("Difference in Filtering Loss"),  fit=fit)) 
}

