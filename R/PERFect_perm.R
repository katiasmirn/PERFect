##########################################
#Use permutations to build the distribution
#of DFL
###########################################
#' Permutation PERFect filtering for microbiome data
#'
#' @description Permutation filtering of the provided OTU table X at a test level alpha.
#'   Each set of j taxa significance is evaluated by fitting the Skew-Normal, Normal, t or Cauchy distribution
#'   to the  sampling distribution obtained by permuted taxa labels.
#'
#' @usage PERFect_perm(X, infocol = NULL, Order = "NP", Order.user = NULL, normalize = "counts",
#'     algorithm = "fast", center = FALSE, quant = c(0.1, 0.25, 0.5), distr = "sn",
#'     alpha = 0.1, lag = 3, direction = "left", pvals_sim = NULL,
#'     k = 10000, dfl_distr = NULL, nbins = 30, hist = FALSE, col = "red",
#'     fill = "green", hist_fill = 0.2, linecol = "blue")
#'
#' @param X OTU table, where taxa are columns and samples are rows of the table.
#' It should be a in data frame format with columns corresponding to taxa names.
#' 
#' @param infocol Index vector of the metadata. We assume user only gives a taxa table,
#' but if the metadata of the samples are included in the columns of the input, this option
#' needs to be specified.
#'
#' @param Order Taxa ordering. The default ordering is the number of occurrences (NP) of the taxa in all samples.
#'  Other types of order are p-value ordering, number of connected taxa and weighted number of connected taxa,
#'  denoted as \code{"pvals"}, \code{"NC"}, \code{"NCw"} respectively. More details about taxa ordering are described in Smirnova et al.
#'  User can also specify their preference order with Order.user.
#'
#' @param normalize Normalizing taxa count. The default option does not normalize taxa count,
#'  but user can convert the OTU table into a proportion table using the option \code{"prop"}
#'  or convert it into a presence/absence table using \code{"pres"}.
#'
#' @param algorithm Algorithm speed. The default is speed is \code{"fast"}, which allows the program to efficiently
#' search for significant taxa without computing all the p-values. User must use the default option \code{"hist = FALSE"}
#' for the fast algorithm. The alternative setting is \code{"full"}, which computes all the taxa's p-values.
#'
#' @param center Centering OTU table. The default option does not center the OTU table.
#'
#' @param quant Quantile values used to fit the distribution to log DFL values.
#'  The number of quantile values corresponds to the number of parameters in the distribution the data is fitted to.
#'  Assuming that at least 50\% of taxa  are not informative, we suggest fitting the log Skew-Normal distribution
#'  by matching the 10\%, 25\% and 50\% percentiles of the log-transformed samples to the Skew-Normal distribution.
#'
#' @param distr The type of distribution to fit log DFL values to. While we suggest using Skew-Normal distribution,
#' and set as the default distribution, other choices are available.
#' \describe{
#' \item{\code{"sn"}}{Skew-Normal distribution with 3 parameters: location xi, scale omega^2 and shape alpha}
#' \item{\code{"norm"}}{Normal distribution with 2 parameters: mean and standard deviation sd}
#' }
#' @param alpha Test level alpha, set to 0.1 by default.
#'
#' @param lag Integer width of the rolling window in rolling average (moving mean), set to 3 by default.
#'
#' @param direction Character specifying whether the index of the result should be left- or right-aligned
#'  or centered compared to the rolling window of observations, set to "left" by default.
#'
#' @param pvals_sim Object resulting from simultaneous PERFect with taxa abundance ordering,
#'  allowing user to perform Simultaneous PERFect with p-values ordering.
#'  Be aware that the choice of distribution for both methods must be the same.
#'
#' @param k The number of permutations, set to 10000 by default.
#'
#' @param nbins Number of bins used to visualize the histogram  of log DFL values, set to 30 by default.
#'
#' @param col Graphical parameter for  color of histogram bars border, set to "red" by default.
#'
#' @param fill Graphical parameter for color of histogram fill, set to "green" by default.
#'
#' @param hist_fill Graphical parameter for intensity of histogram fill, set to 0.2 by default.
#'
#' @param linecol Graphical parameter for the color of the fitted distribution density, set to "blue" by default.
#'
#' @details Filtering is the process of identifying and removing a subset of taxa according to a particular criterion.
#'   As opposed to the the simultaneous filtering approach, we do not assume that all distributions
#'   for each set of taxa  are identical and equal to the distribution of simultaneous filtering.
#'   Function PERFect_perm()  filters  the provided OTU table X and outputs a filtered table that contains signal taxa.
#'   PERFect_perm() calculates differences in filtering loss DFL for each taxon according to the given taxa order.
#'   By default, the function fits Skew-Normal distribution to the log-differences in filtering loss but Normal,
#'   t, or Cauchy distributions can be also used.
#'
#' @return
#'
#' \item{filtX}{Filtered OTU table}
#'
#' \item{pvals}{P-values of the test}
#'
#' \item{DFL}{Differences in filtering loss values}
#'
#' \item{fit}{Fitted values and further goodness of fit details passed from the \code{fitdistr()} function}
#'
#' \item{hist}{Histogram of log differences in filtering loss}
#'
#' \item{est}{Estimated distribution parameters}
#'
#' \item{dfl_distr}{Plot of differences in filtering loss values}
#'
#' @references Azzalini, A. (2005). The skew-normal distribution and related multivariate families. Scandinavian Journal of Statistics, 32(2), 159-188.
#'
#' @references Smirnova, E., Huzurbazar, H., Jafari, F. ``PERFect: permutationfiltration of microbiome data", to be submitted.
#'
#' @author Ekaterina Smirnova
#'
#' @seealso \code{\link{PERFect_sim}}
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
#' # Perform simultaenous filtering of the data
#' res_sim <- PERFect_sim(X=Counts)
#'
#' #order according to p-values
#' pvals_sim <- pvals_Order(Counts, res_sim)
#' system.time(res_perm <- PERFect_perm(X = Prop, Order.user = pvals_sim))
#'
#' #permutation perfect colored by FLu values
#' pvals_Plots(PERFect = res_perm, X = Counts, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.05)
#' @export

# Algorithm with parallel processing
PERFect_perm <- function(X, infocol = NULL, Order = "NP",   Order.user = NULL,
                         normalize = "counts", algorithm = "fast", center = FALSE,
                         quant = c(0.10, 0.25, 0.5),  distr ="sn",
                         alpha = 0.10, lag = 3, direction ="left",
                         pvals_sim = NULL,
                         k=10000, dfl_distr = NULL,
                         nbins =30,  hist = TRUE,
                         col = "red", fill = "green", hist_fill = 0.2, linecol = "blue"){
  
  info <- NULL
  #infocol = index vector of other info
  if(!is.null(infocol)){
    info <- X[,infocol]
    X <- X[,-infocol]
  }
  
  # Check the format of X
  if(!(class(X) %in% c("matrix"))){X <- as.matrix(X)}
  #   stop('X must be a data frame or a matrix')
  # if(!(class(X) == "matrix")){X <- as.matrix(X)}
  
  # Check the format of Order
  if(!(Order %in% c("NP","pvals","NC","NCw")))
    stop('Order argument can only be "NP", "pvals", "NC", or "NCw" ')
  
  # Check the format of normalize
  if(!(normalize %in% c("counts","prop","pres")))
    stop('normalize argument can only be "counts", "prop", or "pres" ')
  
  # Check the format of center
  if(class(center) != "logical")
    stop('center argument must be a logical value')
  
  # Check the format of quant
  if(!is.vector(quant)) stop('quant argument must be a vector')
  
  # Check the format of distr
  if(!(distr %in% c("sn","norm","t","cauchy")))
    stop('normalize argument can only be "sn", "norm", "t", or "cauchy" ')
  
  # Check the format of alpha
  if(!is.numeric(alpha)) stop('alpha argument must be a numerical value')
  
  # Check if pvals_sim object is input correctly
  if(class(pvals_sim) != "NULL" & length(pvals_sim$pvals) == 0)
    stop('pvals_sim object must be a result from simultaneous PERFect with taxa abundance ordering')
  
  #Order columns by importance
  if(is.null(Order.user)){
    if(Order == "NP") {Order.vec <- NP_Order(X)}
    if(Order == "pvals") {Order.vec <- pvals_Order(X, pvals_sim)}
    if(Order == "NC"){Order.vec <- NC_Order(X)}
    if(Order == "NCw"){Order.vec <- NCw_Order(X)}
  } else {
    Order.vec <- Order.user #user-specified ordering of columns of X
  }
  
  #properly order columns of X
  X <- X[,Order.vec]
  
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
  
  if(algorithm == "fast"){
    
    #index for iteration
    n <- sum_n(dim(X)[2])$idx-1
    
    #initiate a vector to store p-values
    pvals <- rep(NA,p-1)
    
    #calculate DFL values and FL values
    Order_Ind <- rep(1:length(Order.vec)) #convert to numeric indicator values
    DFL <- DiffFiltLoss(X = X, Order_Ind, Plot = TRUE, Taxa_Names = Order.vec)
    FL <- FiltLoss(X = X, Order.user = Order.vec, type =  "Ind", Plot = TRUE)$FL
    
    #name p-values
    names(pvals) <- names(DFL$DFL)
    
    #For each taxon j, create a distribution of its DFL's by permuting the labels
    if(is.null(dfl_distr)){dfl_distr <- sampl_distr(X = X, k=k)}
    
    #convert to log differences
    lfl <- lapply(dfl_distr, function(x) log(x[!x==0]))
    
    #initiate a vector for OTUs to be double checked
    check = c()
    
    # Calculate the number of cores
    no_cores <- parallel::detectCores()-1
    
    # Initiate cluster, start parrallel processing
    cl <- parallel::makeCluster(no_cores)
    
    if(distr == "norm"){
      # load packages for each core in order to use function qmedist and sn distribution
      parallel::clusterEvalQ(cl,{
        library(fitdistrplus)
      })
      
      # extract the lfl of checked OTUs
      lfl_main = lfl[n]
      
      #check quantile for normal distribution fit
      if(length(quant) > 2){quant <- quant[(length(quant) - 1):length(quant)]
      print("Warning: more than 2 quantile values are given. \nLargest 2 quantiles are used.")}
      if(length(quant) < 2){stop("At least two quantile values must be specified.")}
      
      #load variables for each core
      parallel::clusterExport(cl,c("distr","quant"),envir=environment())
      
      #find the first significant OTU among OTU of index n
      #fit <- parallel::parLapply(cl, lfl_main, function(x) fitdist(x, distr, method="qme", probs=quant))
      fit <- parallel::parLapply(cl, lfl_main, function(x) qmedist(x, distr, probs=quant))
      est <- lapply(fit, function(x) x$estimate)
      for(i in 1:length(n)){
        pvals[n[i]] <- pnorm(q=log(DFL$DFL[n[i]]),mean = est[[i]][1], sd = est[[i]][2],
                             lower.tail = FALSE, log.p = FALSE)
        #calculate a few p-values before that first significant OTU
        if(pvals[n[i]] < alpha){
          stop_idx = n[i-1]
          temp = c((n[i-1]+1):(n[i]-1),n[i]+1, n[i]+2)
          lfl_temp = lfl[temp]
          #fit <- parallel::parLapply(cl,lfl_temp, function(x) fitdist(x, distr, method = "qme", probs=quant))
          fit <- parallel::parLapply(cl,lfl_temp, function(x) qmedist(x, distr, probs=quant))
          est <- lapply(fit, function(x) x$estimate)
          for(j in 1:length(temp)){
            pvals[temp[j]] <- pnorm(q=log(DFL$DFL[temp[j]]), mean = est[[j]][1], sd = est[[j]][2],
                                    lower.tail = FALSE, log.p = FALSE)
            #if(pvals[temp[j]] < 0.1) {break}
          }
          break}
      }
      # find the potential cut off index
      potential <- which(pvals == min(pvals, na.rm = T))
      
      # check the 10% most recent interval for unsual DFL values
      for(i in round((match(stop_idx,n)*9/10)):match(stop_idx,n)){
        check = c(check, names(which(DFL$DFL[n[i]:n[(i+1)]] > max(DFL$DFL[n[i+1]],DFL$DFL[n[i]]) &
                                       FL[n[i]:n[(i+1)]] > max(FL[n[i+1]],FL[n[i]])
        )
        )
        )
      }
      #Identify OTU with DFL higher than the potential cut off taxon 
      check = c(check, names(which(DFL$DFL[1:(potential-1)] > DFL$DFL[potential] &
                                     FL[1:(potential-1)] > FL[potential])))
      
      if (length(check) !=0){
        check_idx = which(names(pvals)%in% check)
        # Calculate their pvalues
        lfl_check = lfl[check_idx]
        #fit <- parallel::parLapply(cl,lfl_check, function(x) fitdist(x, distr,method = "qme", probs=quant))
        fit <- parallel::parLapply(cl,lfl_check, function(x) qmedist(x, distr, probs=quant))
        est <- lapply(fit, function(x) x$estimate)
        for(j in 1:length(check_idx)){
          pvals[check_idx[j]] <- pnorm(q=log(DFL$DFL[check_idx[j]]),mean = est[[j]][1], sd = est[[j]][2],
                                       lower.tail = FALSE, log.p = FALSE)
        }
      }
    }
    
    #fit the distribution
    if(distr == "sn"){
      
      # load packages for each core in order to use function qmedist and sn distribution
      parallel::clusterEvalQ(cl,{
        library(fitdistrplus)
        library(sn)
      })
      
      # starting values to fit the skew normal distribution
      lp <- list(xi = mean(lfl[[1]]), omega = sd(lfl[[1]]), alpha = 1.5)
      
      #get the sampled values from OTU of index n
      lfl_main = lfl[n]
      
      #load variables for each core
      parallel::clusterExport(cl,c("distr","quant","lp"),envir=environment())
      
      #find the first significant OTU among OTU of index n
      #fit <- parallel::parLapply(cl, lfl_main, function(x) fitdist(x, distr,method = "qme", probs=quant, start=lp))
      suppressWarnings(fit <- parallel::parLapply(cl, lfl_main, function(x) qmedist(x, distr, probs=quant, start=lp)))
      est <- lapply(fit, function(x) x$estimate)
      for(i in 1:length(n)){
        pvals[n[i]] <- 1 - psn(x=log(DFL$DFL[n[i]]),xi = est[[i]][1], omega = est[[i]][2], alpha = est[[i]][3])
        #calculate a few p-values before that first significant OTU
        if(pvals[n[i]] < alpha){
          stop_idx = n[i-1]
          temp = c((n[i-1]+1):(n[i]-1),n[i]+1, n[i]+2)
          lfl_temp = lfl[temp]
          #fit <- parallel::parLapply(cl,lfl_temp, function(x) fitdist(x, distr, method = "qme", probs=quant, start=lp))
          suppressWarnings(fit <- parallel::parLapply(cl,lfl_temp, function(x) qmedist(x, distr, probs=quant, start=lp)))
          est <- lapply(fit, function(x) x$estimate)
          for(j in 1:length(temp)){
            pvals[temp[j]] <- 1 - psn(x=log(DFL$DFL[temp[j]]), xi = est[[j]][1], omega = est[[j]][2], alpha = est[[j]][3])
            #if(pvals[temp[j]] < alpha) {break}
          }
          break}
      }
      # find the potential cut off index
      potential <- which(pvals == min(pvals, na.rm = T))
      
      # check the 10% most recent interval for unsual DFL values
      for(i in round((match(stop_idx,n)*9/10)):match(stop_idx,n)){
        check = c(check, names(which(DFL$DFL[n[i]:n[(i+1)]] > max(DFL$DFL[n[i+1]],DFL$DFL[n[i]]) &
                                       FL[n[i]:n[(i+1)]] > max(FL[n[i+1]],FL[n[i]])
        )
        )
        )
      }
      #Identify OTU with DFL higher than the potential cut off taxon 
      check = c(check, names(which(DFL$DFL[1:(potential-1)] > DFL$DFL[potential] &
                                     FL[1:(potential-1)] > FL[potential])
      )
      )
      if (length(check) != 0){
        check_idx = which(names(pvals)%in% check)
        # Calculate their pvalues
        lfl_check = lfl[check_idx]
        #fit <- parallel::parLapply(cl,lfl_check, function(x) fitdist(x, distr, method = "qme", probs=quant, start=lp))
        suppressWarnings(fit <- parallel::parLapply(cl,lfl_check, function(x) qmedist(x, distr, probs=quant, start=lp)))
        est <- lapply(fit, function(x) x$estimate)
        for(j in 1:length(check_idx)){
          pvals[check_idx[j]] <- 1 - psn(x=log(DFL$DFL[check_idx[j]]),
                                         xi = est[[j]][1], omega = est[[j]][2], alpha = est[[j]][3])
        }
      }
    }
    
    # End the parallel processing
    parallel::stopCluster(cl)
    
    #smooth p-values
    non_na_ind <- which(!is.na(pvals))
    pvals_avg <- pvals
    pvals_avg[non_na_ind] <- zoo::rollapply(pvals[non_na_ind], width =lag, mean, align=direction, fill=NA, na.rm = T)
    #replace na's with original values
    pvals_avg[non_na_ind][is.na(pvals_avg[non_na_ind])] <- pvals[non_na_ind][is.na(pvals_avg[non_na_ind])]
    names(pvals_avg) <- names(pvals)[(length(pvals)-length(pvals_avg)+1):length(pvals)]
    
    #select the first significant taxon
    Ind = which(pvals_avg <=alpha & pvals_avg > 0)[1]
    
    if (length(Ind !=0)) {Ind <- min(Ind)}
    else{Ind <- dim(X)[2]-1
    warning("no taxa are significant at a specified alpha level")}
    
    #if jth DFL is significant, then throw away all taxa 1:j
    filtX <- X.orig[,-(1:Ind)]
    
    return(list(filtX = filtX, info = info, pvals = pvals_avg))
  }  else if(algorithm == "full"){
    pvals <- rep(0,p-1)
    hist_list <- lapply(1:(p-1),function(x) NULL)
    
    #calculate DFL values
    Order_Ind <- rep(1:length(Order.vec)) #convert to numeric indicator values
    DFL <- DiffFiltLoss(X = X, Order_Ind, Plot = TRUE, Taxa_Names = Order.vec)
    
    #name p-values
    names(pvals) <- names(DFL$DFL)
    
    #For each taxon j, create a distribution of its DFL's by permuting the labels
    if(is.null(dfl_distr)){dfl_distr <- sampl_distr(X = X, k=k)}
    
    #convert to log differences
    lfl <- lapply(dfl_distr, function(x) log(x[!x==0]))
    
    #fit the distribution
    if(distr == "norm"){
      if(length(quant) > 2){quant <- quant[(length(quant) - 1):length(quant)]
      print("Warning: more than 2 quantile values are given. \nLargest 2 quantiles are used.")}
      if(length(quant) < 2){stop("At least two quantile values must be specified.")}
      #fit <- lapply(lfl, function(x) fitdist(x, distr, method = "qme", probs=quant))
      fit <- lapply(lfl, function(x) qmedist(x, distr, probs=quant))
      est <- lapply(fit, function(x) x$estimate)
      for(i in 1:(p-1)){
        pvals[i] <- pnorm(q=log(DFL$DFL[i]), mean = est[[i]][1], sd = est[[i]][2],
                          lower.tail = FALSE, log.p = FALSE)}
    }
    
    if(distr == "sn"){
      #lp <- lapply(lfl, function(x) list(xi = mean(x), omega = sd(x), alpha = 1.5))
      lp <- list(xi = mean(lfl[[1]]), omega = sd(lfl[[1]]), alpha = 1.5)
      #fit <- lapply(lfl, function(x) fitdist(x, distr, method = "qme", probs=quant, start=lp))
      suppressWarnings(fit <- lapply(lfl, function(x) qmedist(x, distr, probs=quant, start=lp)))
      est <- lapply(fit, function(x) x$estimate)
      for(i in 1:(p-1)){
        pvals[i] <- 1 - psn(x=log(DFL$DFL[i]),
                            xi = est[[i]][1], omega = est[[i]][2], alpha = est[[i]][3])}
    }
    
    if(hist == TRUE){
      lfl <- lapply(lfl, function(x) data.frame(x))
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
          ggtitle("") + xlab("log differences in filtering loss") + ylab("Density")
        
        #estimate using normal
        if(distr == "norm"){
          #add density line to the plot
          hist <- hist + stat_function(fun = dnorm,
                                       args = list(mean = est[[i]][1], sd = est[[i]][2]), colour=linecol)
          hist_list[[i]] <- hist
        }
        
        #estimate using skew normal
        if(distr == "sn"){
          hist <- hist + stat_function(fun = dsn,
                                       args = list(xi = est[[i]][1], omega = est[[i]][2], alpha = est[[i]][3]), colour=linecol)
          hist_list[[i]] <- hist
        }
      }
    }#end if hist == TRUE
    #smooth p-values
    pvals_avg <-zoo::rollapply(pvals, width =lag, mean, align=direction, fill=NA, na.rm = T)
    #replace na's with original values
    pvals_avg[is.na(pvals_avg)] <- pvals[is.na(pvals_avg)]
    names(pvals_avg) <- names(pvals)
    
    #select taxa that are kept in the data set at significance level alpha
    Ind <- which(pvals_avg <=alpha)
    if (length(Ind !=0)) {Ind <- min(Ind)}
    else{Ind <- dim(X)[2]-1
    warning("no taxa are significant at a specified alpha level")}
    #if jth DFL is significant, then throw away all taxa 1:j
    filtX <- X.orig[,-(1:Ind)]
    return(list(filtX = filtX,info = info, fit = fit, hist = hist_list,
                est =est, dfl_distr=dfl_distr, pvals = pvals_avg ))
  }# end if(algorithm == "full")
  
}
