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
#' @usage PERFect_perm(X, Order = "NP", Order.user = NULL, normalize = "counts",
#'     center = FALSE, quant = c(0.1, 0.25, 0.5), distr = "sn",
#'     alpha = 0.1, lag = 3, direction = "left", pvals_sim = NULL,
#'     k = 1000, dfl_distr = NULL, nbins = 30, hist = FALSE, col = "red",
#'     fill = "green", hist_fill = 0.2, linecol = "blue")
#'
#' @param X OTU table, where taxa are columns and samples are rows of the table.
#' It should be a in data frame format with columns corresponding to taxa names.
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
#' \item{\code{"t"}}{Student t-distribution with 2 parameters: n degrees of freedom and noncentrality ncp}
#' \item{\code{"cauchy"}}{Cauchy distribution with 2 parameters: location and scale}
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
PERFect_perm <- function(X,  Order = "NP",   Order.user = NULL,
                         normalize = "counts", center = FALSE,
                         quant = c(0.10, 0.25, 0.5),  distr ="sn",
                         alpha = 0.10, lag = 3, direction ="left",
                         pvals_sim = NULL,
                         k=1000, dfl_distr = NULL,
                         nbins =30,  hist = FALSE,
                         col = "red", fill = "green", hist_fill = 0.2, linecol = "blue"){
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

  #convert X to matrix
  if(!(class(X) == "matrix")){X <- as.matrix(X)}

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

  pvals <- rep(0,p-1)
  hist_list <- lapply(1:(p-1),function(x) NULL)

  #calculate loss
  #calculate DFL values
  Order_Ind <- rep(1:length(Order.vec))#convert to numeric indicator values
  DFL <- DiffFiltLoss(X = X, Order_Ind, Plot = TRUE, Taxa_Names = Order.vec)
  #name p-values
  names(pvals) <- names(DFL$DFL)
  #For each taxon j, create a distribution of its DFL's by permuting the labels
  if(is.null(dfl_distr)){dfl_distr <- sampl_distr(X = X, k=k)}
  #convert to log differences
  lfl <- lapply(dfl_distr, function(x) log(x[!x==0]))
  #fit the distribution
  if(distr == "norm"){
    fit <- lapply(lfl, function(x) qmedist(x, distr, probs=quant))
    est <- lapply(fit, function(x) x$estimate)
    for(i in 1:(p-1)){
        pvals[i] <- pnorm(q=log(DFL$DFL[i]), mean = est[[i]][1], sd = est[[i]][2],
                      lower.tail = FALSE, log.p = FALSE)}
  }

  if(distr == "sn"){
    #lp <- lapply(lfl, function(x) list(xi = mean(x), omega = sd(x), alpha = 1.5))
    lp <- list(xi = mean(lfl[[1]]), omega = sd(lfl[[1]]), alpha = 1.5)
    fit <- lapply(lfl, function(x) qmedist(x, distr, probs=quant, start=lp))
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
  pvals_avg <- rollmean(pvals, k=lag, align=direction,  fill=NA )

  #replace na's with original values
  pvals_avg[is.na(pvals_avg)] <- pvals[is.na(pvals_avg)]

  #select taxa that are kept in the data set at significance level alpha
  Ind <- which(pvals_avg <=alpha)
  if (length(Ind !=0)) {Ind <- min(Ind)}
  else{Ind <- dim(X)[2]-1
       warning("no taxa are significant at a specified alpha level")}
  #if jth DFL is significant, then throw away all taxa 1:j
  filtX <- X.orig[,-(1:Ind)]
  return(list(filtX = filtX, pvals = pvals_avg, fit = fit, hist = hist_list,
              est =est,   dfl_distr=dfl_distr ))
}

