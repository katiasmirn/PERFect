##############################
#define skew normal distribution to be used in fitristplus
#############################
# require(sn);
# require(fitdistrplus);
# dsmnorm<<-function(x, mean = 0, sd = 1, xi = 1.5, log = FALSE){return(dsnorm(x, mean = 0, sd = 1, xi = 1.5, log = FALSE))}
# psmnorm<<-function(q, mean = 0, sd = 1, xi = 1.5){return(psnorm(q, mean = 0, sd = 1, xi = 1.5))}
# qsmnorm<<-function(p, mean = 0, sd = 1, xi = 1.5){return(qsnorm(p, mean = 0, sd = 1, xi = 1.5))}
# rsmnorm<<-function(n, mean = 0, sd = 1, xi = 1.5){return(rsnorm(n, mean = 0, sd = 1, xi = 1.5))}

#available distributions:
#normal - "norm"
#skew-normal - "sn"
#non-central t - "t"
#non-central beta - "cauchy"
######################################
#Fit quantiles to the log data and get p-values
######################################
#' Simulation PERFect filtering for microbiome data
#'
#' @description Simultaneous filtering of the provided OTU table X at a test level alpha. One distribution is fit to taxa simultaneously.
#'
#' @usage PERFect_sim(X,infocol = NULL, Order = "NP", Order.user = NULL, normalize = "counts",
#'          center = FALSE, quant = c(0.1, 0.25, 0.5), distr = "sn",
#'          alpha = 0.1, lag = 3, direction = "left", pvals_sim = NULL,
#'          nbins = 30, col = "red", fill = "green", hist_fill = 0.2,
#'          linecol = "blue")
#'
#' @param X OTU table, where taxa are columns and samples are rows of the table.
#' It should be a in data frame format with columns corresponding to taxa names.
#' It could contains columns of metadata.
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
#' @param Order.user User's taxa ordering. This argument takes a character vector of ordered taxa names.
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
#' Function \code{PERFect_sim()}  filters  the provided OTU table X and outputs a filtered table
#' that contains signal taxa. \code{PERFect_sim()} calculates differences in filtering loss DFL
#' for each taxon according to the given taxa order. By default, the function fits Skew-Normal distribution
#' to the log-differences in filtering loss but Normal, t, or Cauchy distributions can be also used.
#' This is implementation of Algorithm 1 described in Smirnova et al.
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
#' \item{pDFL}{Plot of differences in filtering loss values}
#'
#' @references Azzalini, A. (2005). The skew-normal distribution and related multivariate families. Scandinavian Journal of Statistics, 32(2), 159-188.
#'
#' @references Smirnova, E., Huzurbazar, H., Jafari, F. ``PERFect: permutationfiltration of microbiome data", to be submitted.
#'
#' @author Ekaterina Smirnova
#'
#' @seealso \code{\link{PERFect_perm}}
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
#'
#' @import phyloseq
#' @import ggplot2
#' @importFrom Matrix nnzero
#' @importFrom sn qsn psn dsn
#' @importFrom fitdistrplus qmedist
#' @importFrom psych tr
#' @importFrom zoo rollmean rollapply
#' @importFrom stats dcauchy dnorm dt pcauchy pnorm pt quantile sd
#' @export
PERFect_sim <- function(X, infocol= NULL,  Order = "NP",   Order.user = NULL,
                        normalize = "counts", center = FALSE,
                        quant = c(0.10, 0.25, 0.5),  distr ="sn",
                        alpha = 0.10, lag = 3, direction ="left",
                        pvals_sim = NULL,
                        nbins =30,
                        col = "red", fill = "green", hist_fill = 0.2, linecol = "blue"){

  pDFL <- NULL
  phist <- NULL
  info <- NULL
  #infocol = index vector of other info
  if(!is.null(infocol)){
    info <- X[,infocol]
    X <- X[,-infocol]
  }

  # Check the format of X
  if(!(class(X) %in% c("matrix"))){X <- as.matrix(X)}
  #  stop('X must be a data frame or a matrix')
  #if(!(class(X) == "matrix")){X <- as.matrix(X)}

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
  X <- X[,Order.vec]#properly order columns of X

  #remove all-zero OTU columns
  nzero.otu <- apply(X, 2, Matrix::nnzero) != 0
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
  Order_Ind <- rep(seq_len(length(Order.vec)))#convert to numeric indicator values
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
    fit <- fitdistrplus::qmedist(lfl$DFL, distr, probs=quant)
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
    fit <- fitdistrplus::qmedist(lfl$DFL, distr, probs=quant, start=list(df=2, ncp = mean(lfl$DFL)))
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
    fit <- fitdistrplus::qmedist(lfl$DFL, distr, probs=quant)
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
    suppressWarnings(fit <- fitdistrplus::qmedist(lfl$DFL, distr, probs=quant, start=lp))
    #fit <- fitdist(lfl$DFL, distr, method = "qme", probs=quant, start=lp)
    est <- fit$estimate
    hist <- hist + stat_function(fun = dsn, args = list(xi = est[1], omega = est[2], alpha = est[3]), colour=linecol)
    #calculate p-values
    pvals <- 1- psn(x=lfl$DFL, xi = est[1], omega = est[2], alpha = est[3])
  }

  #select taxa that are kept in the data set at significance level alpha
  names(pvals) <- names(DFL$DFL)

  #smooth p-values
  pvals_avg <- zoo::rollmean(pvals, k=lag, align=direction,  fill=NA )
  #replace na's with original values
  pvals_avg[is.na(pvals_avg)] <- pvals[is.na(pvals_avg)]

  Ind <- which(pvals_avg <=alpha)
  if (length(Ind !=0)) {Ind <- min(Ind)}
  else{Ind <- dim(X)[2]-1
  warning("no taxa are significant at a specified alpha level")}
  #if jth DFL is significant, then throw away all taxa 1:j
  filtX <- X.orig[,-seq_len(Ind)]

  return(list(filtX = filtX, info = info, pvals = round(pvals_avg,5), DFL = DFL$DFL, fit=fit, hist = hist, est = est,
              pDFL = DFL$p + ylab("Difference in Filtering Loss")))
}

