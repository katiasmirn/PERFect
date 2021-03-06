#' Permutation PERFect filtering for microbiome data
#'
#' @description This function  filters  the provided OTU table X at a test level alpha given a fitted object
#' perfect_perm obtained by running \code{PERFect_perm()} function. \code{PERFect_perm_reorder()} reavaluates taxa
#' significance p-values for a different taxa ordering.
#'
#' @usage PERFect_perm_reorder(X, Order = "NP", Order.user = NULL, res_perm, normalize = "counts",
#'     center = FALSE, alpha = 0.1, distr = "sn", rollmean = TRUE, direction = "left",
#'     pvals_sim = NULL)
#'
#' @param X OTU table, where taxa are columns and samples are rows of the table.
#' It should be a in data frame format with columns corresponding to taxa names.
#'
#' @param Order Taxa ordering. The default ordering is the number of occurrences (NP) of the taxa in all samples.
#'  Other types of order are p-value ordering, number of connected taxa and weighted number of connected taxa,
#'  denoted as \code{"pvals"}, \code{"NC"}, \code{"NCw"} respectively. More details about taxa ordering are described in Smirnova et al.
#'  User can also specify their preference order with Order.user.
#'
#' @param Order.user User's taxa ordering. This argument takes a character vector of ordered taxa names.
#'
#' @param res_perm Output of \code{PERFect_perm()} function.
#'
#' @param normalize Normalizing taxa count. The default option does not normalize taxa count,
#'  but user can convert the OTU table into a proportion table using the option \code{"prop"}
#'  or convert it into a presence/absence table using \code{"pres"}.
#'
#' @param center Centering OTU table. The default option does not center the OTU table.
#'
#' @param alpha Test level alpha, set to 0.1 by default.
#'
#' @param distr The type of distribution used in \code{PERFect_perm()} function to obtain res_perm object.
#' \describe{
#' \item{\code{"sn"}}{Skew-Normal distribution with 3 parameters: location xi, scale omega^2 and shape alpha}
#' \item{\code{"norm"}}{Normal distribution with 2 parameters: mean and standard deviation sd}
#' \item{\code{"t"}}{Student t-distribution with 2 parameters: n degrees of freedom and noncentrality ncp}
#' \item{\code{"cauchy"}}{Cauchy distribution with 2 parameters: location and scale}
#' }
#'
#' @param rollmean Binary TRUE/FALSE value. If TRUE, rolling average (moving mean) of p-values will be calculated,
#'  with the lag window set to 3 by default.
#'
#' @param direction Character specifying whether the index of the result should be left- or right-aligned
#'  or centered compared to the rolling window of observations, set to "left" by default.
#'
#' @param pvals_sim Object resulting from simultaneous PERFect with taxa abundance ordering,
#'  allowing user to perform Simultaneous PERFect with p-values ordering.
#'  Be aware that the choice of distribution for both methods must be the same.
#'
#' @details This function is designed to save computational time needed to obtain and fit the sampling distribution
#' for each taxon if taxa ordering different from the one used in \code{PERFect_perm()} is used.
#' Note, the distribution and OTU table X should match the distribution used in \code{PERFect_perm()}.
#'
#' @return
#'
#' \item{res_perm}{The perfect_perm object updated according to the alternative taxa ordering.
#' All elements in this list are same as in perfect_perm object given by \code{PERFect()} function.}
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
#'
#' # Proportion data matrix
#' Prop <- mock2$Prop
#'
#' # Counts data matrix
#' Counts <- mock2$Counts
#'
#' #### Uncomment to run algorithm with parallel processing with more than 2 cores
#' # #obtain permutation PERFEct results using NP taxa ordering
#' # system.time(res_perm <- PERFect_perm(X=Prop, k = 1000, algorithm = "fast))
#'
#' # #run PERFEct_sim() function and obtain p-values ordering
#' # res_sim <- PERFect_sim(X=Prop)
#'
#' # #order according to p-values
#' # pvals_sim <- pvals_Order(Counts, res_sim)
#'
#' # #update perfect_perm object according to p-values ordering
#' # res_reorder <- PERFect_perm_reorder(X=Prop,  Order.user = pvals_sim,  res_perm = res_perm)
#'
#' # #permutation perfect colored by FLu values
#' # pvals_Plots(PERFect = res_perm, X = Counts, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.05)
#' @export

PERFect_perm_reorder <- function(X,  Order ="NP",  Order.user = NULL, res_perm,
                                 normalize = "counts", center = FALSE, alpha = 0.10, distr = "sn",
                                 rollmean = TRUE, direction ="left", pvals_sim = NULL){
  # Check the format of X
  if (!(class(X) %in% c("matrix"))) {
    X <- as.matrix(X)
  }
  #   stop('X must be a data frame or a matrix')
  # if(!(class(X) == "matrix")){X <- as.matrix(X)}

  # Check the format of Order
  if (!(Order %in% c("NP", "pvals", "NC", "NCw")))
    stop('Order argument can only be "NP", "pvals", "NC", or "NCw" ')

  # Check the format of res_perm
  if (class(res_perm) != "NULL" & length(res_perm$pvals) == 0)
    stop('res_perm argument must be the output from the function PERFect_perm()')

  # Check the format of normalize
  if (!(normalize %in% c("counts", "prop", "pres")))
    stop('normalize argument can only be "counts", "prop", or "pres" ')

  # Check the format of center
  if (class(center) != "logical")
    stop('center argument must be a logical value')

  # Check the format of distr
  if (!(distr %in% c("sn", "norm", "t", "cauchy")))
    stop('normalize argument can only be "sn", "norm", "t", or "cauchy" ')

  # Check the format of alpha
  if (!is.numeric(alpha))
    stop('alpha argument must be a numerical value')

  #Order columns by importance
  if (Order == "NP") {
    Order.vec <- NP_Order(X)
  }
  if (Order == "pvals") {
    Order.vec <- pvals_Order(X, pvals_sim)
  }
  if (Order == "NC") {
    Order.vec <- NC_Order(X)
  }
  if (Order == "NCw") {
    Order.vec <- NCw_Order(X)
  }
  else if (!is.null(Order.user)) {
    Order.vec = Order.user
  } #user-specified ordering of columns of X

  X <- X[, Order.vec]#properly order columns of X
  #save non-centered, non-normalized X
  X.orig <- X

  #remove all-zero OTU columns
  nzero.otu <- apply(X, 2, Matrix::nnzero) != 0
  X <- X[, nzero.otu]
  p <- dim(X)[2]
  Order.vec <- Order.vec[nzero.otu]

  #normalize the data
  if (normalize == "prop") {
    X <- X / apply(X, 1, sum)
  }
  else if (normalize == "pres") {
    X[X != 0] <- 1
  }

  #center if true
  if (center) {
    X <- apply(X, 2, function(x) {
      x - mean(x)
    })
  }

  Order_Ind <-
    rep(seq_len(length(Order.vec)))#convert to numeric indicator values
  DFL <-
    DiffFiltLoss(X = X,
                 Order_Ind,
                 Plot = TRUE,
                 Taxa_Names = Order.vec)
  #re-evaluate p-values
  pvals <- rep(0, length(DFL$DFL))
  names(pvals) <- names(DFL$DFL)
  for (i in seq_len(DFL$DFL)) {
    if (distr == "sn") {
      pvals[i] <- 1 - psn(
        x = log(DFL$DFL[i]),
        xi = res_perm$est[[i]][1],
        omega = res_perm$est[[i]][2],
        alpha = res_perm$est[[i]][3]
      )
    }
    if (distr == "norm") {
      #calculate p-values
      pvals[i] <-
        pnorm(
          q = log(DFL$DFL[i]),
          mean = res_perm$est[[i]][1],
          sd = res_perm$est[[i]][2],
          lower.tail = FALSE,
          log.p = FALSE
        )
    }

  }

  #re-calculate filtered X
  #smooth p-values
  if (rollmean){
  pvals_avg <-
    zoo::rollmean(pvals,
                  k = 3,
                  align = direction,
                  fill = NA)
  } else {
    pvals_avg <- pvals
  }
  #replace na's with original values
  pvals_avg[is.na(pvals_avg)] <- pvals[is.na(pvals_avg)]
  #select taxa that are kept in the data set at significance level alpha
  Ind <- which(pvals_avg <= alpha)
  if (length(Ind != 0)) {
    Ind <- min(Ind)
  }
  else{
    Ind <- dim(X)[2] - 1
    warning("no taxa are significant at a specified alpha level")
  }
  #if jth DFL is significant, then throw away all taxa 1:j
  res_perm$filtX <- X.orig[, -seq_len(Ind)]
  res_perm$pvals <- pvals_avg #end if !is.null(res_perm)

  return(res_perm)
}
