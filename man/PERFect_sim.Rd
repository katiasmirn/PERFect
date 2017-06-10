\name{PERFect_sim}
\alias{PERFect_sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PERFect_sim(X, Order, nbins = 30, quant = c(0.25, 0.5), distr = c("norm", "sn"), alpha = 0.05, col = "red", fill = "green", hist_fill = 0.2, linecol = "blue")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Order}{
%%     ~~Describe \code{Order} here~~
}
  \item{nbins}{
%%     ~~Describe \code{nbins} here~~
}
  \item{quant}{
%%     ~~Describe \code{quant} here~~
}
  \item{distr}{
%%     ~~Describe \code{distr} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{fill}{
%%     ~~Describe \code{fill} here~~
}
  \item{hist_fill}{
%%     ~~Describe \code{hist_fill} here~~
}
  \item{linecol}{
%%     ~~Describe \code{linecol} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, Order, nbins = 30, quant = c(0.25, 0.5), distr = c("norm", 
    "sn"), alpha = 0.05, col = "red", fill = "green", hist_fill = 0.2, 
    linecol = "blue") 
{
    pDFL <- NULL
    phist <- NULL
    X <- X[, Order]
    p <- dim(X)[2]
    Order_Ind <- rep(1:length(Order))
    DFL <- DiffFiltLoss(X = X, Order_Ind, Plot = TRUE, Taxa_Names = Order)
    Taxa <- Order[-length(Order)]
    lfl <- data.frame(Taxa, log(DFL$DFL))
    names(lfl) <- c("Taxa", "DFL")
    hist <- ggplot(data = lfl, aes(lfl$DFL)) + geom_histogram(bins = nbins, 
        aes(y = ..density..), col = col, fill = fill, alpha = hist_fill) + 
        theme(panel.background = element_rect(fill = "white"), 
            panel.grid.major = element_line(colour = "grey90"), 
            axis.text.x = element_text(size = 10)) + ggtitle("") + 
        xlab("Filtering Loss") + ylab("Density")
    if (distr == "norm") {
        fit <- qmedist(lfl$DFL, distr, probs = quant)
        est <- fit$estimate
        hist <- hist + stat_function(fun = dnorm, args = list(mean = est[1], 
            sd = est[2]), colour = linecol)
        pvals <- pnorm(q = lfl$DFL, mean = est[1], sd = est[2], 
            lower.tail = FALSE, log.p = FALSE)
    }
    if (distr == "sn") {
        lp <- list(xi = mean(lfl$DFL), omega = sd(lfl$DFL), alpha = 1.5)
        fit <- qmedist(lfl$DFL, distr, probs = quant, start = lp)
        est <- fit$estimate
        hist <- hist + stat_function(fun = dsn, args = list(xi = est[1], 
            omega = est[2], alpha = est[3]), colour = linecol)
        pvals <- 1 - psn(x = lfl$DFL, xi = est[1], omega = est[2], 
            alpha = est[3])
    }
    names(pvals) <- names(DFL$DFL)
    Ind <- which(pvals <= alpha)
    if (length(Ind != 0)) {
        Ind <- min(Ind)
    }
    else {
        Ind <- dim(X)[2] - 1
        warning("no taxa are significant at a specified alpha level")
    }
    filtX <- X[, -(1:Ind)]
    return(list(hist = hist, est = est, pvals = round(pvals, 
        5), filtX = filtX, DFL = DFL$DFL, pDFL = DFL$p + ylab("Difference in Filtering Loss"), 
        fit = fit))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
