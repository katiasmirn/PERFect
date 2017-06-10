\name{PERFect_perm_reorder}
\alias{PERFect_perm_reorder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PERFect_perm_reorder(X, Order_alt, res_perm, alpha = 0.05, distr = "sn")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Order_alt}{
%%     ~~Describe \code{Order_alt} here~~
}
  \item{res_perm}{
%%     ~~Describe \code{res_perm} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{distr}{
%%     ~~Describe \code{distr} here~~
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
function (X, Order_alt, res_perm, alpha = 0.05, distr = "sn") 
{
    X2 <- X
    X2 <- X2[, Order_alt]
    Order_Ind <- rep(1:length(Order))
    DFL <- DiffFiltLoss(X = X2, Order_Ind, Plot = TRUE, Taxa_Names = Order_alt)
    pvals <- rep(0, length(DFL$DFL))
    names(pvals) <- names(DFL$DFL)
    for (i in 1:length(DFL$DFL)) {
        if (distr == "sn") {
            pvals[i] <- 1 - psn(x = log(DFL$DFL[i]), xi = res_perm$est[[i]][1], 
                omega = res_perm$est[[i]][2], alpha = res_perm$est[[i]][3])
        }
        if (distr == "norm") {
            pvals[i] <- pnorm(q = log(DFL$DFL[i]), mean = res_perm$est[[i]][1], 
                sd = res_perm$est[[i]][2], lower.tail = FALSE, 
                log.p = FALSE)
        }
    }
    Ind <- which(pvals <= alpha)
    if (length(Ind != 0)) {
        Ind <- min(Ind)
    }
    else {
        Ind <- dim(X2)[2] - 1
        warning("no taxa are significant at a specified alpha level")
    }
    res_perm$filtX <- X2[, -(1:Ind)]
    res_perm$pvals <- pvals
    return(res_perm)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
