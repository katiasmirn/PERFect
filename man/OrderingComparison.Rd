\name{OrderingComparison}
\alias{OrderingComparison}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
OrderingComparison(X, Order, res_perm, nbins = 30, col = "red", fill = "green", alpha = 0.05, distr = "sn", quant = c(0.1, 0.25, 0.5), hist_fill = 0.2, linecol = "blue")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Order}{
%%     ~~Describe \code{Order} here~~
}
  \item{res_perm}{
%%     ~~Describe \code{res_perm} here~~
}
  \item{nbins}{
%%     ~~Describe \code{nbins} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{fill}{
%%     ~~Describe \code{fill} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{distr}{
%%     ~~Describe \code{distr} here~~
}
  \item{quant}{
%%     ~~Describe \code{quant} here~~
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
function (X, Order, res_perm, nbins = 30, col = "red", fill = "green", 
    alpha = 0.05, distr = "sn", quant = c(0.1, 0.25, 0.5), hist_fill = 0.2, 
    linecol = "blue") 
{
    res_sim <- PERFect_sim(X = X, Order = Order, nbins = 30, 
        col = "red", fill = "green", alpha = 0.05, distr = "sn", 
        quant = c(0.1, 0.25, 0.5), hist_fill = 0.2, linecol = "blue")
    res_perm <- PERFect_perm_reorder(X = X, Order_alt = Order, 
        res_perm = res_perm, alpha = alpha, distr = distr)
    PERFect_sim_taxa <- names(res_sim$filtX)
    PERFect_perm_taxa <- names(res_perm$filtX)
    Loss_PERFect_sim <- FL_J(X = X, J = PERFect_sim_taxa, leave = TRUE)
    Loss_PERFect_perm <- FL_J(X = X, J = PERFect_perm_taxa, leave = TRUE)
    Res <- cbind(rbind(length(PERFect_sim_taxa), length(PERFect_perm_taxa)), 
        rbind(Loss_PERFect_sim, Loss_PERFect_perm))
    colnames(Res) <- c("N_taxa", "FL")
    rownames(Res) <- c("PERFect_sim", "PERFect_perm")
    return(list(Res = Res, res_sim = res_sim, res_perm = res_perm))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
