\name{NCw_Order}
\alias{NCw_Order}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
NCw_Order(Counts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Counts}{
%%     ~~Describe \code{Counts} here~~
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
function (Counts) 
{
    Netw <- t(as.matrix(Counts)) \%*\% as.matrix(Counts)
    Netw2 <- Netw
    diag(Netw2) <- 0
    NC_val <- sort(apply(Netw2, 2, nnzero), decreasing = FALSE)
    NC <- names(NC_val)
    n <- dim(Counts)[1]
    n_Tilde <- apply(Counts[, NC], 2, nnzero)
    NCW_val <- sort((n_Tilde * NC_val)/n, decreasing = FALSE)
    NCW <- names(NCW_val)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
