\name{TraditR2}
\alias{TraditR2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TraditR2(X, Ab_min = 0.001, prop = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Ab_min}{
%%     ~~Describe \code{Ab_min} here~~
}
  \item{prop}{
%%     ~~Describe \code{prop} here~~
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
function (X, Ab_min = 0.001, prop = TRUE) 
{
    if (prop == FALSE) {
        X <- sweep(Counts, MARGIN = 1, apply(Counts, 1, sum), 
            "/")
    }
    n <- dim(X)[1]
    Abund <- apply(X, 2, max)
    Nsamp <- apply(X, 2, nnzero)
    Psamp <- Nsamp/n
    if (is.null(names(Abund))) {
        names(Abund) <- 1:length(Abund)
    }
    Taxa <- names(Abund[Abund > Ab_min])
    Abund <- Abund[Taxa]
    Nsamp <- Nsamp[Taxa]
    Psamp <- Psamp[Taxa]
    C1 <- names(Abund[Abund > 0.01 & Nsamp >= 1])
    C2 <- names(Abund[Abund > 0.001 & Psamp >= 0.02])
    C3 <- names(Abund[Psamp >= 0.05])
    Taxa <- unique(c(C1, C2, C3))
    filtX <- X[, Taxa]
    return(filtX)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
