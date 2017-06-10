\name{FL_J}
\alias{FL_J}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FL_J(X, J, leave = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{J}{
%%     ~~Describe \code{J} here~~
}
  \item{leave}{
%%     ~~Describe \code{leave} here~~
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
function (X, J, leave = TRUE) 
{
    Ind <- which(names(X) \%in\% J)
    if (leave == TRUE) {
        X_R <- X[, Ind]
    }
    else {
        X_R <- X[, -Ind]
    }
    Netw <- t(as.matrix(X)) \%*\% as.matrix(X)
    Netw_R <- t(as.matrix(X_R)) \%*\% as.matrix(X_R)
    FL <- 1 - tr(t(Netw_R) \%*\% Netw_R)/tr(t(Netw) \%*\% Netw)
    return(FL)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
