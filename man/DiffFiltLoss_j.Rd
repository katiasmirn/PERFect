\name{DiffFiltLoss_j}
\alias{DiffFiltLoss_j}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DiffFiltLoss_j(Perm_Order, Netw, j)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Perm_Order}{
%%     ~~Describe \code{Perm_Order} here~~
}
  \item{Netw}{
%%     ~~Describe \code{Netw} here~~
}
  \item{j}{
%%     ~~Describe \code{j} here~~
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
function (Perm_Order, Netw, j) 
{
    J_jp1 <- Perm_Order[1:(j + 1)]
    DFL <- (2 * t(Netw[J_jp1[max(NROW(J_jp1))], -J_jp1]) \%*\% 
        Netw[J_jp1[max(NROW(J_jp1))], -J_jp1] + Netw[J_jp1[max(NROW(J_jp1))], 
        J_jp1[max(NROW(J_jp1))]]^2)
    return(DFL)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
