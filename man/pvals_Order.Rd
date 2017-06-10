\name{pvals_Order}
\alias{pvals_Order}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pvals_Order(Counts, res_sim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Counts}{
%%     ~~Describe \code{Counts} here~~
}
  \item{res_sim}{
%%     ~~Describe \code{res_sim} here~~
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
function (Counts, res_sim) 
{
    vec <- names(res_sim$pvals)
    Order <- c(setdiff(names(Counts), vec), vec)
    Order_pvals <- Order[c(1, sort.int(res_sim$pvals, index.return = TRUE, 
        decreasing = TRUE)$ix + 1)]
    return(Order_pvals)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
