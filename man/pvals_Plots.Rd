\name{pvals_Plots}
\alias{pvals_Plots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots of PERFect p-values 
}
\description{
Graphical representation of p-values obtained by running PERFect_sim() or  PERFect_perm() for jth taxon colored by quantile values of individual filtering loss 
}
\usage{
pvals_Plots(PERFect, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PERFect}{Output of PERFect_sim() or PERFect_perm() function}
  \item{quantiles}{Quantile values  for coloring, these are set to 25\%, 50\%, 0.8\% and 0.9\% percentiles of the ndividual filtering  loss values}
  \item{alpha}{alpha level of the test}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  \item{p_pvals}{Plot of p-values}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Ekaterina Smirnova
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also as \code{\link{PERFect_sim}}, ~\code{\link{PERFect_perm}}
}
%\examples{

%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
