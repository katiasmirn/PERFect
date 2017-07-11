\name{pvals_Order}
\alias{pvals_Order}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Taxa importance ordering by PERFect p-values
}
\description{
This function orders taxa by increasing significance of simultaneous PERFect p-values
}
\usage{
pvals_Order(Counts, res_sim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Counts}{
OTU COUNTS table, where taxa are columns and samples are rows of the table. Should be in data frame format with columns corresponding to taxa names. Must be a data frame, where taxa are column names. 
}
  \item{res_sim}{
Output of PERFect_sim() function
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
 \item{Order_pvals}{Taxa names in increasing order of p-values significance}
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
See Also as \code{\link{PERFect_sim}}, \code{\link{PERFect_perm}}
}
\examples{
require(Matrix)
require(ggplot2)
require(sn)
require(fitdistrplus)
require(psych)

data(mock2)
Prop <- mock2$Prop
Counts <- mock2$Counts

#arrange counts in order of increasing number of samples taxa are present in 
NP <- NP_Order(Counts)
Counts <- Counts[,NP]

res_sim <- PERFect_sim(X=Prop,  Order=NP,  nbins = 30, col = "red", fill = "green", alpha = 0.05,distr = "sn", quant = c(0.1, 0.25, 0.5), hist_fill =0.2, linecol = "blue")

#order according to p-values
pvals_sim <- pvals_Order(Counts, res_sim)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
