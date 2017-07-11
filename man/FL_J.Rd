\name{FL_J}
\alias{FL_J}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filtering Loss for a set of filtered taxa J
}
\description{
This function calculates filtering loss for the candidate set J of filtered out taxa 
}
\usage{
FL_J(X, J)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
An n x p OTU table, where p taxa are columns and n samples are rows of the table. Should be in data frame format with columns corresponding to taxa names.
}
  \item{J}{
The set J of filtered out taxa, which must be subset of column names of X.
}

}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  \item{FL}{Filtering loss value}
}
\references{
Smirnova, E., Huzurbazar, H., Jafari, F. ``PERFect: permutation  filtration of microbiome data", to be submitted.
}
\author{
Ekaterina Smirnova
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also as \code{\link{FiltLoss}}
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

#Presence absence data matrix
Pres <- Counts
Pres[Pres!=0]<-1

#Proportions  data matrix
Prop <- Counts/rowSums(Counts)

res_sim <- PERFect_sim(X=Prop,  Order=NP,  nbins = 30, col = "red", fill = "green", alpha = 0.05,distr = "sn", quant = c(0.1, 0.25, 0.5), hist_fill =0.2, linecol = "blue")

#permutation perfect colored by FLu values
pvals_Plots(PERFect = res_sim, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.05)

#calculate filtering loss due to removing these taxa
Loss_PERFect_perm <- FL_J(X = Prop, J = names(res_reorder$filtX))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
