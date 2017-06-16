\name{PERFect_perm}
\alias{PERFect_perm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Permutation PERFect filtering for microbiome data}
\description{
Permutation filtering of the provided OTU table X at a test level alpha.  Each set of j taxa significance is evaluated by fitting the Skew-Normal, Normal, t or Cauchy distribution to the  sampling distribution obtained by permuted taxa labels.  
}
\usage{
PERFect_perm(X, Order, quant = c(0.1, 0.25, 0.5), distr = "sn", alpha = 0.05, k = 10000, nbins = 30, col = "red", fill = "green", hist_fill = 0.2, linecol = "blue")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{OTU table, where taxa are columns and samples are rows of the table. Should be in data frame format with columns corresponding to taxa names.}

  \item{Order}{Taxa importance ordering according. We suggest using the simultaneous PERFEct p-values.}

  \item{quant}{Quantile values used to fit the distribution to log DFL values. The number of quantile values corresponds to the number of parameters in the distribution the data is fitted to. Assuming that at least 50\% of taxa  are not informative, we suggest fitting the  log Skew-Normal distribution by matching the  10\%, 25\% and 50\% percentiles of the log-transformed samples to the Skew-Normal distribution.}
  \item{distr}{The type of distribution to fit log DFL values to. While we suggest using Skew-Normal distribution, and set as the defaul distribution, other  choices are available.
  
  "sn" -- Skew-Normal distribution with  3 parameters: location  xi, scale  omega^2 and shape alpha 
  
   "norm" -- Normal distribution with 2 parameters: mean  and standard deviation sd
   
   "t" -- Student t-distribution with 2 parameters: n degrees of freedom and noncentrality ncp
   
   "cauchy" -- Cauchy distribution with 2 parameters: location  and scale }
  \item{alpha}{Test level alpha, set to 0.05 by default.}
  \item{k}{The number of permutations}
  \item{nbins}{Number of bins used to visualize the histogram  of log DFL values.}
  \item{col}{Graphical parameter for  color of histogram bars border}
  \item{fill}{Graphical parameter for color of histogram fill}
  \item{hist_fill}{Graphical parameter for intensity of histogram fill}
  \item{linecol}{Graphical parameter for the color of the fitted distribution density}
}
\details{
Filtering is the process of identifying and removing a subset of taxa  according to a particular criterion.  As opposed to the the simultaneous filtering approach, we do not assume that all distributions for each set of taxa  are identical and equal to the distribution of simultaneous filtering. Function PERFect_perm()  filters  the provided OTU table X and outputs a filtered table  that contains signal taxa.  PERFect_perm() calculates differences in filtering loss DFL for each taxon according to the given taxa order. By default, the function fits Skew-Normal distribution to the log-differences in filtering loss but Normal, t, or Cauchy distributions can be also used.
}
\value{
  \item{filtX}{filtered OTU table}
  \item{pvals}{p-values of the test}
  \item{fit}{list of fitted values and further goodness of fit details passed from fitdistr() function for each set of taxa}
  \item{hist}{list of histograms of log differences in filtering loss for each set of taxa}
  \item{est}{estimated distribution parameters}
  \item{dfl_distr}{plot of differences in filtering loss values}
}
\references{
Azzalini, A. (2005). The skew-normal distribution and related multivariate families. Scandinavian Journal of Statistics, 32(2), 159–188.
}
\author{
Ekaterina Smirnova
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also as\code{\link{PERFect_sim}}
}
\examples{
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

#order according to p-values
pvals_sim <- pvals_Order(Counts, res_sim)
system.time(res_perm <- PERFect_perm(X=Prop,  Order=pvals_sim,   quant = c(0.1,0.25, 0.5), distr = "sn",
                                     alpha = 0.05, k = 10000, 
                                     nbins =30, col = "red", fill = "green", hist_fill = 0.2, linecol = "blue"))


#permutation perfect colored by FLu values
pvals_Plots(PERFect = res_perm, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.05)


#calculate filtering loss due to removing these taxa
Loss_PERFect_perm <- FL_J(X = Prop, J = names(res_perm$filtX), leave = TRUE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
