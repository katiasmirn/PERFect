\name{PERFect_perm}
\alias{PERFect_perm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Permutation PERFect filtering for microbiome data}
\description{
Permutation filtering of the provided OTU table X at a test level alpha.  Each set of j taxa significance is evaluated by fitting the Skew-Normal, Normal, t or Cauchy distribution to the  sampling distribution obtained by permuted taxa labels.  
}
\usage{
PERFect_perm(X, Order = "NP", Order.user = NULL, normalize = "counts", 
    center = FALSE, quant = c(0.1, 0.25, 0.5), distr = "sn", 
    alpha = 0.1, lag = 3, direction = "left", pvals_sim = NULL, 
    k = 1000, dfl_distr = NULL, nbins = 30, hist = FALSE, col = "red", 
    fill = "green", hist_fill = 0.2, linecol = "blue")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{OTU table, where taxa are columns and samples are rows of the table. It should be a in data frame format with columns corresponding to taxa names.}
  
  \item{Order}{Taxa ordering. The default ordering is the number of occurrences (NP) of the taxa in all samples. Other types of order are p-value ordering, number of connected taxa and weighted number of connected taxa, denoted as "pvals","NC","NCw" respectively. More details about taxa ordering are described in Smirnova et al. User can also specify their preference order with Order.user.}
  
  \item{normalize}{Normalizing taxa count. The default option does not normalize taxa count, but user can convert the OTU table into a proprotion table using the option "prop" or convert it into a presence/absence table using "pres".}
  
  \item{center}{Centering OTU table. The default option does not center the OTU table.}  
  
  \item{quant}{Quantile values used to fit the distribution to log DFL values. The number of quantile values corresponds to the number of parameters in the distribution the data is fitted to. Assuming that at least 50\% of taxa  are not informative, we suggest fitting the  log Skew-Normal distribution by matching the  10\%, 25\% and 50\% percentiles of the log-transformed samples to the Skew-Normal distribution.}
  
  \item{distr}{The type of distribution to fit log DFL values to. While we suggest using Skew-Normal distribution, and set as the defaul distribution, other  choices are available.
  
  "sn" -- Skew-Normal distribution with  3 parameters: location  xi, scale  omega^2 and shape alpha 
  
   "norm" -- Normal distribution with 2 parameters: mean  and standard deviation sd
   
   "t" -- Student t-distribution with 2 parameters: n degrees of freedom and noncentrality ncp
   
   "cauchy" -- Cauchy distribution with 2 parameters: location  and scale }
   
  \item{alpha}{Test level alpha, set to 0.1 by default.}
  
  \item{lag}{Integer width of the rolling window in rolling average (moving mean), set to 3 by default.}
  
  \item{direction}{Character specifying whether the index of the result should be left- or right-aligned or centered compared to the rolling window of observations, set to "left" by default.}
  
  \item{pvals_sim}{Object resulting from simultaneous PERFect with taxa abundance ordering, allowing user to perform Simultaneous PERFect with p-values ordering. Be aware that the choice of distribution for both methods must be the same.}
  
  \item{k}{The number of permutations, set to 10000 by default.}  
  \item{nbins}{Number of bins used to visualize the histogram  of log DFL values, set to 30 by default.}
  \item{col}{Graphical parameter for  color of histogram bars border, set to "red" by default.}
  \item{fill}{Graphical parameter for color of histogram fill, set to "green" by default.}
  \item{hist_fill}{Graphical parameter for intensity of histogram fill, set to 0.2 by default.}
  \item{linecol}{Graphical parameter for the color of the fitted distribution density, set to "blue" by default.}
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
Azzalini, A. (2005). The skew-normal distribution and related multivariate families. Scandinavian Journal of Statistics, 32(2), 159-188.

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
\code{\link{PERFect_sim}}
}
\examples{

data(mock2)

# Proportion data matrix
Prop <- mock2$Prop

# Counts data matrix
Counts <- mock2$Counts

# Perform simultaenous filtering of the data
res_sim <- PERFect_sim(X=Counts)

#order according to p-values
pvals_sim <- pvals_Order(Counts, res_sim)
system.time(res_perm <- PERFect_perm(X=Prop,  Order=pvals_sim))

#permutation perfect colored by FLu values
pvals_Plots(PERFect = res_perm, X = Counts, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.05)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
