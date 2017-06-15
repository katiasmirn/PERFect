\name{PERFect_sim}
\alias{PERFect_sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation PERFect filtering for microbiome data}

\usage{
PERFect_sim(X, Order, nbins = 30, quant = c(0.25, 0.5), distr = c("sn", "norm","t", "cauchy"), alpha = 0.05, col = "red", fill = "green", hist_fill = 0.2, linecol = "blue")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{OTU table, where taxa are columns and samples are rows of the table. Should be in data frame format with columns corresponding to taxa names.}
  \item{Order}{Taxa importance ordering according. We suggest using the number of occurrences of the taxa in the $n$ samples $NP$.}
  \item{nbins}{Number of bins used to visualize the histogram  of log $DFL$ values}
  \item{quant}{Quantile values used to fit the distribution to log $DFL$ values. The number of quantile values corresponds to the number of parameters in the distribution the data is fitted to. Assuming that at least $50\%$ of taxa  are not informative, we suggest fitting the  log Skew-Normal distribution by matching the  $10\%$, $25\%$ and $50\%$ percentiles of the log-transformed samples to the Skew-Normal distribution.}
  \item{distr}{The type of distribution to fit log $DFL$ values to. While we suggest using Skew-Normal distribution, and set as the defaul distribution, other  choices are available.\\
  "sn" -- Skew-Normal distribution with  location parameter $\xi$, scale parameter $\omega^2$ and shape parameter $\alpha$ denoted by $\SN(\xi,\omega^2, \alpha)$ \\
   "norm" -- Normal distribution with mean $\mu$ and variance $\sigma^2$ denoted by $N(\mu, \sigma^2)$\\
   "t" -- Student t-distribution with $n$ degrees of freedom and noncentrality parameter $ncp$\\
   "cauchy" -- Cauchy distribution with location parameter $\mu$ and scale parameter $\sigma^2$}
  \item{alpha}{Test level $\alpha$, set to $0.05$ by default.}
  \item{col}{Graphical parameter for  color of histogram bars border}
  \item{fill}{Graphical parameter for color of histogram fill}
  \item{hist_fill}{Graphical parameter for intensity of histogram fill}
  \item{linecol}{Graphical parameter for the color of the fitted distribution density}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  \item{filtX }{filtered OTU table}
  \item{pvals }{p-values of the test}
  \item{DFL}{differences in filtering loss values}
  \item{fit}{fitted values and further goodness of fit details passed from fitdistr() function}
  \item{hist}{histogram of log differences in filtering loss}
  \item{est}{estimated distribution parameters}
  \item{pDFL}{plot of differences in filtering loss values}
}
\description{
Filtering is the process of identifying and removing a subset of taxa $X_J=\{\x_j\}_{j\in J}$, where $J\subset\{1,\ldots,p\}$, according to a particular criterion. The original data matrix $X$ can be written (after re-arranging some columns) as $X=(X_{J},X_{-J})$, where $X_J$ is the $n\times |J|$ dimensional matrix containing the taxa that are removed and $X_{-J}$ is the $n\times (p-|J|)$ dimensional matrix containing the taxa that are retained for further analysis. Here $|A|$ denotes the cardinality of the set of indices $A$. Function PERFect_sim()  filters  the provided OTU table X and outputs a filtered table $X_{-J}$ that contains signal taxa.  PERFect_sim() calculates differences in filtering loss $DFL$ for each taxon according to the given taxa order. By default, the function fits Skew-Normal distribution to the log-differences in filtering loss but Normal, t, or Cauchy distributions can be also used.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ekaterina Smirnova
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also as \code{\link{PERFect_perm}}
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

res_sim <- PERFect_sim(X=Prop,  Order=NP,  nbins = 30, col = "red", fill = "green", alpha = 0.05,distr = "sn",
                    quant = c(0.1, 0.25, 0.5), hist_fill =0.2, linecol = "blue")
#order according to p-values
pvals_sim <- pvals_Order(Counts, res_sim)
system.time(res_perm <- PERFect_perm(X=Prop,  Order=pvals_sim,   quant = c(0.1,0.25, 0.5), distr = "sn",
                                     alpha = 0.05, k = 10000, 
                                     nbins =30, col = "red", fill = "green", hist_fill = 0.2, linecol = "blue"))

save.image(file=paste(pathtoplots, "res_Data2.RData", sep=""))

#permutation perfect colored by FLu values
pvals_Plots(PERFect = res_perm, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.05)

#reorder
res_reorder <- PERFect_perm_reorder(X=Prop,  Order_alt = NP,  res_perm = res_perm, alpha = 0.05, distr = "sn")

#calculate filtering loss due to removing these taxa
Loss_PERFect_perm <- FL_J(X = Prop, J = names(res_reorder$filtX), leave = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
