\name{PERFect_perm_reorder}
\alias{PERFect_perm_reorder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Permutation PERFect filtering for microbiome data}
\description{
This function  filters  the provided OTU table X at a test level alpha given a fitted object  perfect_perm obtained by running PERFect_perm() function. PERFect_perm_reorder() reavaluates taxa significance p-values for a different taxa ordering. 
}
\usage{
PERFect_perm_reorder(X, Order = "NP", Order.user = NULL, res_perm, normalize = "counts", 
    center = FALSE, alpha = 0.1, distr = "sn", lag = 3, direction = "left", 
    pvals_sim = NULL) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{OTU table, where taxa are columns and samples are rows of the table. It should be a in data frame format with columns corresponding to taxa names.}
  
  \item{Order}{Taxa ordering. The default ordering is the number of occurrences (NP) of the taxa in all samples. Other types of order are p-value ordering, number of connected taxa and weighted number of connected taxa, denoted as "pvals","NC","NCw" respectively. More details about taxa ordering are described in Smirnova et al. User can also specify their preference order with Order.user.}
  
  \item{res_perm}{Output of PERFect_perm() function.}
  
  \item{normalize}{Normalizing taxa count. The default option does not normalize taxa count, but user can convert the OTU table into a proprotion table using the option "prop" or convert it into a presence/absence table using "pres".}
  
  \item{center}{Centering OTU table. The default option does not center the OTU table.}  
  
  \item{alpha}{Test level alpha, set to 0.1 by default.}
  
  \item{distr}{The type of distribution used in PERFect_perm() function to obtain res_perm object.
  
  "sn" -- Skew-Normal distribution with  3 parameters: location  xi, scale  omega^2 and shape alpha 
  
   "norm" -- Normal distribution with 2 parameters: mean  and standard deviation sd
   
   "t" -- Student t-distribution with 2 parameters: n degrees of freedom and noncentrality ncp
   
   "cauchy" -- Cauchy distribution with 2 parameters: location  and scale 
}
   \item{lag}{Integer width of the rolling window in rolling average (moving mean), set to 3 by default.}
   
    \item{direction}{Character specifying whether the index of the result should be left- or right-aligned or centered compared to the rolling window of observations, set to "left" by default.}
     
  \item{pvals_sim}{Object resulting from simultaneous PERFect with taxa abundance ordering, allowing user to perform Simultaneous PERFect with p-values ordering. Be aware that the choice of distribution for both methods must be the same.}

}
\details{
This function  is designed to save computational time needed to obtain and fit the sapling distribution for each taxon if taxa ordering different from the one used in PERFect_perm() is used. Note, the distribution and OTU table X should match the distribution used in PERFect_perm(). 
}
\value{
  \item{res_perm}{The perfect_perm object updated according to the alternative taxa ordering. All elements in this list are same as in perfect_perm object given by PERFect() function}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ekaterina Smirnova
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{PERFect_perm}}
}
\examples{

data(mock2)

# Proportion data matrix
Prop <- mock2$Prop

# Counts data matrix
Counts <- mock2$Counts

#obtain permutation PERFEct results using NP taxa ordering
system.time(res_perm <- PERFect_perm(X=Prop, k =2))

#run PERFEct_sim() function and obtain p-values ordering
res_sim <- PERFect_sim(X=Prop)

#order according to p-values
pvals_sim <- pvals_Order(Counts, res_sim)

#update perfect_perm object according to p-values ordering
res_reorder <- PERFect_perm_reorder(X=Prop,  Order.user = pvals_sim,  res_perm = res_perm)

#permutation perfect colored by FLu values
pvals_Plots(PERFect = res_perm, X = Counts, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.05)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
