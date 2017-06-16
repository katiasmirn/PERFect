\name{PERFect_perm_reorder}
\alias{PERFect_perm_reorder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Permutation PERFect filtering for microbiome data}
\description{
This function  filters  the provided OTU table X at a test level alpha given a fitted object  perfect_perm obtained by running PERFect_perm() function. PERFect_perm_reorder() reavaluates taxa significance p-values for a different taxa ordering. 
}
\usage{
PERFect_perm_reorder(X, Order_alt, res_perm, alpha = 0.05, distr = "sn")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{OTU table, where taxa are columns and samples are rows of the table. Should be in data frame format with columns corresponding to taxa names}
  \item{Order_alt}{Taxa importance ordering}
  \item{res_perm}{Output of PERFect_perm() function}
  \item{alpha}{Test level alpha, set to 0.05 by default}
  \item{distr}{The type of distribution used  in PERFect_perm() function to obtain res_perm object.
  
  "sn" -- Skew-Normal distribution with  3 parameters: location  xi, scale  omega^2 and shape alpha 
  
   "norm" -- Normal distribution with 2 parameters: mean  and standard deviation sd
   
   "t" -- Student t-distribution with 2 parameters: n degrees of freedom and noncentrality ncp
   
   "cauchy" -- Cauchy distribution with 2 parameters: location  and scale 
}
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
See Also as \code{\link{PERFect_perm}}, ~~~
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

#obtain permutation PERFEct results using NP taxa ordering
system.time(res_perm <- PERFect_perm(X=Prop,  Order=NP,   quant = c(0.1,0.25, 0.5), distr = "sn",
                                     alpha = 0.05, k = 10000, 
                                     nbins =30, col = "red", fill = "green", hist_fill = 0.2, linecol = "blue"))

#run PERFEct_sim() function and obtain p-values ordering
res_sim <- PERFect_sim(X=Prop,  Order=NP,  nbins = 30, col = "red", fill = "green", alpha = 0.05,distr = "sn", quant = c(0.1, 0.25, 0.5), hist_fill =0.2, linecol = "blue")

#order according to p-values
pvals_sim <- pvals_Order(Counts, res_sim)

#update perfect_perm object according to p-values ordering
#reorder
res_reorder <- PERFect_perm_reorder(X=Prop,  Order_alt = pvals_sim,  res_perm = res_perm, alpha = 0.05, distr = "sn")

#permutation perfect colored by FLu values
pvals_Plots(PERFect = res_perm, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.05)

#calculate filtering loss due to removing these taxa
Loss_PERFect_perm <- FL_J(X = Prop, J = names(res_reorder$filtX), leave = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
