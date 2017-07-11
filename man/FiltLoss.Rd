\name{FiltLoss}
\alias{FiltLoss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filtering Loss
}
\description{
Sequential filtering loss calculation for removing a set of J_j taxa for J= 1, ..., p.
}
\usage{
FiltLoss(X, Order, type = c("Ind", "Cum"), Plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
An n x p OTU table, where p taxa are columns and n samples are rows of the table. Should be in data frame format with columns corresponding to taxa names.
}
  \item{Order}{
Taxa column order  corresponding to taxa importance arrangment. Can be either character vector corresponding to taxa names recorded in column names of X or a numeric vector with  column indices.
}
  \item{type}{
Type of filtering loss calculation. 

"Ind" -- individual taxon's filtering loss FL_u(j)

"Cum" -- cumulative filtering loss FL(J) due to removing a set of taxa J
}
  \item{Plot}{
A binary TRUE/FALSE value. If TRUE, the function returns plot of sequential differences in filtering loss.
}
}
\details{
The individual filtering loss due to removing one taxon j is defined as 

FL_u(j)= 1- (||X^T_{-j} X_{-j}||_F^2/||X^TX||_F^2),

where X_{-j} is the matrix X without column corresponding to jth taxon and ||Z||_F is the Frobenious norm of a matrix Z.

The cumulative filtering loss due to removing a set of taxa 

FL(J)= 1- (||X^T_{-J}  X_{-J}||_F^2\||X^TX||_F^2),

where X_{-J} is the n x (p-|J|) dimensional matrix obtained by removing the columns indexed by the set J from the data matrix X.

The cumulative filtering loss is calculated sequentially for each set of taxa J_j, j=1, ..., p.
}
\value{
  \item{FL}{Description of 'comp1'}
 \item{p_FL}{Description of 'comp2'}
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
See Also as \code{\link{DiffFiltLoss}}
}
\examples{
Prop <- mock2$Prop
Counts <- mock2$Counts
#arrange counts in order of increasing number of samples taxa are present in 
NP <- NP_Order(Counts)
#Calculate cumulative filtering loss
FL <- FiltLoss(X=Prop, Order = Order, type = "Cum", Plot = TRUE)
#Differences in filtering loss values
FL$FL
#Plot of the differences in filtering loss
FL$p_FL
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
