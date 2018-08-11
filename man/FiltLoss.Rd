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
FiltLoss(X, Order = "NP", Order.user = NULL, type = "Cumu", Plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{OTU table, where taxa are columns and samples are rows of the table. It should be a in data frame format with columns corresponding to taxa names.}
  \item{Order}{Taxa column order corresponding to taxa importance arrangment. The default ordering is the number of occurrences (NP) of the taxa in all samples. Other types of order are p-value ordering, number of connected taxa and weighted number of connected taxa, denoted as "pvals","NC","NCw" respectively. More details about taxa ordering are described in Smirnova et al. User can also specify their preference order with Order.user.}
  \item{type}{Type of filtering loss calculation. 

"Ind" -- individual taxon's filtering loss FL_u(j)

"Cumu" -- cumulative filtering loss FL(J) due to removing a set of taxa J
}
  \item{Plot}{A binary TRUE/FALSE value. If TRUE, the function returns plot of sequential differences in filtering loss.}
}
\details{
The individual filtering loss due to removing one taxon j is defined as 

FL_u(j)= 1- (||X^T_{-j} X_{-j}||_F^2/||X^TX||_F^2),

where X_{-j} is the matrix X without column corresponding to jth taxon and ||Z||_F is the Frobenious norm of a matrix Z.

The cumulative filtering loss due to removing a set of taxa is defined as

FL(J)= 1- (||X^T_{-J}  X_{-J}||_F^2\||X^TX||_F^2),

where X_{-J} is the n x (p-|J|) dimensional matrix obtained by removing the columns indexed by the set J from the data matrix X.

The cumulative filtering loss is calculated sequentially for each set of taxa J_j, j=1, ..., p.
}
\value{
  \item{FL}{Filtering loss values}
 \item{p_FL}{Plot of filtering loss values}
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
\code{\link{DiffFiltLoss}}
}
\examples{

data(mock2)

# Proportion data matrix
Prop <- mock2$Prop

# Counts data matrix
Counts <- mock2$Counts

#arrange counts in order of increasing number of samples taxa are present in 
NP <- NP_Order(Counts)

#obtain numeric column order corresponding to taxa importance arrangment
Order_Ind <- match(NP, names(Prop))

#Calculate cumulative filtering loss
FL <- FiltLoss(X=Prop, Order = "NP", type = "Cumu", Plot = TRUE)

#Differences in filtering loss values
FL$FL

#Plot of the differences in filtering loss
FL$p_FL
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
