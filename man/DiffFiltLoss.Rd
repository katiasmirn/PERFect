\name{DiffFiltLoss}
\alias{DiffFiltLoss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Difference in filtering loss
}
\description{
This function calculates differences in filtering loss due to removing a set of J taxa sequentially.
}
\usage{
DiffFiltLoss(X, Order_Ind, Plot = TRUE, Taxa_Names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{OTU table, where taxa are columns and samples are rows of the table. It should be a in data frame format with columns corresponding to taxa names.}
  \item{Order_Ind}{Numeric column order corresponding to taxa importance arranegment.}
  \item{Plot}{A binary TRUE/FALSE value. If TRUE, the function returns plot of sequential differences in filtering loss.}
  \item{Taxa_Names}{Optional taxa labels corresponding to the columns ordering given by Order_Ind.}
}
\details{
This function calculates and plots (if Plot = TRUE) differences in filtering loss sequentially  for removing the first j taxa as DFL(j+1) = FL(J_{j+1}) - FL(J_j) for taxa j=1, ..., p. 

}
\value{
  \item{DFL = DFL}{Differences in filtering loss values}
  \item{p_FL= p_FL}{Plot of the differences in filtering loss}

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
\code{\link{FiltLoss}}
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
DFL <- DiffFiltLoss(X=Prop, Order_Ind = Order_Ind, Plot = TRUE, Taxa_Names = NP)

#Differences in filtering loss values
DFL$DFL

#Plot of the differences in filtering loss
DFL$p_FL
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
