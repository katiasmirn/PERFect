\name{NP_Order}
\alias{NP_Order}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Taxa importance ordering by the number of occurrences of the taxa in the n samples
}
%\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
%}
\usage{
NP_Order(Counts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Counts}{
OTU COUNTS table, where taxa are columns and samples are rows of the table. It should be a in data frame format with columns corresponding to taxa names.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
 \item{NP}{Taxa names in increasing  order of the number of samples taxa are present in.}
}
\references{
Smirnova, E., Huzurbazar, H., Jafari, F. ``PERFect: permutation  filtration of microbiome data", to be submitted.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%See Also as \code{\link{help}}, ~~~
%}
\examples{

data(mock2)
# Proportion data matrix
Prop <- mock2$Prop

# Counts data matrix
Counts <- mock2$Counts

#arrange counts in order of increasing number of samples taxa are present in 
NP <- NP_Order(Counts)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
