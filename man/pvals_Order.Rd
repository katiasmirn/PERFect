\name{pvals_Order}
\alias{pvals_Order}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Taxa importance ordering by PERFect p-values
}
\description{
This function orders taxa by increasing significance of simultaneous PERFect p-values.
}
\usage{
pvals_Order(Counts, res_sim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Counts}{OTU COUNTS table, where taxa are columns and samples are rows of the table. It should be a in data frame format with columns corresponding to taxa names.}
  \item{res_sim}{Output of PERFect_sim() function.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
 \item{Order_pvals}{Taxa names in increasing order of p-values significance.}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Ekaterina Smirnova
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{PERFect_sim}}, \code{\link{PERFect_perm}}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
