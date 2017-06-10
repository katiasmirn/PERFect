\name{ResultsComparison}
\alias{ResultsComparison}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ResultsComparison(X, Counts, Ab_min = 0.001, rel = FALSE, thresh = 5, prop = FALSE, res_sim, res_perm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Counts}{
%%     ~~Describe \code{Counts} here~~
}
  \item{Ab_min}{
%%     ~~Describe \code{Ab_min} here~~
}
  \item{rel}{
%%     ~~Describe \code{rel} here~~
}
  \item{thresh}{
%%     ~~Describe \code{thresh} here~~
}
  \item{prop}{
%%     ~~Describe \code{prop} here~~
}
  \item{res_sim}{
%%     ~~Describe \code{res_sim} here~~
}
  \item{res_perm}{
%%     ~~Describe \code{res_perm} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, Counts, Ab_min = 0.001, rel = FALSE, thresh = 5, 
    prop = FALSE, res_sim, res_perm) 
{
    Filt_R1 <- TraditR1(X = Counts, rel = rel, thresh = thresh)
    Filt_R2 <- TraditR2(X = Counts, Ab_min = Ab_min, prop = prop)
    Trad_R1 <- names(Filt_R1)
    Trad_R2 <- names(Filt_R2)
    Loss_Trad_R1 <- FL_J(X = X, J = Trad_R1, leave = TRUE)
    Loss_Trad_R2 <- FL_J(X = X, J = Trad_R2, leave = TRUE)
    PERFect_sim_taxa <- names(res_sim$filtX)
    PERFect_perm_taxa <- names(res_perm$filtX)
    Loss_PERFect_sim <- FL_J(X = X, J = PERFect_sim_taxa, leave = TRUE)
    Loss_PERFect_perm <- FL_J(X = X, J = PERFect_perm_taxa, leave = TRUE)
    Res <- cbind(rbind(length(Trad_R1), length(Trad_R2), length(PERFect_sim_taxa), 
        length(PERFect_perm_taxa)), rbind(Loss_Trad_R1, Loss_Trad_R2, 
        Loss_PERFect_sim, Loss_PERFect_perm))
    colnames(Res) <- c("N_taxa", "FL")
    rownames(Res) <- c("Trad_R1", "Trad_R2", "PERFect_sim", "PERFect_perm")
    return(list(Res = Res, Filt_R1 = Filt_R1, Filt_R2 = Filt_R2, 
        res_sim = res_sim, res_perm = res_perm))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
