\name{DiffFiltLoss}
\alias{DiffFiltLoss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DiffFiltLoss(X, Order_Ind, Plot = TRUE, Taxa_Names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Order_Ind}{
%%     ~~Describe \code{Order_Ind} here~~
}
  \item{Plot}{
%%     ~~Describe \code{Plot} here~~
}
  \item{Taxa_Names}{
%%     ~~Describe \code{Taxa_Names} here~~
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
function (X, Order_Ind, Plot = TRUE, Taxa_Names = NULL) 
{
    p <- dim(X)[2]
    DFL <- rep(1, p - 1)
    X <- as.matrix(X)
    Netw <- t(X) \%*\% X
    for (j in 1:(p - 1)) {
        DFL[j] <- DiffFiltLoss_j(Order_Ind, Netw, j)
    }
    DFL <- DFL/tr(t(Netw) \%*\% Netw)
    if (Plot == TRUE) {
        df <- data.frame(Order_Ind[-1], rep(1:length(DFL)), DFL)
        names(df)[2] <- "x"
        Lab <- 1:length(Order_Ind[-1])
        df <- cbind(Lab, df)
        p_FL <- ggplot(df) + geom_line(aes(x = Lab, y = DFL, 
            group = 1), colour = "dodgerblue3") + theme(panel.background = element_rect(fill = "white"), 
            panel.grid.major = element_line(colour = "grey90"), 
            axis.text.x = element_text(size = 10, colour = "black", 
                angle = 90, hjust = 1)) + ggtitle("") + ylab("Differences in Filtering Loss") + 
            xlab("Taxa") + xlim(0, max(df$Lab))
    }
    if (!is.null(Taxa_Names)) {
        names(DFL) <- Taxa_Names[-1]
    }
    return(list(DFL = DFL, p_FL = p_FL))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
