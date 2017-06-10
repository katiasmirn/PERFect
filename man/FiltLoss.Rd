\name{FiltLoss}
\alias{FiltLoss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FiltLoss(X, Order, type = c("Ind", "Cum"), Plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Order}{
%%     ~~Describe \code{Order} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{Plot}{
%%     ~~Describe \code{Plot} here~~
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
function (X, Order, type = c("Ind", "Cum"), Plot = TRUE) 
{
    p <- dim(X)[2]
    Norm_Ratio <- rep(1, p)
    X <- as.matrix(X)
    X <- X[, Order]
    Order_Ind <- 1:length(Order)
    Netw <- t(X) \%*\% X
    for (i in 1:p) {
        if (type == "Cum") {
            Ind <- Order_Ind[-(1:i)]
        }
        else {
            Ind <- Order_Ind[-i]
        }
        Netw_R <- Netw[Ind, Ind]
        Norm_Ratio[i] <- tr(t(Netw_R) \%*\% Netw_R)
    }
    FL <- 1 - Norm_Ratio/tr(t(Netw) \%*\% Netw)
    if (Plot == TRUE) {
        df <- data.frame(Order, rep(1:length(FL)), FL)
        names(df)[2] <- "x"
        Lab <- 1:length(Order)
        df <- cbind(Lab, df)
        p_FL <- ggplot(df) + geom_line(aes(x = Lab, y = FL, group = 1), 
            colour = "dodgerblue3") + theme(panel.background = element_rect(fill = "white"), 
            panel.grid.major = element_line(colour = "grey90"), 
            axis.text.x = element_text(size = 10, colour = "black", 
                angle = 90, hjust = 1)) + ggtitle("") + ylab("Filtering Loss") + 
            xlab("Taxa") + xlim(0, max(df$Lab))
    }
    names(FL) <- colnames(X)
    return(list(FL = FL, p_FL = p_FL))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
