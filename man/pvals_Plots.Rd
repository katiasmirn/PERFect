\name{pvals_Plots}
\alias{pvals_Plots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pvals_Plots(PERFect, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PERFect}{
%%     ~~Describe \code{PERFect} here~~
}
  \item{quantiles}{
%%     ~~Describe \code{quantiles} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
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
function (PERFect, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha = 0.05) 
{
    pvals <- PERFect$pvals
    Order_pvals <- names(pvals)
    taxa_filt <- Order_pvals[!(Order_pvals \%in\% names(PERFect$filtX))]
    taxa <- max(which(taxa_filt \%in\% Order_pvals))
    res_FLu <- FiltLoss(X = Prop[, Order_pvals], Order = Order_pvals, 
        type = "Ind", Plot = TRUE)$FL
    FLu_vals <- res_FLu
    breaks <- quantile(res_FLu, quantiles)
    FLu_vals[FLu_vals <= breaks[1]] <- 1
    FLu_vals[FLu_vals > breaks[1] & FLu_vals <= breaks[2]] <- 2
    FLu_vals[FLu_vals > breaks[2] & FLu_vals <= breaks[3]] <- 3
    FLu_vals[FLu_vals > breaks[3] & FLu_vals <= breaks[4]] <- 4
    FLu_vals[!(FLu_vals \%in\% c(1, 2, 3, 4))] <- 5
    FLu_vals[FLu_vals == 1] <- "0-25\%"
    FLu_vals[FLu_vals == 2] <- "25 - 50\%"
    FLu_vals[FLu_vals == 3] <- "50 - 80\%"
    FLu_vals[FLu_vals == 4] <- "80 - 90\%"
    FLu_vals[FLu_vals == 5] <- "90 - 100\%"
    df <- data.frame(seq(1:length(pvals)), pvals, res_FLu[which(names(res_FLu) \%in\% 
        names(pvals))], FLu_vals[which(names(FLu_vals) \%in\% names(pvals))])
    names(df) <- c("Taxa", "p_value", "FLu", "Quantiles")
    p_pvals <- ggplot(df) + geom_point(aes(x = Taxa, y = p_value, 
        color = Quantiles)) + ggtitle("Permutation PERFect p-values") + 
        theme(panel.background = element_rect(fill = "white"), 
            panel.grid.major = element_line(colour = "grey90"), 
            axis.text.x = element_text(size = 10, colour = "black", 
                angle = 90, hjust = 1)) + guides(color = guide_legend(title = "FLu Quantiles"))
    p_pvals <- p_pvals + geom_hline(yintercept = alpha, color = "red", 
        linetype = "dashed")
    p_pvals <- p_pvals + ggtitle("") + geom_vline(xintercept = taxa, 
        color = "purple", linetype = "dashed")
    return(p_pvals)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
