## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(warning=FALSE, message=FALSE, echo = FALSE,fig.align='center')

## ----required_libraries--------------------------------------------------
rm(list=ls())

#library(Matrix)
#library(devtools)
#library(ggplot2)
#library(sn)
#library(fitdistrplus)
#library(psych) 
library(PERFecttest)
library(dirmult)
library(knitr)
library(gridExtra)
library(grid)
library(reshape2)
#library(zoo)
library(kableExtra)

set.seed(12341)
#setwd("C:/Users/Kwee/Dropbox/Quy/PERFect/PERFect-master/R")#change me
#setwd("~/Dropbox/PERFect/RCode/")


## ------------------------------------------------------------------------
data(mock2)
Counts <- mock2$Counts
taxa <- c("Bradyrhizobiaceae.cluster49", "Propionibacterium.acnes", "Gemella.OTU86",
          "Comamonadaceae.cluster57", "Streptococcus.gordonii", "Caulobacter.leidyi",
          "Finegoldia.magna", "Aerococcus.christensenii", "Agromyces.cluster54",
          "Mycoplasma.orale_faucium", "Enterobacteriaceae.cluster31",
          "Lactobacillus.jensenii", "Fusobacterium.cluster48", "Methylophilus.cluster11",
          "Coriobacteriaceae.OTU27", "Prevotella.cluster2",
          "Delftia.acidovorans_lacustris_tsuruhatensis", "Hyphomicrobium.methylovorum",
          "Bifidobacterium.longum_infantis_suis", "Streptococcus.parasanguinis")
Ind <- apply(Counts[,taxa], 2, function(x) which(x!=0))
sampleID <- c(31,78, 39, 5, 9, 76, 22, 105, 34, 82)
data <- Counts[sampleID,c(taxa, mock2$TrueTaxa)]
data <- data[,apply(data, 2, nnzero) >0]
data<- data[,-3]
#sort by abundance
data <- data[, NP_Order(data)]
#rename columns by noise and signal
SignTaxa <- colnames(data) %in% mock2$TrueTaxa
colnames(data)[SignTaxa] <-paste0("S~", 1:sum(SignTaxa),"~")
colnames(data)[!SignTaxa] <-paste0("N~", 1:(dim(data)[2] -sum(SignTaxa)),"~")
rownames(data) <- paste("Sample ", 1:nrow(data))


kable(data[,1:13], 
      format='html', 
      digits=2, 
      row.names=TRUE,
      escape = FALSE,
      padding=-3L) %>% kable_styling(font_size = 10) 

kable(data[,14:20], 
      format='html',
      row.names=TRUE,
      escape = FALSE,
      padding=-3L) %>% kable_styling(font_size = 10)

## ------------------------------------------------------------------------
X <- as.matrix(data)
Netw <- t(as.matrix(X))%*%as.matrix(X)
den <- tr(t(Netw)%*%Netw)
#removing 10th noise taxon
X_R <- X[,-which(colnames(X) == "N~10~")]
#calculate corresponding norm
Netw_R <- t(X_R)%*%X_R
num <-  tr(t(Netw_R)%*%Netw_R)
FL <-  1 - (num/den)
#df <- data.frame(num = num, den = den, FL = FL, log_FL = round(log(FL), 3))
df <- data.frame(FL = FL, log_FL = round(log(FL), 3))

#removing first least abundant signal taxon
X_R <- X[,-which(colnames(X) == "S~1~")]
Netw_R <- t(X_R)%*%X_R
num <-  tr(t(Netw_R)%*%Netw_R)
FL <-  1 - (num/den)
#df <- rbind(df, c(num, den, FL, round(log(FL), 3)))
df <- rbind(df, c(FL, round(log(FL), 3)))

#removing all noise taxa
X_R <- X[,-which(substr(colnames(X), start = 1, stop = 2) %in% "N~")]
Netw_R <- t(X_R)%*%X_R
num <-  tr(t(Netw_R)%*%Netw_R)
FL <-  1 - (num/den)
#df <- rbind(df, c(num, den, FL, round(log(FL), 3)))
df <- rbind(df, c(FL, round(log(FL), 3)))

#removing all noise + first signal taxon
X_R <- X[,-c(which(substr(colnames(X), start = 1, stop = 2) %in% "N~"),
             which(colnames(X) == "S~1~"))]
Netw_R <- t(X_R)%*%X_R
num <-  tr(t(Netw_R)%*%Netw_R)
FL <-  1 - (num/den)
#df <- rbind(df, c(num, den, FL, round(log(FL), 3)))
df <- rbind(df, c(FL, round(log(FL), 3)))

rownames(df) <- c("N~10~", "S~1~", "Noise taxa", "Noise taxa and S~1~")
colnames(df) <- c("Filtering Loss","Log Filtering Loss")
df[,1] <- format(df[,1], width = 1, digits = 5, format = "g")

kable(df, 
      format='markdown',
      row.names=TRUE,
      escape = FALSE,
      padding=-3L) %>% kable_styling("striped",font_size = 12)

## ------------------------------------------------------------------------
#function to output results for each simulation run
resSummary <-function(X, filtX,  time = NA){
  
  rank_pvals = NULL 
  rank_pres = NULL
  ntotal <- dim(X)[2]
  npres <- dim(filtX)[2] #number of taxa left
  
  pfilt <- (ntotal - npres)/ntotal #percentage of taxa filtered

  #combine into results vector
  res <- c(ntotal, npres, pfilt)
  names(res) <- c("ntotal", "npres", "pfilt")
  
  return(list(res = res,  time = time))
  
}

## ------------------------------------------------------------------------

#########################
#quantiles from fit a
##########################
start <- Sys.time()
res_sim <- PERFect_sim(X=data)
end <-  Sys.time()-start
summary_sim <- resSummary(X = data, filtX = res_sim$filtX, time = end)


## ------------------------------------------------------------------------
#add p-values column to the data
pvals  <- c(1, round(res_sim$pvals, 2))
pvals.sort  <- sort.int(pvals, decreasing = TRUE, index.return = TRUE)
#sort data by p-values
data <- data[, pvals.sort$ix]
#save X for illustrating method steps
X = as.matrix(data)
#add p-values to the data
data <- rbind(pvals.sort$x, data)
rownames(data)[1] <- "Sim pvalues"
data[1,1] <- NA


kable(data[1,1:13], 
      format='html', 
      digits=2, 
      row.names=TRUE,
      escape = FALSE,
      padding=-3L) %>% kable_styling(font_size = 10) 

kable(data[1,14:20], 
      format='html',
      row.names=TRUE,
      escape = FALSE,
      padding=-3L) %>% kable_styling(font_size = 10)

## ------------------------------------------------------------------------
Order.vec <- pvals_Order(data, res_sim)
Order_Ind <- rep(1:length(Order.vec))#convert to numeric indicator values
DFL <- DiffFiltLoss(X = X, Order_Ind, Plot = TRUE, Taxa_Names = Order.vec)

#add log DFL values to the data
data <- rbind(c(NA, round(log(DFL$DFL), 2)), data)
rownames(data)[1] <- "log(DFL)"
#rearrange rows to have DFL below sim-pvalues
data <- data[c(2,1, 3:nrow(data)),]
#data

#display DFL and log(DFL) values
df <- rbind(c(NA,DFL$DFL),c(NA, log(DFL$DFL)))
df[1,] <- format(df[1,], width = 1, digits = 3, format = "g")
df[2,] <- round(as.numeric(df[2,]), 2)
rownames(df) <- c("DFL", "log(DFL)")
colnames(df)[1]<- "N~1~"
kable(df[,1:8], 
      format='html',
      row.names=TRUE,
      escape = FALSE,
      padding=-3L) %>% kable_styling(font_size = 12)
kable(df[,9:13], 
      format='html',
      row.names=TRUE,
      escape = FALSE,
      padding=-3L) %>% kable_styling(font_size = 12)
kable(df[,14:20], 
      format='html',
      row.names=TRUE,
      escape = FALSE,
      padding=-3L) %>% kable_styling(font_size = 12)

## ------------------------------------------------------------------------

Perm_j_s_il <- function(j, Netw, k,p, p2 = NULL){
  #Netw = X'X - data matrix with taxa in columns and samples in rows
  #k - number of permutations used  
  #create a list of k*p arraments of orders
  if(is.null(p2)){p2 <- p}
  labl <- sapply(1:k,function(x) NULL)
  labl <- lapply(labl,function(x)  sample(1:p,p2))
  FL_j <- sapply(labl,DiffFiltLoss_j, Netw = Netw, j=j)
  return(list(labl = labl, FL_j = FL_j))
}

#we illustrate computation of the first k=6 permutations for each taxon
p <- dim(X)[2] 
Netw <- t(X)%*%X
full_norm <- tr(t(Netw)%*%Netw)#this is full norm value
#For each taxon j, create a distribution of its DFL's by permuting the labels 
res_all <- lapply(1:(p-1),function(x) x)
FL <- lapply(res_all, function(x) Perm_j_s_il(j = x, Netw =Netw, k=6, p =p, p2 = x+1))
#extract every 2nd subelement for each element of list FL to get FL_j values
FL_j <- lapply(FL, "[[", 2)
#divide by the full matrix norm values to get DFL*(j+1) values
res_pres <- lapply(FL_j, function(x) {x/full_norm})

#extract every 1st subelement for each element of list of selected labels
labl <- lapply(FL, "[[", 1)
labl<- lapply(labl, function(x) matrix(unlist(x), ncol = 6, byrow = FALSE))


#take sets of 5 taxa
T5 <- apply(labl[[4]], 2, function(x) Order.vec[x])
df <- T5
df <- rbind(df,round(log(res_pres[[4]]),2))
rownames(df) <- c(paste0("Taxa ", 1:5),"log(DFL*)")
colnames(df) <- paste0("Permutation ",1:6)

kable(df, 
      format='html',
      row.names=TRUE,
      escape = FALSE,
      padding=-3L) %>% kable_styling(font_size = 12)

## ------------------------------------------------------------------------
#distribution based on k = 1000 permutations
dfl_distr <- sampl_distr(X = X, k=1000) #, sample = FALSE

#distribution fit to skew-normal 
lfl <- lapply(dfl_distr, function(x) log(x[!x==0]))
lp <- list(xi = mean(lfl[[1]]), omega = sd(lfl[[1]]), alpha = 1.5)
fit <- lapply(lfl, function(x) qmedist(x, "sn", probs=c(0.10, 0.25, 0.5), start=lp))  
est <- lapply(fit, function(x) x$estimate)

#histogram of log(DFL^*) values for the 5th taxon
i=4
lfl <- data.frame(log(dfl_distr[[i]][!dfl_distr[[i]]==0]))
    if(length(dfl_distr[[i]][dfl_distr[[i]]==0])>0){
      print(paste("taxon", i, "number of zeroes = ", 
                  length(dfl_distr[[i]][dfl_distr[[i]]==0])))}
    names(lfl) <- c("DFL")
    #plot histogram
    x = "DFL"
    ord_map = aes_string(x = x)
    #hist <- ggplot(lfl, ord_map) + geom_histogram()
    hist <- ggplot(lfl, ord_map) + geom_histogram(bins = 30, aes(y=..density..), 
                                                  col = "red", fill = "green", alpha =0.2)+
      theme(panel.background = element_rect(fill = "white"), 
            panel.grid.major = element_line(colour = "grey90"),
            axis.text.x  = element_text( size=10))+
      ggtitle("") + xlab("log differences in filtering loss") + ylab("Density")

#add skew-normal fit
hist <- hist + stat_function(fun = dsn, 
            args = list(xi = est[[i]][1], omega = est[[i]][2], alpha = est[[i]][3]), colour="blue")
hist   

## ------------------------------------------------------------------------
#all p-values
    for(i in 1:(p-1)){
      pvals[i] <- 1 - psn(x=log(DFL$DFL[i]), 
                      xi = est[[i]][1], omega = est[[i]][2], alpha = est[[i]][3])}

data <- rbind(c(NA, round(pvals,2)), data)
rownames(data)[1] <- "Perm pvalues"


## ------------------------------------------------------------------------
pvals_avg <- rollmean(pvals, k=3, align="left",  fill=NA )
#replace na's with original values
pvals_avg[is.na(pvals_avg)] <- pvals[is.na(pvals_avg)]
  
data <- rbind(c(NA, round(pvals_avg, 2)), data)
rownames(data)[1] <- "Avg Perm pvalues"

#rearrange rows
data <- data[c(3,4,2,1, 5:nrow(data)),]

#nice display
kable(data[1:4,1:13], 
      format='html', 
      digits=2, 
      row.names=TRUE,
      escape = FALSE,
      padding=-3L) %>% kable_styling(font_size = 10) 

kable(data[1:4,14:20], 
      format='html',
      row.names=TRUE,
      escape = FALSE,
      padding=-3L) %>% kable_styling(font_size = 10)

## ------------------------------------------------------------------------
 #select taxa that are kept in the data set at significance level alpha
   Ind <- min(which(pvals_avg <=0.1))
   #if jth DFL is significant, then throw away all taxa 1:j 
   filtX <- X[,-(1:Ind)]
   
# kable(filtX, 
#       format='html',
#       row.names=TRUE,
#       escape = FALSE,
#       padding=-3L) %>% kable_styling(font_size = 12)
#   

