## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(warning=FALSE, message=FALSE, echo = FALSE,fig.align='center')

## ----required_libraries--------------------------------------------------

rm(list=ls())
#packages for salter data
library(decontam)
library(phyloseq)
library(reshape2)#;packageVersion("reshape2")

#library(Matrix)
#library(devtools)
#require(ggplot2)
#require(sn)
#require(fitdistrplus)
#require(psych) 
library(PERFecttest)
library(dirmult)
library(HMP)
library(knitr)
library(kableExtra)
library(gridExtra)
library(grid)
library(reshape2)
#library(zoo)


set.seed(12341)
#setwd("~/Dropbox/PERFect/RCode/")
#setwd("C:/Users/Kwee/Dropbox/Quy/PERFect/PERFect-master/R")
#source("~/Dropbox/PERFect/RCode/SimData.R")
#pathtodata <- "~/Dropbox/PERFect/RCode/Responses/Rdata/"

## ------------------------------------------------------------------------
path.ampli <- "C:/Users/Quy/Dropbox/Quy/PERFect/PERFecttest/data" # CHANGE ME
load(file.path(path.ampli, "salter1.rda"))
df.ampli <- salter1
rownames(df.ampli) <- df.ampli$Run.accession

## ------------------------------------------------------------------------
load(file.path(path.ampli, "st.ampli.rda"))
load(file.path(path.ampli, "tax.ampli.rda"))
ft.ampli <- sweep(st.ampli, 1, rowSums(st.ampli), "/")
df.ampli$Dilution.number[df.ampli$Dilution.number == "0 (original culture)"] <- "0"
df.ampli$Dilution.number[df.ampli$Dilution.number == "Negative control"] <- "Neg"
conc.dict <- c("0"=1e3, "0 (original culture)"=1e3, "1"=1e2, "2"=1e1, "3"=1, "4"=1, "5"=1, "Neg"=1)
df.ampli$conc <- conc.dict[df.ampli$Dilution.number]
#identical(rownames(df.ampli), rownames(st.ampli)) # TRUE
ps.ampli <- phyloseq(otu_table(st.ampli, taxa_are_rows=FALSE), tax_table(tax.ampli),sample_data(df.ampli))

## ------------------------------------------------------------------------
#all non-Salmonealla reads in each sample are contaminants
#these are the First 3  real SVs from the S. bongori strain
true <- "Salmonella" == unname(tax.ampli)[,6]
#make a vector of taxa cont and  id's
type <-  rep("contam", dim(tax.ampli)[1])
type[1:3] <- "true"

taxaInfo <- data.frame(unname(tax.ampli))
taxa <- rownames(taxaInfo)
#change to numeric id  
names(taxaInfo) <- c("kindgdom", "phylum", "class", "order",
                    "family", "genus")
taxaInfo <- cbind(taxa, type, taxaInfo)

#arrange taxa in prevalence order NP
#head is least dominant, tail is most dominant
Salter.counts <- data.frame(otu_table(ps.ampli))
names(Salter.counts) <- taxa
NP <- NP_Order(Salter.counts)
taxaInfo <- taxaInfo[match(NP, taxaInfo$taxa),]
taxaInfo$taxa <- factor(taxaInfo$taxa, levels = NP)
#all 3 true features are the most dominant in the data

## ------------------------------------------------------------------------
#dim(Salter.counts)
Counts <- Salter.counts

## ------------------------------------------------------------------------
#function to output results for each simulation run
resSummary <-function(X, filtX, taxaInfo,  time = NA){
  rank_pvals = NULL 
  rank_pres = NULL
  ntotal <- dim(X)[2]
  npres <- dim(filtX)[2] #number of taxa left
  pfilt <- (ntotal - npres)/ntotal #percentage of taxa filtered
  ntrue <- sum(colnames(filtX) %in% c(1:3))#compare with true taxa
  ncont <- length(taxaInfo$taxa) - 3
  perccont <- (npres - ntrue)/ncont
  #combine into results vector
  res <- c(ntotal, npres, pfilt, ntrue, perccont)
  names(res) <- c("ntotal", "npres", "pfilt", "ntrue", "perccont")
  return(list(res = res,  time = time))
}

## ------------------------------------------------------------------------

#########################
#quantiles from fit a
##########################
start <- Sys.time()
res_sim_sn_a <- PERFect_sim(X=Counts,  Order="NP",  nbins = 30, col = "red", 
                       fill = "green", alpha = 0.1, distr = "sn", 
                       quant = c(0.05,0.10, 0.25), hist_fill =0.2, linecol = "blue",
                       lag = 3, direction ="left")
end <-  Sys.time()-start
summary_sim_sn_a <- resSummary(X = Counts, filtX = res_sim_sn_a$filtX, 
                          taxaInfo = taxaInfo,  time = end)

#apply two more p-values thresholds
filtX <- filt_pval(X = Counts, pvals =res_sim_sn_a$pvals, alpha = 0.05)
summary_sim_sn_a_0.05 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo, time = NA)

filtX <- filt_pval(X = Counts, pvals =res_sim_sn_a$pvals, alpha = 0.15)
summary_sim_sn_a_0.15 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)


#########################
#quantiles from fit b
##########################
start <- Sys.time()
res_sim_sn_b <- PERFect_sim(X=Counts,  Order="NP",  nbins = 30, col = "red", 
                       fill = "green", alpha = 0.1, distr = "sn", 
                       quant = c(0.1,0.25, 0.40), hist_fill =0.2, linecol = "blue",
                       lag = 3, direction ="left")
end <-  Sys.time()-start
#results for sumultaneous PERFect with NP ordering
summary_sim_sn_b <- resSummary(X = Counts, filtX = res_sim_sn_b$filtX, 
                          taxaInfo = taxaInfo,  time = NA)

#apply two more p-values thresholds
filtX <- filt_pval(X = Counts, pvals =res_sim_sn_b$pvals, alpha = 0.05)
summary_sim_sn_b_0.05 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)

filtX <- filt_pval(X = Counts, pvals =res_sim_sn_b$pvals, alpha = 0.15)
summary_sim_sn_b_0.15 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)
#########################
#quantiles from fit c
##########################
start <- Sys.time()
res_sim_sn_c <- PERFect_sim(X=Counts,  Order="NP",  nbins = 30, col = "red", 
                       fill = "green", alpha = 0.1, distr = "sn", 
                       quant = c(0.1,0.25, 0.50), hist_fill =0.2, linecol = "blue",
                       lag = 3, direction ="left")
end <-  Sys.time()-start
#results for sumultaneous PERFect with NP ordering
summary_sim_sn_c <- resSummary(X = Counts, filtX = res_sim_sn_c$filtX, 
                          taxaInfo = taxaInfo,  time = NA)

#apply two more p-values thresholds
filtX <- filt_pval(X = Counts, pvals =res_sim_sn_c$pvals, alpha = 0.05)
summary_sim_sn_c_0.05 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)

filtX <- filt_pval(X = Counts, pvals =res_sim_sn_c$pvals, alpha = 0.15)
summary_sim_sn_c_0.15 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)

#########################
#quantiles from fit d
##########################
start <- Sys.time()
res_sim_sn_d <- PERFect_sim(X=Counts,  Order="NP",  nbins = 30, col = "red", 
                       fill = "green", alpha = 0.1, distr = "sn", 
                       quant = c(0.20,0.30, 0.60), hist_fill =0.2, linecol = "blue",
                       lag = 3, direction ="left")
end <-  Sys.time()-start
#results for sumultaneous PERFect with NP ordering
summary_sim_sn_d <- resSummary(X = Counts, filtX = res_sim_sn_d$filtX, 
                          taxaInfo = taxaInfo,  time = NA)


#apply two more p-values thresholds
filtX <- filt_pval(X = Counts, pvals =res_sim_sn_d$pvals, alpha = 0.05)
summary_sim_sn_d_0.05 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)

filtX <- filt_pval(X = Counts, pvals =res_sim_sn_d$pvals, alpha = 0.15)
summary_sim_sn_d_0.15 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)

## ------------------------------------------------------------------------

#########################
#quantiles from fit a
##########################
start <- Sys.time()
res_sim_sn_a_pvals <- PERFect_sim(X=Counts,  Order="pvals",  nbins = 30, col = "red", 
                       fill = "green", alpha = 0.1, distr = "sn", pvals_sim = res_sim_sn_a,
                       quant = c(0.05,0.10, 0.25), hist_fill =0.2, linecol = "blue",
                       lag = 3, direction ="left")
end <-  Sys.time()-start
summary_sim_sn_a_pvals <- resSummary(X = Counts, filtX = res_sim_sn_a_pvals$filtX, 
                          taxaInfo = taxaInfo,  time = end)

#apply two more p-values thresholds
filtX <- filt_pval(X = Counts, pvals =res_sim_sn_a_pvals$pvals, alpha = 0.05)
summary_sim_sn_a_0.05_pvals <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)

filtX <- filt_pval(X = Counts, pvals =res_sim_sn_a_pvals$pvals, alpha = 0.15)
summary_sim_sn_a_0.15_pvals <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)


#########################
#quantiles from fit b
##########################
start <- Sys.time()
res_sim_sn_b_pvals <- PERFect_sim(X=Counts,  Order="pvals",  nbins = 30, col = "red", 
                       fill = "green", alpha = 0.1, distr = "sn", pvals_sim = res_sim_sn_b,
                       quant = c(0.1,0.25, 0.40), hist_fill =0.2, linecol = "blue",
                       lag = 3, direction ="left")
end <-  Sys.time()-start
#results for sumultaneous PERFect with NP ordering
summary_sim_sn_b_pvals <- resSummary(X = Counts, filtX = res_sim_sn_b_pvals$filtX, 
                          taxaInfo = taxaInfo,  time = NA)

#apply two more p-values thresholds
filtX <- filt_pval(X = Counts, pvals =res_sim_sn_b_pvals$pvals, alpha = 0.05)
summary_sim_sn_b_0.05_pvals <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)

filtX <- filt_pval(X = Counts, pvals =res_sim_sn_b_pvals$pvals, alpha = 0.15)
summary_sim_sn_b_0.15_pvals <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)
#########################
#quantiles from fit c
##########################
start <- Sys.time()
res_sim_sn_c_pvals <- PERFect_sim(X=Counts,  Order="pvals",  nbins = 30, col = "red", 
                       fill = "green", alpha = 0.1, distr = "sn", pvals_sim = res_sim_sn_c,
                       quant = c(0.1,0.25, 0.50), hist_fill =0.2, linecol = "blue",
                       lag = 3, direction ="left")
end <-  Sys.time()-start
#results for sumultaneous PERFect with NP ordering
summary_sim_sn_c_pvals <- resSummary(X = Counts, filtX = res_sim_sn_c_pvals$filtX, 
                          taxaInfo = taxaInfo,  time = NA)

#apply two more p-values thresholds
filtX <- filt_pval(X = Counts, pvals =res_sim_sn_c_pvals$pvals, alpha = 0.05)
summary_sim_sn_c_0.05_pvals <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)

filtX <- filt_pval(X = Counts, pvals =res_sim_sn_c_pvals$pvals, alpha = 0.15)
summary_sim_sn_c_0.15_pvals <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)

#########################
#quantiles from fit d
##########################
start <- Sys.time()
res_sim_sn_d_pvals <- PERFect_sim(X=Counts,  Order="pvals",  nbins = 30, col = "red", 
                       fill = "green", alpha = 0.1, distr = "sn", pvals_sim = res_sim_sn_b,
                       quant = c(0.20,0.30, 0.60), hist_fill =0.2, linecol = "blue",
                       lag = 3, direction ="left")
end <-  Sys.time()-start
#results for sumultaneous PERFect with NP ordering
summary_sim_sn_d_pvals <- resSummary(X = Counts, filtX = res_sim_sn_d_pvals$filtX, 
                          taxaInfo = taxaInfo,  time = NA)


#apply two more p-values thresholds
filtX <- filt_pval(X = Counts, pvals =res_sim_sn_d_pvals$pvals, alpha = 0.05)
summary_sim_sn_d_0.05_pvals <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)

filtX <- filt_pval(X = Counts, pvals =res_sim_sn_d_pvals$pvals, alpha = 0.15)
summary_sim_sn_d_0.15_pvals <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)

## ------------------------------------------------------------------------
start <- Sys.time()
res_perm_a <- PERFect_perm(X=Counts,  Order="NP",  nbins = 30, col = "red", k = 2,
                       fill = "green", alpha = 0.1, distr = "sn", 
                       quant = c(0.05,0.10, 0.25), hist_fill =0.2, linecol = "blue",
                        lag = 3, direction ="left")
end <-  Sys.time()-start
#save results
#saveRDS(res_perm_a, file = paste0(pathtodata, "res_perm_Salter_a.RDS"))
#saveRDS(end, file = paste0(pathtodata, "time_Salter_a.RDS"))
#read results for faster processing
#res_perm_a <- readRDS(file = paste0(pathtodata, "res_perm_Salter_a.RDS"))
#end <- readRDS(file = paste0(pathtodata, "time_Salter_a.RDS"))
#very normal looking for these data
#res_perm_a$hist[[200]]
#results for permutation PERFect with NP ordering
summary_perm_np_a <- resSummary(X = Counts, filtX = res_perm_a$filtX, 
                          taxaInfo = taxaInfo,  time = end)

#apply two more p-values thresholds
filtX <- filt_pval(X = Counts, pvals =res_perm_a$pvals, alpha = 0.05)
summary_perm_np_a_0.05 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)

filtX <- filt_pval(X = Counts, pvals =res_perm_a$pvals, alpha = 0.15)
summary_perm_np_a_0.15 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo, time = NA)
#########################
#quantiles from fit b
##########################
start <- Sys.time()
res_perm_b <- PERFect_perm(X=Counts,  Order="NP",  nbins = 30, col = "red", 
                       fill = "green", alpha = 0.1, distr = "sn", 
                       quant = c(0.1,0.25, 0.40), hist_fill =0.2, linecol = "blue",
                       dfl_distr = res_perm_a$dfl_distr,
                       lag = 3, direction ="left")
end <-  Sys.time()-start
#saveRDS(res_perm_b, file = paste0(pathtodata, "res_perm_Salter_b.RDS"))
#read results for faster processing
#res_perm_b <- readRDS(file = paste0(pathtodata, "res_perm_Salter_b.RDS"))
#results for sumultaneous PERFect with NP ordering
summary_perm_np_b <- resSummary(X = Counts, filtX = res_perm_b$filtX, 
                          taxaInfo = taxaInfo,  time = NA)

#apply two more p-values thresholds
filtX <- filt_pval(X = Counts, pvals =res_perm_b$pvals, alpha = 0.05)
summary_perm_np_b_0.05 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)

filtX <- filt_pval(X = Counts, pvals =res_perm_b$pvals, alpha = 0.15)
summary_perm_np_b_0.15 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)
#########################
#quantiles from fit c
##########################
start <- Sys.time()
res_perm_c <- PERFect_perm(X=Counts,  Order="NP",  nbins = 30, col = "red", 
                       fill = "green", alpha = 0.1, distr = "sn", 
                       quant = c(0.1,0.25, 0.50), hist_fill =0.2, linecol = "blue",
                       dfl_distr = res_perm_a$dfl_distr, lag = 3, direction ="left")
end <-  Sys.time()-start
#saveRDS(res_perm_c, file = paste0(pathtodata, "res_perm_Salter_c.RDS"))

#read results for faster processing
#res_perm_c <- readRDS(file = paste0(pathtodata, "res_perm_Salter_c.RDS"))
#results for sumultaneous PERFect with NP ordering
summary_perm_np_c <- resSummary(X = Counts, filtX = res_perm_c$filtX, 
                          taxaInfo = taxaInfo,  time = NA)

#apply two more p-values thresholds
filtX <- filt_pval(X = Counts, pvals =res_perm_c$pvals, alpha = 0.05)
summary_perm_np_c_0.05 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)

filtX <- filt_pval(X = Counts, pvals =res_perm_c$pvals, alpha = 0.15)
summary_perm_np_c_0.15 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)
#########################
#quantiles from fit d
##########################
start <- Sys.time()
res_perm_d <- PERFect_perm(X=Counts,  Order="NP",  nbins = 30, col = "red", 
                       fill = "green", alpha = 0.1, distr = "sn", 
                       quant = c(0.20,0.30, 0.60), hist_fill =0.2, linecol = "blue",
                       dfl_distr = res_perm_a$dfl_distr, lag = 3, direction ="left")
end <-  Sys.time()-start

#saveRDS(res_perm_d, file = paste0(pathtodata, "res_perm_Salter_d.RDS"))


#read results for faster processing
#res_perm_d <- readRDS(file = paste0(pathtodata, "res_perm_Salter_d.RDS"))

#results for sumultaneous PERFect with NP ordering
summary_perm_np_d <- resSummary(X = Counts, filtX = res_perm_d$filtX, 
                          taxaInfo = taxaInfo,  time = NA)

#apply two more p-values thresholds
filtX <- filt_pval(X = Counts, pvals =res_perm_d$pvals, alpha = 0.05)
summary_perm_np_d_0.05 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)

filtX <- filt_pval(X = Counts, pvals =res_perm_d$pvals, alpha = 0.15)
summary_perm_np_d_0.15 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)

## ------------------------------------------------------------------------

start <- Sys.time()
res_perm_pvals_a <- PERFect_perm_reorder(X=Counts,  Order = "pvals",  
                                       pvals_sim = res_sim_sn_a,
                                       res_perm = res_perm_a, alpha = 0.1, distr = "sn",
                                       lag = 3, direction ="left")
end <-  Sys.time()-start
#results for permutation PERFect with NP ordering
summary_perm_pvals_a <- resSummary(X = Counts, filtX = res_perm_pvals_a$filtX, 
                          taxaInfo = taxaInfo,  time = NA)

#apply two more p-values thresholds
filtX <- filt_pval(X = Counts, pvals =res_perm_pvals_a$pvals, alpha = 0.05)
summary_perm_pvals_a_0.05 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)

filtX <- filt_pval(X = Counts, pvals =res_perm_pvals_a$pvals, alpha = 0.15)
summary_perm_pvals_a_0.15 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)
#########################
#quantiles from fit b
##########################

start <- Sys.time()
res_perm_pvals_b <- PERFect_perm_reorder(X=Counts,  Order = "pvals",  
                                       pvals_sim = res_sim_sn_b,  
                                       res_perm = res_perm_b, alpha = 0.1, distr = "sn",
                                       lag = 3, direction ="left")
end <-  Sys.time()-start
#results for permutation PERFect with NP ordering
summary_perm_pvals_b <- resSummary(X = Counts, filtX = res_perm_pvals_b$filtX, 
                          taxaInfo = taxaInfo,  time = NA)

#apply two more p-values thresholds
filtX <- filt_pval(X = Counts, pvals =res_perm_pvals_b$pvals, alpha = 0.05)
summary_perm_pvals_b_0.05 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)

filtX <- filt_pval(X = Counts, pvals =res_perm_pvals_b$pvals, alpha = 0.15)
summary_perm_pvals_b_0.15 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo, time = NA)
#########################
#quantiles from fit c
##########################

start <- Sys.time()
res_perm_pvals_c <- PERFect_perm_reorder(X=Counts,  Order = "pvals",  
                                       pvals_sim = res_sim_sn_c, 
                                       res_perm = res_perm_c, alpha = 0.1, distr = "sn",
                                       lag = 3, direction ="left")
end <-  Sys.time()-start
#results for permutation PERFect with NP ordering
summary_perm_pvals_c <- resSummary(X = Counts, filtX = res_perm_pvals_c$filtX, 
                          taxaInfo = taxaInfo,  time = NA)

#apply two more p-values thresholds
filtX <- filt_pval(X = Counts, pvals =res_perm_pvals_c$pvals, alpha = 0.05)
summary_perm_pvals_c_0.05 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)

filtX <- filt_pval(X = Counts, pvals =res_perm_pvals_c$pvals, alpha = 0.15)
summary_perm_pvals_c_0.15 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)

#########################
#quantiles from fit d
##########################

start <- Sys.time()
res_perm_pvals_d <- PERFect_perm_reorder(X=Counts,  Order = "pvals",  
                                       pvals_sim = res_sim_sn_d,
                                       res_perm = res_perm_d, alpha = 0.1, distr = "sn",
                                       lag = 3, direction ="left")
end <-  Sys.time()-start
#results for permutation PERFect with NP ordering
summary_perm_pvals_d <- resSummary(X = Counts, filtX = res_perm_pvals_d$filtX, 
                          taxaInfo = taxaInfo,  time = NA)

#apply two more p-values thresholds
filtX <- filt_pval(X = Counts, pvals =res_perm_pvals_d$pvals, alpha = 0.05)
summary_perm_pvals_d_0.05 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)

filtX <- filt_pval(X = Counts, pvals =res_perm_pvals_d$pvals, alpha = 0.15)
summary_perm_pvals_d_0.15 <- resSummary(X = Counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = NA)


## ------------------------------------------------------------------------
start <- Sys.time()
ampli.min <- isContaminant(ps.ampli, method="frequency", conc="conc", batch="Processing.Institute", batch.combine="minimum", normalize=TRUE)
end <-  Sys.time()-start
ampli.pool <- isContaminant(ps.ampli, method="frequency", conc="conc", normalize=TRUE)

## ------------------------------------------------------------------------
#Plot the removal of contaminants as a function of the classification threshold: a) 0.05; b) 0.1; c) 0.1; d) 0.2. 
#contaminant taxa according to decontam package
#which(ampli.min$contaminant == TRUE) this corresponds to threshold value of t = 0.1

# Determine the total number of contaminant (i.e. non-Salmonealla) reads in each sample
threshs <- c(0.05, 0.1, 0.2, 0.3) #here 0.1 threshold is the default

decontam_pool <- list()
decontam_min <- list()

#otu table defined previously and renamed according to taxa ids instead of DNA fragments
#Salter.counts <- data.frame(otu_table(ps.ampli))
#names(Salter.counts) <- taxa

for (i in 1:length(threshs)){
  
t <- threshs[i]
cont <- ampli.pool$p<t
cont[is.na(cont)] <- FALSE
filtX <- Salter.counts[,!cont]  
pvals = ampli.pool$p
names(pvals) <- taxa
#results for decontam filtering with NP ordering
decontam_pool[[i]] <- resSummary(X = Salter.counts, filtX = filtX, 
                          taxaInfo = taxaInfo,  time = end)

cont <- ampli.min$p<t
cont[is.na(cont)] <- FALSE
filtX <- Salter.counts[,!cont]  
pvals = ampli.min$p
names(pvals) <- taxa
#results for decontam filtering with NP ordering
decontam_min[[i]] <- resSummary(X = Salter.counts, filtX = filtX, 
                          taxaInfo = taxaInfo, time = end)

}


## ------------------------------------------------------------------------
start <- Sys.time()
#traditional filtering
res_trad <- TraditR1(Counts, thresh =5)
end <-  Sys.time()-start
#results for traditional filtering with NP ordering
summary_trad_r1 <- resSummary(X = Counts, filtX = res_trad, 
                          taxaInfo = taxaInfo, time= end)


## ------------------------------------------------------------------------

#traditional filtering
start <- Sys.time()
res_trad <- TraditR2(Counts)
end <-  Sys.time()-start
#results for traditional filtering with NP ordering
summary_trad_r2 <- resSummary(X = Counts, filtX = res_trad, 
                          taxaInfo = taxaInfo,  time = end)


## ------------------------------------------------------------------------
#Manual fit to display all 4 quantile fits on one graph 
hist <- res_sim_sn_a$hist + stat_function(fun = dsn, aes(linetype = "Q1"),
      args = list(xi = res_sim_sn_a$est[1], omega = res_sim_sn_a$est[2], 
                  alpha = res_sim_sn_a$est[3]), colour="blue")+
      stat_function(fun = dsn, aes(linetype = "Q2"),
      args = list(xi = res_sim_sn_b$est[1], omega = res_sim_sn_b$est[2], 
                  alpha = res_sim_sn_b$est[3]), colour="blue")+ 
      stat_function(fun = dsn, aes(linetype = "Q3"),
      args = list(xi = res_sim_sn_c$est[1], omega = res_sim_sn_c$est[2], 
                  alpha = res_sim_sn_c$est[3]), colour="blue")+ 
      stat_function(fun = dsn, aes(linetype = "Q4"),
      args = list(xi = res_sim_sn_d$est[1], omega = res_sim_sn_d$est[2], 
                  alpha = res_sim_sn_d$est[3]), colour="blue")+
      scale_linetype_manual("", values = c("solid", "dashed", "dotted", "twodash"),
                            labels = c("5%, 10%, 25%", "10%, 25%, 40%", 
                                       "10%, 25%, 50%", "20%, 30%, 60%"))+
      theme(legend.position=c(0.6,0.8))
hist
#ggsave("~/Dropbox/PERFect/RCode/Responses/Plots/Quantiles_Salter.pdf")

## ------------------------------------------------------------------------
df <- rbind(summary_sim_sn_a_0.15$res,summary_sim_sn_b_0.15$res,
            summary_sim_sn_c_0.15$res,summary_sim_sn_d_0.15$res,
            summary_sim_sn_a$res,summary_sim_sn_b$res,
            summary_sim_sn_c$res,summary_sim_sn_d$res,
            summary_sim_sn_a_0.05$res,summary_sim_sn_b_0.05$res,
            summary_sim_sn_c_0.05$res,summary_sim_sn_d_0.05$res,
            summary_sim_sn_a_0.15_pvals$res,summary_sim_sn_b_0.15_pvals$res,
            summary_sim_sn_c_0.15_pvals$res,summary_sim_sn_d_0.15_pvals$res,
            summary_sim_sn_a_pvals$res,summary_sim_sn_b_pvals$res,
            summary_sim_sn_c_pvals$res,summary_sim_sn_d_pvals$res,
            summary_sim_sn_a_0.05_pvals$res,summary_sim_sn_b_0.05_pvals$res,
            summary_sim_sn_c_0.05_pvals$res,summary_sim_sn_d_0.05_pvals$res,
            summary_perm_np_a_0.15$res, summary_perm_np_b_0.15$res,
            summary_perm_np_c_0.15$res, summary_perm_np_d_0.15$res,
            summary_perm_np_a$res, summary_perm_np_b$res,
            summary_perm_np_c$res, summary_perm_np_d$res,
            summary_perm_np_a_0.05$res, summary_perm_np_b_0.05$res,
            summary_perm_np_c_0.05$res, summary_perm_np_d_0.05$res,
            summary_perm_pvals_a_0.15$res,summary_perm_pvals_b_0.15$res,
            summary_perm_pvals_c_0.15$res,summary_perm_pvals_d_0.15$res,
            summary_perm_pvals_a$res,summary_perm_pvals_b$res,
            summary_perm_pvals_c$res,summary_perm_pvals_d$res,
            summary_perm_pvals_a_0.05$res,summary_perm_pvals_b_0.05$res,
            summary_perm_pvals_c_0.05$res,summary_perm_pvals_d_0.05$res,
            summary_trad_r1$res, summary_trad_r2$res,
            decontam_pool[[1]]$res, decontam_pool[[2]]$res,
            decontam_pool[[3]]$res, decontam_pool[[4]]$res,
            decontam_min[[1]]$res, decontam_min[[2]]$res,
            decontam_min[[3]]$res, decontam_min[[4]]$res)
# rownames(df) <- c("PERFect_sim_a_0.15", "PERFect_sim_b_0.15", 
#                   "PERFect_sim_c_0.15", "PERFect_sim_d_0.15",
#                   "PERFect_sim_a_0.10", "PERFect_sim_b_0.10", 
#                   "PERFect_sim_c_0.10", "PERFect_sim_d_0.10",
#                   "PERFect_sim_a_0.05", "PERFect_sim_b_0.05", 
#                   "PERFect_sim_c_0.05", "PERFect_sim_d_0.05",
#                   "PERFect_sim_a_pvals_0.15", "PERFect_sim_b_pvals_0.15", 
#                   "PERFect_sim_c_pvals_0.15", "PERFect_sim_d_pvals_0.15",
#                   "PERFect_sim_a_pvals_0.10", "PERFect_sim_b_pvals_0.10", 
#                   "PERFect_sim_c_pvals_0.10", "PERFect_sim_d_pvals_0.10",
#                   "PERFect_sim_a_pvals_0.05", "PERFect_sim_b_pvals_0.05", 
#                   "PERFect_sim_c_pvals_0.05", "PERFect_sim_d_pvals_0.05",
#                   "PERFect_perm_np_a_0.15","PERFect_perm_np_b_0.15",
#                   "PERFect_perm_np_c_0.15","PERFect_perm_np_d_0.15",
#                   "PERFect_perm_np_a_0.10","PERFect_perm_np_b_0.10",
#                   "PERFect_perm_np_c_0.10","PERFect_perm_np_d_0.10",
#                   "PERFect_perm_np_a_0.05","PERFect_perm_np_b_0.05",
#                   "PERFect_perm_np_c_0.05","PERFect_perm_np_d_0.05",
#                   "PERFect_perm_pvals_a_0.15","PERFect_perm_pvals_b_0.15",
#                   "PERFect_perm_pvals_c_0.15","PERFect_perm_pvals_d_0.15",
#                   "PERFect_perm_pvals_a_0.10","PERFect_perm_pvals_b_0.10",
#                   "PERFect_perm_pvals_c_0.10","PERFect_perm_pvals_d_0.10",
#                   "PERFect_perm_pvals_a_0.05","PERFect_perm_pvals_b_0.05",
#                   "PERFect_perm_pvals_c_0.05","PERFect_perm_pvals_d_0.05",
#                   "Trad_R1", "Trad_R2",
#                   "decontam_pooled_0.05", "decontam_pooled_0.1",
#                   "decontam_pooled_0.2", "decontam_pooled_0.3",
#                   "decontam_batched_0.05", "decontam_batched_0.1",
#                   "decontam_batched_0.2", "decontam_batched_0.3")
df[,c(3,5)] <- round(df[,c(3,5)], 4)*100
df <- as.data.frame(df)

## ------------------------------------------------------------------------
df$pval <- c(rep(c(rep(0.15,4),rep(0.10,4),rep(0.05,4)),4),"Rule 1","Rule 2",rep(c(0.05,0.1,0.2,0.3),2))
df$setting <- c(rep(c("5%,10%,25%","10%,25%,40%","10%,25%,50%","20%,30%,60%"),12),rep(NA,10))
df$method <- c(rep(c("Simultaneous <br/> PERFect <br/> abundance <br/> ordering"),12),
               rep(c("Simultaneous <br/> PERFect <br/> p-values <br/> ordering"),12),
               rep(c("Permutation <br/> PERFect <br/> abundance <br/> ordering"),12),
               rep(c("Permutation <br/> PERFect <br/> p-values <br/> ordering"),12),
               rep(c("Traditional"),2),
               rep(c("Decontam prevalence"),8)
)
df <- df[,c("method","pval","setting","npres","pfilt","perccont")]
#alpha='\u03b1'
colnames(df) <- c("Method",paste0("Significance <br/> level <br/> alpha"),"Setting for <br> quantile <br> matching","# taxa <br/> preserved","% filtered","% contaminants <br/> preserved")
#kable(df,format="latex",booktabs=TRUE,caption="Sourcetracker results: estimated proportion of taxa from each environment in NICU samples (Knights et al, 2011)") %>% kable_styling(latex_options = c("striped", "hold_position"))
# names_spaced <- c(
#    'Total number<br/>of taxa', '# of taxa <br/> preserved', 
#   'Percentage of <br/> filtered taxa', 
#   '# of true taxa',
#   'Percentage of <br/> contaminant taxa')

kable(df,format='html', digits=2, align='lccccc',escape = FALSE)  %>% 
  kable_styling(full_width = F) %>%
  column_spec(1, bold = T) %>%
  collapse_rows(columns = 1:2, latex_hline = "major", valign = "top")

res_Salter <- df
#saveRDS(res_Salter, file = paste0(pathtodata, "test_Salter.RDS"))

## ------------------------------------------------------------------------
df <- data.frame(summary_sim_sn_a$time, summary_perm_np_a$time, 
                 summary_trad_r1$time, summary_trad_r2$time, decontam_pool[[1]]$time)
names(df) <- c("Simultaneous PERFect", "Permutation PERFect", 
               "Traditional Rule 1", "Traditional Rule 2", "decontam")

names_spaced <- c(
   'PERFect<br/> Simultaneous', 'PERFect<br/> Permutation', 
    'Traditional<br/> Rule 1', 'Traditional<br/> Rule 2','decontam')

kable(df, 
      format='html', 
      digits=4, 
      row.names=FALSE, 
      align='ccccc', 
      col.names = names_spaced,
      escape = FALSE) 



