#mock data set 2
#@source \url{https://bmcmicrobiol.biomedcentral.com/articles/10.1186/s12866-015-0351-6}
load("mock2.RData")
#list named mock2 with taxa in columns and samples in rows
#Counts -- raw counts data
#Prop -- proportions data, counts divided by total sum of taxa observed in the sample
#TrueTaxa -- true taxa present in mock samples

#ravel vaginal microbiome data set
#@source \url{http://www.pnas.org/content/108/Supplement_1/4680.abstract?tab=ds}
load("ravel.RData")
#list named ravel with taxa in columns and samples in rows
#Counts -- raw counts data
#Prop -- proportions data, counts divided by total sum of taxa observed in the sample
#Sample_Data