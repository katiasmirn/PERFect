#' Ravel vaginal microbiome data
#'
#' The vaginal microbiome study of a cohort of 396 reproductive age women previously published in (Ravel et al., 2011),
#'  where data and details on data collection and preprocessing are also available. The goal of the study was to understand
#'  the role and ultimate function of vaginal microbiota in reducing risk of infectious diseases and to identify factors
#'  leading to disease susceptibility. Microbiome data were obtained by pyrosequencing of barcoded 16S rRNA genes;
#'  in this data set 247 taxa were identified.
#'
#' @docType data
#' @keywords datasets
#' @name ravel
#' @usage data(ravel)
#' @format This file contains the OTU tables in counts and proportion, as well as its metadata. The OTU table has 394 samples
#'  and 247 taxa.
"ravel"

#' Bias experiment data
#'
#' These publicly available data (Brooks et al., 2015) were generated as a part of a study designed to evaluate the bias
#'  at each step of the VCU sequencing protocol, namely,  DNA extraction,  PCR amplification,  sequencing and taxonomic
#'  classification. Mock community samples were created out of 7 vaginally relevant bacteria  by mixing prescribed
#'  quantities of cells, with quantities varying across samples according to an experimental design described in
#'  Brooks et al, 2015. As opposed to the positive controls data, bacteria appear in different proportions across samples.
#'  The number of taxa identified by the sequencing and bioinformatics pipeline was 46.
#'
#' @docType data
#' @keywords datasets
#' @name mock2
#' @usage data(mock2)
#' @format This file contains a count OTU table and a proportion OTU table, each with 240 samples and 46 taxa. A list of true
#'  taxa is also given.
NULL


#' Bacterial Diversity in Neonatal Intensive Care Units (NICUs) data
#'
#' These data (Knights et al., 2011) with 30 samples and 1097 taxa was collected to investigat the sources of bacteria
#'  found on surfaces and equipment in NICU. The data taxa was previously analyzed using sourcetracker software to identify
#'  the proportion of bacteria from each environment using published datasets from environments likely to be sources of indoor
#'  contaminants, namely human skin, oral cavities, feces and temperate soils.
#'
#' @docType data
#' @keywords datasets
#' @name knight
#' @usage data(knight)
#' @format This file contains....
"knight"

#' Reagent and laboratory contamination data
#'
#' The data (Salter et al., 2014) was generated from the study of the effect of present contaminants in DNA extraction kits
#'  and other laboratory reagents on sequencing DNA. Mock samples of a pure Salmonella bongori culture had undergone five
#'  rounds of serial ten-fold dilutions to generate a series of high to low biomass samples. To generate a taxa counts table
#'  from this study, we used samples for the Salmonella bongori culture 16S rRNA gene profiling data, which are deposited as
#'  FASTQ files under ENA project accession EMBL: ERP006737 (https://www.ebi.ac.uk/ena/data/view/PRJEB7055),
#'  and processed using the dada2 R-package.
#'
#' @docType data
#' @keywords datasets
#' @name salter
#' @usage data(salter)
#' @format This file contains the OTU table, its metadata and taxomomy. The OTU table has 42 samples and 625 taxa.
NULL

#' Oral artificial communities data
#'
#' This dataset comes from the MBQC project, a collaborative effort designed to comprehensively evaluate sample processing
#'  and computational methods for human microbiome data analysis (Sinha et al., 2015). There are four types of samples in
#'  this project: (1) 11 unique fresh stool samples; (2) seven unique freeze-dried stool samples; (3) two unique chemostat
#'  samples generated from a Robogut; and (4) two artificial colonies representing the gut and oral cavity. These samples
#'  were randomly sequenced at 15 laboratories and then randomly distributed to 9 bioinformatics facilities for microbiome
#'  analyses. Here, we consider the oral artificial communities data which comprised of 22 true taxa from the human
#'  oral cavity.
#'
#' @docType data
#' @keywords datasets
#' @name mock
#' @usage data(mock)
#' @format This phyloseq object contains an OTU table with 1016 samples and 27211 taxa, its metadata and taxomomy.
#'
"mock"






