#' intron, exons RPKM values of 500 genes from 4sU-seq and RNA-seq experiments
#'
#' A dataset containing the values of exonic and intronic RPKMs of 500 genes 
#' both in 4sU-seq and RNA-seq experiments. The variables are as follows:
#'
#' \itemize{
#'   \item rpkms_4su_exons
#'   \item rpkms_4su_introns
#'   \item rpkms_total_exons
#'   \item rpkms_total_introns
#' }
#'
#' @format A list of 4 matrices with 500 rows and 9 columns
#' @name rpkms
NULL

#' An INSPEcT object with evaluated and modeled rates and concentrations
#' 
#' This INSPEcT object contains the evaluated and modeled rates and concentrations
#' of the very first 10 genes of the dataset rpkms
#'
#' @format An INSPEcT object
#' @name mycerIds10
NULL

#' An INSPEcT_model object of 1000 simulated rates
#'
#' The INSPEcT_model object contains 1000 simulated rates that were obtained
#' using the dataset rpkms as reference
#'
#'
#' @format An INSPEcT_model object
#' @name simRates
NULL

#' An INSPEcT object with 1000 simulated rates and concentration and their modeled rates
#'
#' A dataset containing the rates and concentrations obtained from the dataset
#' simRates with 1 replicates and time points corresponding to: 0, 1/6, 1/3, 1/2, 
#' 1, 2, 4, 8, 16 hours. On this dataset rates and concentrations
#' have been modeled with the method modelRates
#'
#' @format An INSPEcT object
#' @name simData1rep
NULL

#' An INSPEcT object with 1000 simulated rates and concentration and their modeled rates
#'
#' A dataset containing the rates and concentrations obtained from the dataset
#' simRates with 3 replicates and time points corresponding to: 0, 1/6, 1/3, 1/2, 
#' 1, 1.5, 2, 4, 8, 12, 16 and 24 hours. On this dataset rates and concentrations
#' have been modeled with the method modelRates
#'
#' @format An INSPEcT object
#' @name simData3rep
NULL


#' A Gtf file containing exons definition of 100 genes
#'
#' @format A tab separated file
#' @name UCSC_mm9_genes_exons
NULL

#' A Gtf file containing introns definition of 100 genes
#'
#' @format A tab separated file
#' @name UCSC_mm9_genes_introns
NULL
