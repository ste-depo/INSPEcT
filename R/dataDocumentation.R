#' A list containing mature and nascent counts for exons and introns, three replicates and 11 time points: 
#' 0,1/6,1/3,1/2,1,1.5,2,4,8,12,16 hours.
#'
#' @format A list of 4 matrices 500 x 33
#' @name allcounts
NULL

#' Contains two variables: "exWdths" and "intWdths" containing the lenght of the exons and introns, 
#' respectively, relative to the genes in "allcounts"
#'
#' @format numeric vector of length 500 
#' @name featureWidths
NULL

#' Contains two variables: "nascentLS" and "totalLS" containing the sequencing depth of nascent and total libraries 
#' respectively, relative to the experiments in "allcounts"
#'
#' @format numeric vector of length 33 
#' @name libsizes
NULL

#' An INSPEcT object with 1000 simulated rates and concentration and their modeled rates
#'
#' A dataset containing the rates and concentrations obtained from the dataset
#' simRates; 3 replicates and time points corresponding to: 0,1/6,1/3,1/2,1,1.5,2,4,8,12,16 hours.
#'
#' @format An INSPEcT object
#' @name simData3rep_Nascent
NULL


#' An INSPEcT object with 1000 simulated rates and concentration and their modeled rates
#'
#' A dataset containing the rates and concentrations obtained from the dataset
#' simRates; 3 replicates and time points corresponding to: 0,1/6,1/3,1/2,1,1.25,1.5,2,3,4,6,8,10,12,16 hours.
#'
#' @format An INSPEcT object
#' @name simData4rep_Nascent
NULL

#' An INSPEcT object with 1000 simulated rates and concentration and their modeled rates
#'
#' A dataset containing the rates and concentrations obtained from the dataset
#' simRates with 1 replicates and time points corresponding to: 0, 1/6, 1/3, 1/2, 
#' 1, 2, 4, 8, 16 hours. On this dataset rates and concentrations
#' have been modeled with the method modelRates
#'
#' @format An INSPEcT object
#' @name simData3rep_NoNascent
NULL

#' An INSPEcT object with 1000 simulated rates and concentration and their modeled rates
#'
#' A dataset containing the rates and concentrations obtained from the dataset
#' simRates with 1 replicates and time points corresponding to: 0, 1/6, 1/3, 1/2, 
#' 1, 2, 4, 8, 16 hours. On this dataset rates and concentrations
#' have been modeled with the method modelRates
#'
#' @format An INSPEcT object
#' @name simData4rep_NoNascent
NULL