#' A list containing mature and nascent counts for exons and introns, three replicates and 11 time points: 0,1/6,1/3,1/2,1,1.5,2,4,8,12,16 hours.
#'
#' @format A list of 4 matrices 500 x 44
#' @name allcounts
NULL

#' An INSPEcT object without modeled rates and concentrations
#' 
#' This INSPEcT object contains the rates first guesses for the entide database in allcounts
#'
#' @format An INSPEcT object
#' @name nascentInspObj
NULL

#' An INSPEcT object with evaluated and modeled rates and concentrations
#' 
#' This INSPEcT object contains the evaluated and modeled rates and concentrations
#' of the very first 10 genes of the dataset allcounts
#'
#' @format An INSPEcT object
#' @name nascentInspObj10
NULL

# #' An object containing four simulated datasets of 1000 simulated rates
# #'
# #' The INSPEcT_model object contains 1000 simulated rates that were obtained
# #' using the dataset rpkms as reference
# #'
# #'
# #' @format An INSPEcT_model object
# #' @name simRates
# NULL

# #' An INSPEcT object with 1000 simulated rates and concentration and their modeled rates
# #'
# #' A dataset containing the rates and concentrations obtained from the dataset
# #' simRates; 3 replicates and time points corresponding to: 0,1/6,1/3,1/2,1,1.5,2,4,8,12,16 hours.
# #'
# #' @format An INSPEcT object
# #' @name simData3rep_Nascent
# NULL


# #' An INSPEcT object with 1000 simulated rates and concentration and their modeled rates
# #'
# #' A dataset containing the rates and concentrations obtained from the dataset
# #' simRates; 3 replicates and time points corresponding to: 0,1/6,1/3,1/2,1,1.25,1.5,2,3,4,6,8,10,12,16 hours.
# #'
# #' @format An INSPEcT object
# #' @name simData4rep_Nascent
# NULL

# #' An INSPEcT object with 1000 simulated rates and concentration and their modeled rates
# #'
# #' A dataset containing the rates and concentrations obtained from the dataset
# #' simRates with 1 replicates and time points corresponding to: 0, 1/6, 1/3, 1/2, 
# #' 1, 2, 4, 8, 16 hours. On this dataset rates and concentrations
# #' have been modeled with the method modelRates
# #'
# #' @format An INSPEcT object
# #' @name simData3rep_NoNascent
# NULL

# #' An INSPEcT object with 1000 simulated rates and concentration and their modeled rates
# #'
# #' A dataset containing the rates and concentrations obtained from the dataset
# #' simRates with 1 replicates and time points corresponding to: 0, 1/6, 1/3, 1/2, 
# #' 1, 2, 4, 8, 16 hours. On this dataset rates and concentrations
# #' have been modeled with the method modelRates
# #'
# #' @format An INSPEcT object
# #' @name simData4rep_NoNascent
# NULL



