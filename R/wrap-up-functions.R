#' Wrapper function from BAM files
#'
#' @description
#' Function to run the whole INSPEcT differential rate analysis procedure with a single line. 
#' The function save the output analysis to file that can be later loaded in the R environment 
#' or in the INSPEcT-GUI.
#' @param txdb A TranscriptDB object for the selected organism
#' @param annotation_table Paths and experimental design associated to bam files. They could be
#' provided directly as a 'data.frame', or as a path to the file containing the information. Possible file
#' formats are csv' (comma-separated-values), 'tsv' (comma-separated-values), or 'xls' (Excel).
#' In case 'annotation_table' has 2 colums named 'condition' and 'total', INSPEcT- analysis is run.
#' In case 'annotation_table' has 3 colums named 'condition', 'total' and 'nascent', INSPEcT+ analysis is run.
#' 'condition' is a colums indicating the experimental condition, a character vector (containing, for example, 'WT'
#' or 'KD') in case of steady-state experiments, or numerical values indicating the time from the unperturbed
#' condition in case of time-course analysis. 'total' and 'nascent' contains the path to totalRNA and nascentRNA
#' BAM files, respectively.
#' @param labeling_time A numeric indicating the time of labeling exposure to the modified nucleotide.
#' To be indicated only in case of INSPEcT+ analysis.
#' @param strandSpecific A numeric indicating the strandness of the BAM files, 0 for non strand-specific,
#' 1 for stranded, 2 for reversely-stranded. 0 by default.
#' @param isPairedEnd A logical indicating if paired-end sequencing have been performed. FALSE by default.
#' @param estimateRatesWith Either "int" or "der". With "int" the degradation and processing
#'    rates are estimated integrating the system between one time point and the following. 
#'    With "der" degradation and processing rates are estimated using the derivative of total
#'    and pre mRNA. (default is "der")
#' @param useSigmoidFun A logical, whether to choose between sigmoid and impulse function 
#'    to fit rates and concentrations. In case not, always impulse function is used. 
#'    (default is TRUE)
#' @param file A character indicating where the output of the analysis will be stored. If not provided 
#' the file name will be created automaticcally and saved on the current folder.
inspectFromBAM <- function(txdb, annotation_table, labeling_time=NULL, 
													 strandSpecific=0, isPairedEnd=FALSE, 
													 estimateRatesWith = 'der', useSigmoidFun = TRUE,
													 file=NULL) {
	
	#########################
	## check txdb object ####
	#########################
	
	if( class(txdb) != 'TxDb' ) 
		stop('inspectFromBAM: "txdb" must be an object of TxDb class.')
	
	################################
	## Read the annotation file ####
	################################
	
	if( class(annotation_table) != "data.frame" ) {
		if( class(annotation_table) == "character" ) {
			if(!file.exists(annotation_table)) {
				stop('annotation_table is nor a data.frame or a path to an existing file')
			} else {
				file_extension <- tail(strsplit(annotation_table, split='\\.')[[1]],1)
				annotation_table <- switch(file_extension
							 , 'csv'=read.csv(annotation_table, header=TRUE, stringsAsFactors = FALSE)
							 , 'tsv'=read.table(annotation_table, sep='\t', header=TRUE, stringsAsFactors = FALSE)
							 , 'xls'=read.xls(annotation_table, stringsAsFactors = FALSE)
				)
			}
		} else {
			stop('annotation_table must be either a data.frame or a character containing the path to the annotation file')
		}
	}
	
	if( all(c('condition','total','nascent') %in% colnames(annotation_table)) ) {
		analysis_with_nascent <- TRUE
		if( is.null(labeling_time) ) {
			stop('In case of INSPEcT+ analysis, labeling_time must be provided')
		} else {
			if( !is.numeric(labeling_time) ) stop('labeling_time must be numeric')
		}
	} else if( all(c('condition','total') %in% colnames(annotation_table)) ) {
		analysis_with_nascent <- FALSE
	} else {
		stop(paste(
			'annotation_table columns (', 
			paste(colnames(annotation_table), collapse = ', '),
			') have not recognized names'))
	}
	
	tpts<-sort(unique(annotation_table$condition))
	totalExpressions <- quantifyExpressionsFromBAMs(txdb = txdb, 
															BAMfiles = as.character(annotation_table$total), 
															experimentalDesign = annotation_table$condition,
															strandSpecific=strandSpecific, isPairedEnd=isPairedEnd
															)
	
	if( analysis_with_nascent ) {
		nascentEpressions <- quantifyExpressionsFromBAMs(txdb = txdb, 
																BAMfiles = as.character(annotation_table$nascent), 
																experimentalDesign = annotation_table$condition,
																strandSpecific=strandSpecific, isPairedEnd=isPairedEnd
																)
		if(!is.numeric(tpts)) {
			nascentInspObj<-newINSPEcT(tpts = as.character(tpts),
																 labeling_time = labeling_time,
																 nascentExpressions = nascentEpressions,
																 matureExpressions = totalExpressions)
			message('Created steady state object.')
		} else {
			nascentInspObj<-newINSPEcT(tpts = tpts,
																 labeling_time = labeling_time,
																 nascentExpressions = nascentEpressions,
																 matureExpressions = totalExpressions)
			nascentInspObj<-modelRates(nascentInspObj, estimateRatesWith = estimateRatesWith, useSigmoidFun = useSigmoidFun)
		}
		if( is.null(file) ) {
			file <- paste0('nascentInspectAnalysis_tL_', signif(labeling_time,2), 'h.rds')
		}
		saveRDS(nascentInspObj, file = file)
		message('nascent INSPEcT analysis completed. File stored at')
		message(file)
	} else {
		if(!is.numeric(tpts)) {
			totalInspObj<-newINSPEcT(tpts=as.character(tpts), matureExpressions=totalExpressions)
			message('Created steady state object.')
		} else {
			totalInspObj<-newINSPEcT(tpts, matureExpressions=totalExpressions)
			totalInspObj<-modelRates(totalInspObj, estimateRatesWith = estimateRatesWith, useSigmoidFun = useSigmoidFun)	
		}
		if( is.null(file) ) {
			file <- paste0('totalInspectAnalysis.rds')
		}
		saveRDS(totalInspObj, file = file)
		message('total INSPEcT analysis completed. File stored at')
		message(file)
	}

}

#' Wrapper function from PCR quantifications
#' 
#' @description
#' Function to run the whole INSPEcT differential rate analysis procedure with a single line. 
#' The function save the output analysis to file that can be later loaded in the R environment 
#' or in the INSPEcT-GUI.
#' @param totalRNA_table Exonic quantification, intronic quantification and experimental design 
#' associated to totalRNA of a single gene quantified by PCR. They could be provided directly as a 
#' 'data.frame', or as a path to the file containing the information. Possible file
#' formats are csv' (comma-separated-values), 'tsv' (comma-separated-values), or 'xls' (Excel).
#' 'totalRNA_table' must have 3 colums named 'condition', 'total_exonic' and 'total_intronic'.
#' 'condition' is a column indicating the experimental condition, a character vector (containing, for example, 'WT'
#' or 'KD') in case of steady-state experiments, or numerical values indicating the time from the unperturbed
#' condition in case of time-course analysis. 'total_exonic' and 'total_intronic' contains abundance of gene
#' measured in its exonic and intronic regions, respectively, in the total RNA fraction.
#' @param nascentRNA_table similar to 'totalRNA_table' but referred to nascent RNA fraction. In this case, 
#' colums names must be 'condition', 'nascent_exonic' and 'nascent_intronic'. In case this infromation 
#' is not provided, INSPEcT- analysis is run. If otherwise this information is present, INSPEcT+ analysis is run.
#' @param labeling_time A numeric indicating the time of labeling exposure to the modified nucleotide.
#' To be indicated only in case of INSPEcT+ analysis.
#' @param estimateRatesWith Either "int" or "der". With "int" the degradation and processing
#'    rates are estimated integrating the system between one time point and the following. 
#'    With "der" degradation and processing rates are estimated using the derivative of total
#'    and pre mRNA. (default is "der")
#' @param useSigmoidFun A logical, whether to choose between sigmoid and impulse function 
#'    to fit rates and concentrations. In case not, always impulse function is used. 
#'    (default is TRUE)
#' @param file A character indicating where the output of the analysis will be stored. If not provided 
#' the file name will be created automaticcally and saved on the current folder.
#' @examples
#' if( Sys.info()["sysname"] != "Windows" ) {
#' totalAnnTabPCR <- system.file(package = 'INSPEcT', 'totalAnnTabPCR.csv')
#' nascentAnnTabPCR <- system.file(package = 'INSPEcT', 'nascentAnnTabPCR.csv')
#' inspectFromPCR(totalAnnTabPCR, nascentAnnTabPCR, labeling_time=1/6)
#' }
inspectFromPCR <- function(totalRNA_table, nascentRNA_table=NULL, labeling_time=NULL, 
													 estimateRatesWith = 'der', useSigmoidFun = TRUE,
													 file=NULL) {
	
	####################################
	## Read the totalRNA_table file ####
	####################################
	
	if( class(totalRNA_table) != "data.frame" ) {
		
		if( class(totalRNA_table) == "character" ) {
			if(!file.exists(totalRNA_table)) {
				stop('totalRNA_table is nor a data.frame or a path to an existing file')
			} else {
				file_extension <- tail(strsplit(totalRNA_table, split='\\.')[[1]],1)
				totalRNA_table <- switch(file_extension
																	 , 'csv'=read.csv(totalRNA_table, header=TRUE, stringsAsFactors = FALSE)
																	 , 'tsv'=read.table(totalRNA_table, sep='\t', header=TRUE, stringsAsFactors = FALSE)
																	 , 'xls'=read.xls(totalRNA_table, stringsAsFactors = FALSE)
				)
			}
		} else {
			stop('totalRNA_table must be either a data.frame or a character containing the path to the annotation file')
		}
		
	}
	
	if( !all(c('condition','total_exonic','total_intronic') %in% colnames(totalRNA_table)) ) {
		stop(paste(
			'totalRNA_table columns (', 
			paste(colnames(totalRNA_table), collapse = ', '),
			') have not recognized names'))
	}
	
	##################################
	## Read nascentRNA_table file ####
	##################################
	
	if( !is.null(nascentRNA_table) ) {
		
		if( class(nascentRNA_table) != "data.frame" ) {
			if( class(nascentRNA_table) == "character" ) {
				if(!file.exists(nascentRNA_table)) {
					stop('nascentRNA_table is nor a data.frame or a path to an existing file')
				} else {
					file_extension <- tail(strsplit(nascentRNA_table, split='\\.')[[1]],1)
					nascentRNA_table <- switch(file_extension
																		 , 'csv'=read.csv(nascentRNA_table, header=TRUE, stringsAsFactors = FALSE)
																		 , 'tsv'=read.table(nascentRNA_table, sep='\t', header=TRUE, stringsAsFactors = FALSE)
																		 , 'xls'=read.xls(nascentRNA_table, stringsAsFactors = FALSE)
					)
				}
			} else {
				stop('nascentRNA_table must be either a data.frame or a character containing the path to the annotation file')
			}
		} 
		
		analysis_with_nascent <- TRUE
		if( is.null(labeling_time) ) {
			stop('In case of INSPEcT+ analysis, labeling_time must be provided')
		} else {
			if( !is.numeric(labeling_time) ) stop('labeling_time must be numeric')
		}
		
		if( !all(c('condition','nascent_exonic','nascent_intronic') %in% colnames(nascentRNA_table)) ) {
			stop(paste(
				'nascentRNA_table columns (', 
				paste(colnames(nascentRNA_table), collapse = ', '),
				') have not recognized names'))
		}	
		
	} else {
		analysis_with_nascent <- FALSE
	}
	
	tpts<-sort(unique(totalRNA_table$condition))
	totalExpressions <- list(
		exonsExpressions = t(tapply(totalRNA_table$total_exonic, totalRNA_table$condition, mean)),
		intronsExpressions = t(tapply(totalRNA_table$total_intronic, totalRNA_table$condition, mean)),
		exonsVariance = t(tapply(totalRNA_table$total_exonic, totalRNA_table$condition, var)),
		intronsVariance = t(tapply(totalRNA_table$total_intronic, totalRNA_table$condition, var))
	)
	rownames(totalExpressions$exonsExpressions) <- rownames(totalExpressions$intronsExpressions) <- 
		rownames(totalExpressions$exonsVariance) <- rownames(totalExpressions$intronsVariance) <- '1'
	
	if( analysis_with_nascent ) {
		tpts_nascent<-sort(unique(nascentRNA_table$condition))
		if(!identical(tpts_nascent, tpts)) {
			stop('all conditions in totalRNA_table must also be in nascentRNA_table')
		}
		nascentEpressions <- list(
			exonsExpressions = t(tapply(nascentRNA_table$nascent_exonic, nascentRNA_table$condition, mean)),
			intronsExpressions = t(tapply(nascentRNA_table$nascent_intronic, nascentRNA_table$condition, mean)),
			exonsVariance = t(tapply(nascentRNA_table$nascent_exonic, nascentRNA_table$condition, var)),
			intronsVariance = t(tapply(nascentRNA_table$nascent_intronic, nascentRNA_table$condition, var))
		)
		rownames(nascentEpressions$exonsExpressions) <- rownames(nascentEpressions$intronsExpressions) <- 
			rownames(nascentEpressions$exonsVariance) <- rownames(nascentEpressions$intronsVariance) <- '1'
		# increase the parameters for minimization (it's only one gene!)
		if(!is.numeric(tpts)) {
			nascentInspObj<-newINSPEcT(tpts = as.character(tpts),
																 labeling_time = labeling_time,
																 nascentExpressions = nascentEpressions,
																 matureExpressions = totalExpressions)
			message('Created steady state object.')
		} else {
			nascentInspObj<-newINSPEcT(tpts = tpts,
																 labeling_time = labeling_time,
																 nascentExpressions = nascentEpressions,
																 matureExpressions = totalExpressions)
			nascentInspObj<-modelRates(nascentInspObj, estimateRatesWith = estimateRatesWith, useSigmoidFun = useSigmoidFun, nInit = 50)	
		}
		if( is.null(file) ) {
			file <- paste0('nascentInspectAnalysisFromPCR_tL_', signif(labeling_time,2), 'h.rds')
		}
		saveRDS(nascentInspObj, file = file)
		message('nascent INSPEcT analysis completed. File stored at')
		message(file)
	} else {
		if(!is.numeric(tpts)) {
			totalInspObj<-newINSPEcT(tpts=as.character(tpts), matureExpressions=totalExpressions)
			message('Created steady state object.')
		} else {
			# increase the parameters for minimization (it's only one gene!)
			totalInspObj<-newINSPEcT(tpts=tpts, matureExpressions=totalExpressions)
			totalInspObj<-modelRates(totalInspObj, estimateRatesWith = estimateRatesWith, useSigmoidFun = useSigmoidFun, nInit = 50)
		}
		if( is.null(file) ) {
			file <- paste0('totalInspectAnalysisFromPCR.rds')
		}
		saveRDS(totalInspObj, file = file)
		message('total INSPEcT analysis completed. File stored at')
		message(file)
	}
	
}



