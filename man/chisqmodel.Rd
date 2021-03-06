% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/chisqtest-methods.R
\name{chisqmodel}
\alias{chisqmodel}
\alias{chisqmodel,INSPEcT-method}
\title{Retrieve results of chi-squared test for the selected models}
\usage{
chisqmodel(object, gc = NULL, tpts = NULL, ...)

\S4method{chisqmodel}{INSPEcT}(object, gc = NULL, tpts = NULL, ...)
}
\arguments{
\item{object}{An object of class INSPEcT or INSPEcT_model}

\item{gc}{Additional arguments for the generic}

\item{tpts}{Additional arguments for the generic}

\item{...}{Additional arguments for the generic}
}
\value{
A vector of chi-squared test results
}
\description{
This method is used to retrieve the chi-squared test results for the models
that have been selected to better represent the behavior of each gene.
}
\examples{
nascentInspObj10 <- readRDS(system.file(package='INSPEcT', 'nascentInspObj10.rds'))
chisqmodel(nascentInspObj10)
}
