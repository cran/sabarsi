#' A real SERS spectrum data set in
#'
#'This data contains the SERS spectra from two technical replicates.
#'
#' @details This data contains the SERS spectra from two technical replicates. Each replicate contains spectra with 500 frequency channels from 500 time points.
#'
#' @docType data
#'
#' @usage data(SERS)
#'
#' @format An object containing the following variables:
#' \describe{
#' \item{\code{R1}}{A data matrix of SERS spectra in replicate 1.}
#' \item{\code{R2}}{A data matrix of SERS spectra in replicate 2.}
#'
#' }
#'
#' @references Wang et. al "A Statistical Approach of Background Removal and Spectrum Identification for SERS Data"
#' @keywords datasets
#'
#' @examples
#' data(SERS)
#' x <- list()
#' x[[1]] <- SERS$R1
#' x[[2]] <- SERS$R2
"SERS"
