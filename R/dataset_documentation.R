#' Example daily meteorology data
#' 
#' Example data set of meteorological input.
#' 
#' @name examplemeteo
#' @docType data
#' 
#' @format
#' A data frame containing daily meteorology of a location in Catalonia (Spain) for year 2001:
#' \describe{
#'   \item{\code{dates}}{Vector of \code{\link{Date}} objects.}
#'    \item{\code{MinTemperature}}{Minimum daily temperature (in degrees Celsius).}
#'    \item{\code{MaxTemperature}}{Maximum daily temperature (in degrees Celsius).}
#'    \item{\code{Precipitation}}{Daily precipitation (in mm of water).}
#'    \item{\code{MinRelativeHumidity}}{Minimum daily relative humidity (in percent).}
#'    \item{\code{MaxRelativeHumidity}}{Maximum daily relative humidity (in percent).}
#'    \item{\code{Radiation}}{Incoming radiation (in MJ/m2).}
#'    \item{\code{WindSpeed}}{Wind speed (in m/s).}
#' }
#' 
#' @source Interpolated from weather station data (Spanish and Catalan meteorology agencies) using package 'meteoland'.
#' 
#' @seealso \code{\link{spwb}}
#' @examples
#'  data(examplemeteo)
#' @keywords data
NULL

#' Data tables with species parameter definition and values for different countries
#' 
#' A data sets of species parameter definition and values, the latter resulting from existing databases, fit to empirical data or expert-based guesses.
#' 
#' @name SpParams
#' @aliases SpParamsDefinition SpParamsMED SpParamsES SpParamsFR SpParamsUS
#' 
#' @docType data 
#' 
#' @format
#' \itemize{
#'   \item{Data frame \code{SpParamsDefinition} has parameters in rows and columns 'ParameterName', 'ParameterGroup', 'Definition', 'Type' and 'Units'.}
#'   \item{Data frames \code{SpParamsMED} (for Catalonia), \code{SpParamsES} (for Spain), \code{SpParamsFR} (for France) and \code{SpParamsUS} (for US) have species or genus as rows and column names equal to parameter names in \code{SpParamsDefinition}.}
#' }
#' @details
#' Plant trait parameter sources are listed in the bibliography section. Details of the procedures used to obtain the species parameter tables can be found in an article at https://emf-creaf.github.io/medfate/. 
#' @examples
#' data(SpParamsDefinition)
#' data(SpParamsMED)
#' @keywords data
NULL