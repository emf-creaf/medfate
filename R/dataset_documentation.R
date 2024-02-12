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