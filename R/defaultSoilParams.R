#' Default soil parameters
#'
#' Creates a data frame with default soil physical description for model functions
#' 
#' @param n An integer with the number of soil layers (between two and five).
#' 
#' @details The function returns a data frame with default physical soil description, with soil layers in rows. 
#' Users can change those that need to be set to other values and use the list as input for function \code{\link{soil}}.
#' 
#' @note While this function is limited to five soil layers, user defined data frames can discretize soils using an unlimited number of soil layers.
#' 
#' @return A data frame with layers in rows and the following columns (and default values):
#' \itemize{
#'   \item{\code{widths (= c(300,700,1000,2000)}: Width of soil layers (in mm).}
#'   \item{\code{clay (= 25)}: Clay percentage for each layer (in \%).}
#'   \item{\code{sand (= 25)}: Sand percentage for each layer (in \%).}
#'   \item{\code{om (= NA)}: Organic matter percentage for each layer (in \%) (optional).}
#'   \item{\code{nitrogen (= NA)}: Sum of total nitrogen (ammonia, organic and reduced nitrogen) for each layer (in g/kg) (optional).}
#'   \item{\code{bd (= 1.5)}: Bulk density for each layer (in g/cm3).}
#'   \item{\code{rfc (= c(20,40,60,85))}: Percentage of rock fragment content (volume basis) for each layer.}
#' }
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso  \code{\link{soil}}, \code{\link{soil_redefineLayers}},  \code{\link{defaultControl}}, \code{\link{SpParamsMED}}
#' 
#' @examples 
#' defaultSoilParams(4)
#' 
defaultSoilParams<-function(n=4) {
  n = min(max(n,2),5)
  return(data.frame(
    widths = c(300,700,1000,2000,4000)[1:n],
    clay = rep(25,n),
    sand = rep(25,n),
    om = rep(NA,n),
    nitrogen = rep(NA,n),
    bd = rep(1.5,n),
    rfc = c(25,45,75,95,98)[1:n]));
}