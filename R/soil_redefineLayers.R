#' Redefine soil layer widths
#'
#' Allows redefining soil layer widths of an input data frame of soil parameters.
#'
#' @param x A data frame of soil parameters (see an example in \code{\link{defaultSoilParams}}) or an object of class \code{\link{soil}}.
#' @param widths A numeric vector indicating the desired layer widths, in mm.
#'
#' @author \enc{Víctor}{Victor} Granda, EMF-CREAF
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, EMF-CREAF
#'
#' @details
#' If an initialized \code{\link{soil}} is supplied, its hydraulic parameters will be recalculated and the value of state variables will be lost.
#' 
#' @return A data frame or \code{\link{soil}} object with soil parameters, depending on the class of \code{x}. 
#' 
#' @seealso  \code{\link{soil}}, \code{\link{defaultSoilParams}}
#' @export
#' 
#' @examples
#' # Define initial soil with 5 layers
#' spar <- defaultSoilParams(5)
#' spar
#' 
#' # Redefine to four layers
#' soil_redefineLayers(spar)
#' 
#' # Same but after soil parameter initialization
#' examplesoil <- soil(spar)
#' examplesoil
#' 
#' soil_redefineLayers(examplesoil)
#'
soil_redefineLayers<-function(x, widths = c(300, 700, 1000, 2000)) {
  
  is_soil <- inherits(x, "soil")
  
  restarget  <-  data.frame(matrix(nrow = length(widths), ncol = 7))
  names(restarget)  <-  c("widths", "clay", "sand", "om", "nitrogen", "bd", "rfc")
  restarget$widths  <-  widths
  
  bottomdepths <- cumsum(widths)
  topdepths <- bottomdepths - widths
  sgbdepths <- cumsum(x$widths)
  sgtdepths <- sgbdepths - x$widths
  for(j in 1:length(widths)) {
    ini <- topdepths[j]
    fin <- bottomdepths[j]
    p1 <- pmin(pmax(fin -sgtdepths,0),x$widths)
    p2 <- pmin(pmax(sgbdepths -ini,0),x$widths)
    p <- pmin(p1,p2)/x$widths
    if(sum(p)==0) p[length(p)]  <-  1
    restarget$clay[j]  <-  sum(x$clay*p)/sum(p)
    restarget$sand[j]  <-  sum(x$sand*p)/sum(p)
    restarget$rfc[j]  <-  sum(x$rfc*p)/sum(p)
    restarget$bd[j]  <-  sum(x$bd*p)/sum(p)
    restarget$om[j]  <-  sum(x$om*p)/sum(p)
    restarget$nitrogen[j]  <-  sum(x$nitrogen*p)/sum(p)
  }
  if(!is_soil) return(restarget)
  return(soil(restarget))
}
