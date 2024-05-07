#' Redefine soil layer widths
#'
#' Allows redefining soil layer widths of an input data frame of soil parameters.
#'
#' @param SoilParams A data frame of soil parameters (see an example in \code{\link{defaultSoilParams}}).
#' @param widths A numeric vector indicating the desired layer widths, in mm.
#'
#' @author \enc{Víctor}{Victor} Granda, EMF-CREAF
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, EMF-CREAF
#'
#' @return A modified data frame of soil parameters
#'
#' @seealso  \code{\link{soil}}, \code{\link{defaultSoilParams}}
#' @export
#'
soil_redefineLayers<-function(SoilParams, widths = c(300, 700, 1000, 2000)) {
  restarget = data.frame(matrix(nrow = length(widths), ncol = 7))
  names(restarget) = c("widths", "clay", "sand", "om", "nitrogen", "bd", "rfc")
  restarget$widths = widths
  
  bottomdepths = cumsum(widths)
  topdepths = bottomdepths - widths
  sgbdepths = cumsum(SoilParams$widths)
  sgtdepths = sgbdepths - SoilParams$widths
  for(j in 1:length(widths)) {
    ini = topdepths[j]
    fin = bottomdepths[j]
    p1 = pmin(pmax(fin -sgtdepths,0),SoilParams$widths)
    p2 = pmin(pmax(sgbdepths -ini,0),SoilParams$widths)
    p = pmin(p1,p2)/SoilParams$widths
    if(sum(p)==0) p[length(p)] = 1
    restarget$clay[j] = sum(SoilParams$clay*p)/sum(p)
    restarget$sand[j] = sum(SoilParams$sand*p)/sum(p)
    restarget$rfc[j] = sum(SoilParams$rfc*p)/sum(p)
    restarget$bd[j] = sum(SoilParams$bd*p)/sum(p)
    restarget$om[j] = sum(SoilParams$om*p)/sum(p)
    restarget$nitrogen[j] = sum(SoilParams$nitrogen*p)/sum(p)
  }
  return(restarget)
}
