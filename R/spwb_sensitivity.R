#' Sensitivity analysis for soil plant water balance simulations
#' 
#' Performs a set of calls to \code{\link{spwb}} with the aim to determine the sensitivity to particular parameters.
#' 
#' @param x An object of class \code{\link{spwbInput}}.
#' @param soil A list containing the description of the soil (see \code{\link{soil}}).
#' @param meteo A data frame with daily meteorological data series (see \code{\link{spwb}}).
#' @param paramType Data frame of \code{x} to modify.
#' @param paramName Name of the parameter to modify.
#' @param cohort Integer with the cohort to modify.
#' @param p_change Numerical vector with percentages of change.
#' @param summary.fun Summary function to be applied to the results of each simulation. 
#' @param simplify Whether the result of \code{summary.fun} should be simplified (see \code{\link{sapply}}). 
#' @param ... Additional parameters to function \code{\link{spwb}}.
#' 
#' @details Due to parameter dependence, modifying some parameters affects others:
#' \itemize{
#'   \item{Setting \code{paramName = "Z50/Z95"} affects \code{belowLayers$V}, \code{belowLayers$VCroot_kmax} and \code{belowLayers$VGrhizo_kmax}.}
#'   \item{Modifying \code{LAI_live} also affects \code{LAI_expanded}.}
#'   \item{Modifying \code{VCroot_kmax} from \code{paramsTranspiration} affects both \code{VCroot_kmax} and \code{belowLayers$VCroot_kmax}.}
#'   \item{Modifying \code{WaterStorage} affects simultaneously \code{Vleaf} and \code{Vsapwood} from \code{paramsWaterStorage}.}
#'   \item{Modifying \code{c} from \code{paramsTranspiration} affects simultaneously \code{VCleaf_c}, \code{VCstem_c} and \code{VCroot_c}.}
#'   \item{Modifying \code{d} from \code{paramsTranspiration} affects simultaneously \code{VCleaf_d}, \code{VCstem_d} and \code{VCroot_d}.}
#'   \item{Modifying \code{Plant_kmax} from \code{paramsTranspiration} affects \code{VCleaf_kmax}, \code{VCstem_kmax}, \code{VCroot_kmax} and \code{belowLayers$VCroot_kmax}.}
#'   \item{Modifying \code{Al2As} from \code{paramsAnatomy} affects \code{Vsapwood} in \code{paramsWaterStorage}, \code{VCstem_kmax} and \code{VCroot_kmax} of \code{paramsTranspiration} and \code{belowLayers$VCroot_kmax}.}
#'   \item{Setting \code{paramName = "Vmax298/Jmax298"} affects both \code{Vmax298} and \code{Jmax298} from \code{paramsTranspiration}.}
#' }
#' 
#' @return If \code{summary.fun = NULL} the function returns a list whose elements are the result of calling \code{\link{spwb}}. Otherwise, the function applies \code{summary.fun} to each simulation result and returns these summaries (actually, a call to \code{\link{sapply}} is done).
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso \code{\link{spwb}}, \code{\link{summary.spwb}}
#' 
#' @examples 
#' \donttest{
#' #Load example data and species parameters
#' data(examplemeteo)
#' data(exampleforestMED)
#' data(SpParamsMED)
#' 
#' #Initialize input
#' examplesoil <- soil(defaultSoilParams(2))
#' control <- defaultControl("Granier")
#' x <- forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)
#' 
#' #Perform sensitivity analysis
#' res <- spwb_sensitivity(x, examplesoil, examplemeteo, latitude = 41, elevation = 100)
#' }
#' 
spwb_sensitivity<-function(x, soil, meteo, 
                           paramType = "above", paramName = "LAI_live", cohort = 1, 
                           p_change = c(-80,-40,-20,0,20,40,80), summary.fun = NULL, simplify=TRUE,...) {
  n = length(p_change)
  l = vector("list", n)
  names(l) = paste0(ifelse(p_change>0,"+", ""),p_change, "%")
  cat("Running spwb simulations: ")
  for(i in 1:n) {
    cat(paste0(names(l)[i]," "))
    xi = x
    xi$control$verbose= FALSE
    f = (1+(p_change[i]/100))
    .multiplyInputParam(xi, paramType, paramName, cohort - 1, f, FALSE)
    resetInputs(xi)
    l[[i]] = spwb(xi, meteo, ...)
  }
  cat("\n")
  if(!is.null(summary.fun)) {
    l = sapply(l, summary.fun, simplify=simplify)
  }
  return(l)
}