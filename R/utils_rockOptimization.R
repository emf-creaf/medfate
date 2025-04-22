#' Optimization of rock fragment content
#'
#' Function \code{utils_rockOptimization} finds optimum rock fragment content in the soil
#' corresponding to given vegetation, weather and target percent loss
#' of conductance (PLC), following a modification of the method proposed by Druel et al. (2023).
#'
#' @param x An object of class \code{\link[medfate]{forest}}.
#' @param SpParams A data frame with species parameters (see \code{\link[medfate]{SpParamsDefinition}} and \code{\link[medfate]{SpParamsMED}}).
#' @param soil An object of class \code{\link{data.frame}} or \code{\link[medfate]{soil}}, containing soil parameters per soil layer.
#' @param control A list with default control parameters (see \code{\link[medfate]{defaultControl}}).
#' @param meteo A data frame with daily meteorological data series (see \code{\link[medfate]{spwb}}).
#' @param PLCquantile Maximum PLC quantile to be calculated across years.
#' @param qPLC_target Target PLC to be achieved (by default 12%).
#' @param qPLC_tol Tolerance of PLC difference to target accepted when finding solution.
#' @param sew_min Minimum soil extractable water (mm) for rock exploration.
#' @param max_rocks Maximum content in coarse fragments allowed for any soil layer.
#' @param verbose A logical value. Print the internal messages of the function?
#' @param ... Additional parameters to function \code{\link[medfate]{spwb}}.
#'
#' @details
#' The function performs a model inversion based on an ecohydrological assumption,
#' consisting in that forest leaf area index is in equilibrium with a low embolism
#' rate under normal conditions. This is translated in that the (by default 90%) interannual quantile of
#' the maximum annual percent loss of conductance (PLC), averaged over plant cohorts,
#' should be close to a target PLC value (by default 12%).
#' 
#' The algorithm first determines the PLC corresponding to the minimum and maximum soil extractable water (SEW). 
#' The minimum SEW (SEW_min) is an input parameter, whereas the maximum SEW (SEW_max) corresponds to no rock fragments in the soil. 
#' 
#' Then three situations are distinguished:
#' \enumerate{
#'   \item{If \code{PLC(SEW_min) < qPLC_target}  and \code{PLC(SEW_max) < qPLC_target}, the function will use \code{\link{uniroot}} to find the root
#' of the function \code{f(x) = PLC(x) - qPLC_target}, where \code{x} is SEW, which corresponds to a factor that multiplies the original rock fragment content.}
#'   \item{If both \code{PLC(SEW_min) < qPLC_target} and \code{PLC(SEW_max) < qPLC_target}, the function cannot find an optimum, because PLC is always too low, and will return the original rock fragment content}
#'   \item{Analogously, if both \code{PLC(SEW_min) > qPLC_target} and \code{PLC(SEW_max) > qPLC_target}, the function cannot find an optimum, because PLC is always too large, and will return the original rock fragment content}
#' } 
#'
#' @return
#' Function \code{utils_rockOptimization} returns a list containing:
#'  \itemize{
#'    \item{\code{RFC}: A vector with the estimated rock fragment content for each soil layer.}
#'    \item{\code{SEW}: Soil extractable water (mm).}
#'    \item{\code{runs}: Number of simulations performed.}
#'    \item{\code{message}: Text message indicating whether optimization could be done (OK) or not.}
#'  }
#'
#'
#' @references
#' Druel, A., Martins, N., Cochard, H., De Caceres, M., Delzon, S., Mencuccini, M., Torres-Ruiz, J., and Ruffault, J.: European forest vulnerability to hydraulic failure: an ecohydrological approach, EGU General Assembly 2023, Vienna, Austria, 24–28 Apr 2023, EGU23-17068, https://doi.org/10.5194/egusphere-egu23-17068, 2023.
#'
#' @author
#' \enc{Arsène}{Arsene} Druel, URFM-INRAE
#'
#' Nicolas Martin-StPaul, URFM-INRAE
#'
#' Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#'
#' @seealso \code{\link{spwb}}, \code{\link{soil}}, \code{\link{utils_ldrOptimization}}
#'
#' @examples
#' \donttest{
#' #Load example daily meteorological data
#' data(examplemeteo)
#'
#' #Load example plot plant data
#' data(exampleforest)
#'
#' #Default species parameterization
#' data(SpParamsMED)
#'
#' #Initialize soil with two layers
#' examplesoil <- defaultSoilParams(4)
#'
#' #Rock fragment content optimization (Granier)
#' utils_rockOptimization(exampleforest, soil = examplesoil,
#'                        SpParams = SpParamsMED, meteo = examplemeteo,
#'                        control = defaultControl("Granier"),
#'                        elevation = 100, latitude = 41.82592)
#' }
#' @export
utils_rockOptimization<- function(x, soil, SpParams, control, meteo,
                                 PLCquantile = 0.9, qPLC_target = 12, qPLC_tol = 0.5,
                                 sew_min = 30, max_rocks = 99, verbose = FALSE, ...){
  
  control$verbose = FALSE
  control$leafCavitationRecovery = "annual"
  control$stemCavitationRecovery = "annual"
  if(!inherits(soil, "soil")) {
    soil <- medfate::soil(soil)
  }
  nlayers <- nrow(soil)
  
  LAI_max_coh <- medfate::plant_LAI(x, SpParams)
  LAI_max <- sum(LAI_max_coh, na.rm = TRUE)
  
  sew_ori <- sum(medfate::soil_waterExtractable(soil, model = control$soilFunctions))
  
  soil_max <- soil
  soil_max$rfc <- rep(0, nlayers)
  sew_max <- sum(medfate::soil_waterExtractable(soil_max, model = control$soilFunctions))
  if(verbose) cat(paste0("SEW original = ", round(sew_ori), " minimum = ", round(sew_min), " maximum = ", round(sew_max),"\n"))

  # Function to be optimized for factor corresponding to sew_target
  f_sew_diff <- function(factor, sew_target) {
    soil_tmp <- soil
    soil_tmp$rfc <- pmax(pmin(soil_tmp$rfc*factor,max_rocks),0)
    sew <- sum(medfate::soil_waterExtractable(soil_tmp, model = control$soilFunctions))
    return(sew_target - sew)
  }
  
  f_PLC_diff<-function(sew_target, ...) {
    soil_new <- soil
    r <- uniroot(f_sew_diff, c(0,10), sew_target)
    soil_new$rfc <- pmax(pmin(soil_new$rfc*r$root,max_rocks),0)
    input_new <- medfate::spwbInput(x, soil = soil_new, SpParams = SpParams, control = control)
    
    # Launch simulation
    S_new <- spwb(x = input_new, meteo = meteo, ...)
    # 90% quantile by species of annual maximum PLC
    PLC_new <- 100*apply(summary(S_new, output="StemPLC", FUN = max),2,quantile, prob = PLCquantile)
    if(length(PLC_new)>0) {
      PLC_av_new <- sum(PLC_new*LAI_max_coh, na.rm = TRUE)/LAI_max
    } else {
      PLC_av_new <- 0
    }
    return(PLC_av_new - qPLC_target)
  }
  
  # Evaluate function at extremes
  plc_sew_min <- f_PLC_diff(sew_min, ...)
  plc_sew_max <- f_PLC_diff(sew_max, ...)
  if(verbose) cat(paste0("PLC(SEW minimum) = ", round(plc_sew_min + qPLC_target, 2), " PLC(SEW maximum) = ", round(plc_sew_max + qPLC_target),"\n"))
  
  runs <- 2
  
  # Normal situation (negative and positive extremes), find root
  message = "OK"
  if((plc_sew_max < 0.0) && (plc_sew_min > 0.0)) {
    a <-uniroot(f_PLC_diff, c(sew_min, sew_max), tol = qPLC_tol, ...)
    sew_target <- a$root
    runs <- runs + a$iter
  } else if(plc_sew_max==0.0) { #Unlikely but possible
    sew_target <- sew_max
  } else if(plc_sew_min==0.0) { #Unlikely but possible
    sew_target <- sew_min
  } else if((plc_sew_max > 0.0) && (plc_sew_min > 0.0)) {
    message = paste0("PLC larger than target for the whole range of SEW. Returning original SEW.")
    sew_target <- sew_ori
  } else if((plc_sew_max < 0.0) && (plc_sew_min < 0.0)) {
    message = paste0("PLC lower than target for the whole range of SEW. Returning original SEW.")
    sew_target <- sew_ori
  }
  
  f_target <- uniroot(f_sew_diff, c(0,10), sew_target)
  rfc_target <- pmax(pmin(soil$rfc*f_target$root,max_rocks),0)
  res <- list(RFC = rfc_target, 
              SEW = sew_target, 
              runs = runs,
              message = message)
  return(res)
}
