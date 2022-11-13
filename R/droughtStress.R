.droughtstress_sim<-function(x, index = "NDD", freq = "years", bySpecies = FALSE) {
  dates = as.Date(rownames(x$Plants$PlantStress))
  ndaysTotal = length(dates)
  date.factor = cut(dates, breaks=freq)
  if(inherits(x, c("spwb","pwb"))) {
    transpMode = x$spwbInput$control$transpirationMode
  } else {
    transpMode = x$growthInput$control$transpirationMode
  }
  if("spwbInput" %in% names(x)) input = x$spwbInput
  else input = x$growthInput
  
  ndd<-function(dds) { # De Caceres et al (2015)
    return(sum(dds>0.5))
  }
  di<-function(dds) { # De Caceres et al (2015)
    return(sum(pmax(0,(dds-0.5)/0.5))/length(dds))
  }
  wsi<-function(lwp) { # Myers (1988)
    c = max(lwp, na.rm=T)
    return(abs(sum(lwp-c, na.rm=T)))
  }
  if(index=="NDD") {
    M <- apply(x$Plants$PlantStress,2,tapply, INDEX=date.factor, ndd)
  } else if(index=="DI") {
    M <- apply(x$Plants$PlantStress,2,tapply, INDEX=date.factor, di)
  } else if(index=="ADS") {
    M <- apply(x$Plants$PlantStress,2,tapply, INDEX=date.factor, function(x) {return(mean(x, na.rm=T))})
  } else if(index=="MDS") {
    M <- apply(x$Plants$PlantStress,2,tapply, INDEX=date.factor, function(x) {return(max(x, na.rm=T))})
  } else if(index=="WSI") {
    if(transpMode=="Granier") {
      M <- apply(x$Plants$PlantPsi,2,tapply, INDEX=date.factor, wsi)
    } else {
      M <- apply(x$Plants$LeafPsi,2,tapply, INDEX=date.factor, wsi)
    }
  }
  if(is.vector(M)) {
    M = t(as.matrix(M))
    rownames(M) <- levels(date.factor)
  }
  if(bySpecies) {
    cohlai = apply(x$Plants$LAI,2,max, na.rm=T)
    cohsp = as.character(input$cohorts$Name)
    lai1 = tapply(cohlai, cohsp, sum, na.rm=T)
    m1 = t(apply(sweep(M,2,cohlai,"*"),1,tapply, cohsp, sum, na.rm=T))
    if(length(lai1)==1) {
      m1 = t(m1)
    } 
    M = sweep(m1,2,lai1,"/")
    colnames(M) = names(lai1)
  }
  ncases = table(date.factor)
  M = M[ncases>0, ,drop = FALSE]
  return(M)
}

#' Drought stress indicators
#' 
#' Calculates plant drought stress indices, at different temporal scales, from simulation results.
#' 
#' @param x An object of class \code{\link{spwb}}, \code{\link{pwb}}, \code{\link{growth}} or \code{\link{fordyn}}.
#' @param index A string with the index to be calculated, either \code{"DI"}, \code{"NDD"}, \code{"ADS"}, \code{"MDS"} or \code{"WSI"} (see details).
#' @param freq Frequency of stress statistics (see \code{\link{cut.Date}}). Normally, either \code{"years"} or \code{"months"} for yearly-based or monthly-based indices.
#' @param bySpecies Allows aggregating output by species.
#' @param draw A boolean flag to indicate that a plot should be returned.
#' 
#' @details The currently available drought stress indices are:
#' \itemize{
#'   \item{\code{"ADS"}:}{ Average of daily drought stress values for the period considered.}
#'   \item{\code{"MDS"}:}{ Maximum daily drought stress during the period considered.}
#'   \item{\code{"DI"}:}{ Drought intensity, as defined in De \enc{Cáceres}{Caceres} et al. (2015).}
#'   \item{\code{"NDD"}:}{ Number of drought days, as defined in De \enc{Cáceres}{Caceres} et al. (2015).}
#'   \item{\code{"WSI"}:}{ Water stress integral, as defined in Myers (1988).}
#' }
#' 
#' @return A data frame with periods (e.g., years or months) in rows and plant cohorts (or species) in columns. 
#' Values are the calculated stress index. If \code{draw=TRUE} a ggplot is returned instead.
#' 
#' @references 
#' 
#' De \enc{Cáceres}{Caceres} M, \enc{Martínez}{Martinez}-Vilalta J, Coll L, Llorens P, Casals P, Poyatos R, Pausas JG, Brotons L. (2015) Coupling a water balance model with forest inventory data to predict drought stress: the role of forest structural changes vs. climate changes. Agricultural and Forest Meteorology 213: 77-90 (doi:10.1016/j.agrformet.2015.06.012).
#' 
#' Myers BJ (1988) Water stress integral - a link between short-term stress and long-term growth. Tree Physiology 4: 315–323 (doi: 10.1093/treephys/4.4.315)
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso \code{\link{summary.spwb}}, \code{\link{waterUseEfficiency}}
#' 
droughtStress<-function(x, index = "NDD", freq = "years", bySpecies = FALSE, draw = TRUE) {
  if(!inherits(x,c("spwb","pwb","growth", "fordyn"))) {
    stop("'x' should be of class 'spwb', 'pwb', 'growth' or 'fordyn'")
  }
  index = match.arg(index,c("NDD","DI", "ADS", "MDS","WSI"))  
  
  if(inherits(x, c("spwb","pwb", "growth"))) {
     M = .droughtstress_sim(x = x, index = index, freq = freq, bySpecies = bySpecies)
  } else {
    vec<-vector("list", length(x$GrowthResults))
    for(i in 1:length(x$GrowthResults)) {
      vec[[i]] <- .droughtstress_sim(x = x$GrowthResults[[i]], index = index, freq = freq, bySpecies = bySpecies)
    }
    M = .mergeVectorOfMatrices(vec)
  }
  if(draw) {
    if(index=="NDD") {
      g<-.multiple_dynamics(M, ylab = "Number of days with DS > 0.5")
    }
    else if(index=="DI") {
      g<-.multiple_dynamics(M, ylab = "Average drought stress over 0.5")
    }
    else if(index=="ADS") {
      g<-.multiple_dynamics(M, ylab = "Average drought stress")
    }
    else if(index=="MDS") {
      g<-.multiple_dynamics(M, ylab = "Maximum drought stress")
    }
    else if(index=="WSI") {
      g<-.multiple_dynamics(M, ylab = "Water stress integral")
    }
    return(g)
  } else {
    return(M)
  }
}
