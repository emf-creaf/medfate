spwb_sensitivity<-function(x, soil, meteo, 
                           paramType = "above", paramName = "LAI_live", cohort = NA, 
                           p_change = c(-80,-40,-20,0,20,40,80), summary.fun = NULL, simplify=TRUE,...) {
  n = length(p_change)
  l = vector("list", n)
  names(l) = paste0(ifelse(p_change>0,"+", ""),p_change, "%")
  fracRootResistance = x$control$fracRootResistance
  if(is.na(cohort)) cohort = 1:nrow(x$cohorts)
  cat("Running spwb simulations: ")
  for(i in 1:n) {
    cat(paste0(names(l)[i]," "))
    xi = x
    xi$control$verbose= FALSE
    f = (1+(p_change[i]/100))
    if(paramName=="Z50/Z95") {
      d = soil$dVec
      xi$below$Z50[cohort] = xi$below$Z50[cohort]*f
      xi$below$Z95[cohort] = xi$below$Z95[cohort]*f
      for(ci in cohort) {
        v = root_ldrDistribution(xi$below$Z50[ci], xi$below$Z95[ci], d)
        xi$below$V[ci, ] = v
        if(xi$control$transpirationMode=="Sperry") {
          xi$below$VCroot_kmax[ci,] = xi$paramsTransp$VCroot_kmax[ci]*root_xylemConductanceProportions(v, d)
          xi$below$VGrhizo_kmax[ci, ] = v*sum(xi$below$VGrhizo_kmax[ci, ])
        }
      }
    } 
    else if(paramName=="LAI_live") {
      xi$above$LAI_live[cohort] =xi$above$LAI_live[cohort]*f
      xi$above$LAI_expanded[cohort] =xi$above$LAI_expanded[cohort]*f
    } 
    else if(paramName=="Al2As") {
      xi$paramsAnatomy$Al2As[cohort] =xi$paramsAnatomy$Al2As[cohort]*f
      xi$paramsWaterStorage$Vsapwood[cohort] =xi$paramsWaterStorage$Vsapwood[cohort]/f
      xi$paramsTransp$VCstem_kmax[cohort] = xi$paramsTransp$VCstem_kmax[cohort]/f
      xi$paramsTransp$VCroot_kmax[cohort] = xi$paramsTransp$VCroot_kmax[cohort]/f
      xi$below$VCroot_kmax[cohort,] = xi$below$VCroot_kmax[cohort,]/f
      #Update plant kmax
      xi$paramsTransp$Plant_kmax[cohort] = 1/((1/xi$paramsTransp$VCleaf_kmax[cohort])+(1/xi$paramsTransp$VCstem_kmax[cohort])+(1/xi$paramsTransp$VCroot_kmax[cohort]))
    }
    else if(paramName=="WaterStorage") {
      xi$paramsWaterStorage$Vsapwood[cohort] =xi$paramsWaterStorage$Vsapwood[cohort]*f
      xi$paramsWaterStorage$Vleaf[cohort] =xi$paramsWaterStorage$Vleaf[cohort]*f
    } 
    else if(paramName=="VCroot_kmax") {
      xi$paramsTransp$VCroot_kmax[cohort] = xi$paramsTransp$VCroot_kmax[cohort]*f
      xi$below$VCroot_kmax[cohort, ] = xi$below$VCroot_kmax[cohort, ]*f
    }
    else if(paramName=="Plant_kmax") {
      xi$paramsTransp$Plant_kmax[cohort] = xi$paramsTransp$Plant_kmax[cohort]*f
      xi$below$VCroot_kmax[cohort, ] = xi$below$VCroot_kmax[cohort, ]*f
      xi$paramsTransp$VCleaf_kmax[cohort] = xi$paramsTransp$VCleaf_kmax[cohort]*f
      xi$paramsTransp$VCstem_kmax[cohort] = xi$paramsTransp$VCstem_kmax[cohort]*f
      xi$paramsTransp$VCroot_kmax[cohort] = xi$paramsTransp$VCroot_kmax[cohort]*f
    } 
    else if(paramName=="c") {
      xi$paramsTransp$VCleaf_c[cohort] = xi$paramsTransp$VCleaf_c[cohort]*f
      xi$paramsTransp$VCstem_c[cohort] = xi$paramsTransp$VCstem_c[cohort]*f
      xi$paramsTransp$VCroot_c[cohort] = xi$paramsTransp$VCroot_c[cohort]*f
    } 
    else if(paramName=="d") {
      xi$paramsTransp$VCleaf_d[cohort] = xi$paramsTransp$VCleaf_d[cohort]*f
      xi$paramsTransp$VCstem_d[cohort] = xi$paramsTransp$VCstem_d[cohort]*f
      xi$paramsTransp$VCroot_d[cohort] = xi$paramsTransp$VCroot_d[cohort]*f
    } 
    else if(paramName=="Vmax298/Jmax298") {
      xi$paramsTransp$Vmax298[cohort] = xi$paramsTransp$Vmax298[cohort]*f
      xi$paramsTransp$Jmax298[cohort] = xi$paramsTransp$Jmax298[cohort]*f
    }
    else {
      xi[[paramType]][[paramName]][cohort] = xi[[paramType]][[paramName]][cohort]*f
    }
    spwb_resetInputs(xi, soil)
    l[[i]] = spwb(xi, soil, meteo, ...)
  }
  cat("\n")
  if(!is.null(summary.fun)) {
    l = sapply(l, summary.fun, simplify=simplify)
  }
  return(l)
}