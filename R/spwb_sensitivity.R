spwb_sensitivity<-function(x, soil, meteo, 
                           paramType = "above", paramName = "LAI_live", cohort = NA, 
                           p_change = c(-80,-40,-20,0,20,40,80), summary.fun = NULL, simplify=TRUE,...) {
  n = length(p_change)
  l = vector("list", n)
  names(l) = paste0(ifelse(p_change>0,"+", ""),p_change, "%")
  if(is.na(cohort)) cohort = 1:nrow(x$cohorts)
  cat("Running spwb simulations: ")
  for(i in 1:n) {
    cat(paste0(names(l)[i]," "))
    xi = x
    xi$control$verbose= FALSE
    f = (1+(p_change[i]/100))
    if(paramName=="LAI_live") {
      xi$above$LAI_live[cohort] =xi$above$LAI_live[cohort]*f
      xi$above$LAI_expanded[cohort] =xi$above$LAI_expanded[cohort]*f
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