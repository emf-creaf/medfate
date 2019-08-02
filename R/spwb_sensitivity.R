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
    xi[[paramType]][[paramName]][cohort] = xi[[paramType]][[paramName]][cohort]*f
    if(paramName=="LAI_live") xi[[paramType]][["LAI_expanded"]][cohort] = xi[[paramType]][["LAI_expanded"]][cohort]*f
    else if(paramName=="VCroot_kmax") {
      xi$below$VCroot_kmax[cohort, ] = xi$below$VCroot_kmax[cohort, ]*f
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