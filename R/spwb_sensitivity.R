spwb_sensitivity<-function(x, paramType = "above", paramName = "LAI_live", 
                           p_change = c(-80,-40,-20,0,20,40,80), summary.fun = NULL, simplify=TRUE,...) {
  n = length(p_change)
  l = vector("list", n)
  names(l) = paste0(ifelse(p_change>0,"+", ""),p_change, "%")
  cat("Running simulations: ")
  for(i in 1:n) {
    cat(paste0(names(l)[i]," "))
    xi = x
    xi$control$verbose= FALSE
    xi[[paramType]][[paramName]] = xi[[paramType]][[paramName]]*p_change[i]/100
    if(paramName=="LAI_live") xi[[paramType]][["LAI_expanded"]] = xi[[paramType]][["LAI_expanded"]]*p_change[i]/100
    l[[i]] = spwb(xi,...)
  }
  cat("\n")
  if(!is.null(summary.fun)) {
    l = sapply(l, summary.fun, simplify=simplify)
  }
  return(l)
}