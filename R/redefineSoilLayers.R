redefineSoilLayers<-function(SoilParams, widths = c(300, 700, 1000, 2000)) {
  restarget = data.frame(matrix(nrow = length(widths), ncol = 6))
  names(restarget) = c("widths", "clay", "sand", "om", "bd", "rfc")
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
  }
  return(restarget)
}