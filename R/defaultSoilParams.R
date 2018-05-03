defaultSoilParams<-function(n=4) {
  return(list(
    widths = c(300,700,1000,2000,4000)[1:n],
    clay = rep(25,n),
    sand = rep(25,n),
    om = rep(NA,n),
    bd = rep(1.5,n),
    macro = rep(0.1,n),
    rfc = c(20,40,60,75,90)[1:n],
    Gsoil= 0.5,
    Ksoil=0.05));
}