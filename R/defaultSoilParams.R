defaultSoilParams<-function(n=4) {
  n = min(max(n,2),5)
  return(data.frame(
    widths = c(300,700,1000,2000,4000)[1:n],
    clay = rep(25,n),
    sand = rep(25,n),
    om = rep(NA,n),
    bd = rep(1.5,n),
    rfc = c(25,45,75,95,98)[1:n]));
}