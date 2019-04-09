hydrology_infiltrationAmount<-function(input, Ssoil){
  I = rep(0, length(input))
  sel = input>(0.2*Ssoil)
  I[sel] = input[sel]-((input[sel]-0.2*Ssoil)^2/(input[sel]+0.8*Ssoil))
  I[!sel] = input[!sel]
  return(I)
}
hydrology_rainInterception<-function(Rainfall, Cm, p, ER=0.05, method="Gash1995"){
  METHODS <- c("Liu2001","Gash1995")
  method <- match.arg(method, METHODS)
  if(length(ER)==1) ER =rep(ER, length(Rainfall))
  if(length(Cm)==1) Cm =rep(Cm, length(Rainfall))
  if(length(p)==1) p =rep(p, length(Rainfall))  
  
  if(method=="Gash1995") {
    PG = (-Cm/(ER*(1-p)))*log(1-ER) #Rainfall need to saturate the canopy
    PG[Cm==0 | p==1]=0 #Avoid NAs
    sel = Rainfall > PG #Days where the canopy becomes saturated
    I = rep(NA,length(Rainfall))
    I[sel] = (1-p[sel])*PG[sel] + (1-p[sel])*ER[sel]*(Rainfall[sel]-PG[sel]) 
    I[!sel] = (1-p[!sel])*Rainfall[!sel]  
  } else if(method=="Liu2001") {
    I = Cm*(1-exp(-1*(Rainfall)*((1-p)/Cm)))*(1-(ER/(1-p)))+(ER*Rainfall)
  }
  return(I)
}
