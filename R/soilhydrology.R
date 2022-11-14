
hydrology_infiltrationAmount<-function(input, Ssoil){
  I = rep(0, length(input))
  sel = input>(0.2*Ssoil)
  I[sel] = input[sel]-((input[sel]-0.2*Ssoil)^2/(input[sel]+0.8*Ssoil))
  I[!sel] = input[!sel]
  return(I)
}

