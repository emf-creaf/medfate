# Function to calculate Z50 knowing Z95 and the cumulated root proportion V at a given depth Z
.root.ldrZ50 <- function(V,Z,Z95){
  if(sum(V >= 0.9497887)>0) {stop("The function is not defined for V >= 0.9497887")}
  if(sum(Z == Z95)>0) {stop("The function is not defined for Z = Z95")}
  a <- log(V/(1-V))/2.94
  Z50 <- (Z/Z95^a)^(1/(1-a))
  return(Z50)
}

# Convenience function for finding the inverse solution of a given function
.inverse <- function (f, lower, upper) {
  function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper, extendInt = "yes")[1]
}

swb.ldrOptimization<-function(x, soil, meteo, psi_crit,
                              RZmin = 301, RZmax = 4000, V1min = 0.01, V1max = 0.94, resolution = 20, 
                              heat_stop = 0, transformation = "identity", explore_out = FALSE, verbose = FALSE) {
  
  # define the days to keep in the analysis
  op_days <- (heat_stop+1):nrow(meteo)
  
  # Define the values of Z50 and Z90 to be explored
  trans <- function(x) do.call(transformation, list(x))
  inverse_trans <- .inverse(trans, lower = 0.01, upper = 100) # inverse of the function used for the transformation
  
  if(RZmax > soil$SoilDepth){
    if(verbose) cat("\n RZmax is larger than soil depth\n")
  }
  
  RZ_trans <- seq(trans(RZmin), trans(RZmax), length.out = resolution)
  RZ <- unlist(sapply(RZ_trans, FUN = inverse_trans))
  
  # the case where RZ = Z1 will create problems when using the LDR model -> remove if it exists
  Z1 <- soil$dVec[1]
  if(sum(RZ == Z1) > 0){
    if(verbose) cat("\nThe function to derive the root proportion in each soil layer is not defined for RZ = Z1 (depth of the first soil layer)\n",
        "This value is removed\n", 
        paste("Resolution is now\n", resolution-1))
    RZ <- RZ[-which(RZ == Z1)]
  }
  
  if(V1max >= 0.9497887){
    if(verbose) cat("\nThe function to derive the root proportion in each soil layer is only defined for V1 c ]0,0.949[\nV1max is set to 0.94\n")
    V1max <- 0.94  
  }
  if(V1min <= 0){
    if(verbose) cat("\nThe function to derive the root proportion in each soil layer is only defined for V1 c ]0,0.949[\nV1min is set to 0.01\n")
    V1min <- 0.001  
  }
  
  V1 <- seq(V1min,V1max,length.out = length(RZ)) # the proportion of root in the first soil layer
  # Create a matrix with V1 as rows and RZ as column, filled with logical values to indicate the parameter combinations to explore
  mExplore <- matrix(T, nrow = length(V1), ncol = length(RZ), dimnames = list(V1 = V1, RZ = RZ))
  # mExplore[lower.tri(mExplore, diag = T)] <- F
  
  # Calculate Z50
  Z50 <- .root.ldrZ50(V = array(V1,dim = dim(mExplore)), Z = array(Z1, dim = dim(mExplore)), Z95 = t(array(RZ, dim = dim(mExplore))))
  dimnames(Z50) <- dimnames(mExplore)
  # Prepare array for V 
  V <- array(dim = c(length(soil$dVec),length(V1), length(RZ)), 
             dimnames = list(layer = 1:length(soil$dVec), V1 = V1, RZ = RZ))
  
  # Sum LAI of all species
  x$above$LAI_live <- sum(x$above$LAI_live) 
  x$above$LAI_expanded <- sum(x$above$LAI_expanded) 
  x$above$LAI_dead <- sum(x$above$LAI_dead)
  
  # Data outputs
  E <- PsiMin <- array(dim = c(nrow(x$above), length(V1), length(RZ)), dimnames = list(SP = x$cohorts$SP, V1 = V1, RZ = RZ))
  
  # Start loop
  cc <- which(mExplore == T, arr.ind = T)
  for(sp in 1:nrow(x$above)){
    SP <- x$cohorts$SP[sp]
    cat(paste("Exploring root distribution of species", SP,"(", x$cohorts$Name[sp],"):\n"))
    
    x_1sp <- x
    x_1sp$above <- x_1sp$above[sp,]
    x_1sp$paramsBase <- x_1sp$paramsBase[sp,] 
    x_1sp$paramsTransp <- x_1sp$paramsTransp[sp,] 
    x_1sp$Transpiration <- x_1sp$Transpiration[sp] 
    x_1sp$Photosynthesis <- x_1sp$Photosynthesis[sp] 
    x_1sp$ProportionCavitated <- x_1sp$ProportionCavitated[sp] 
    x_1sp$control$verbose <- F
    
    pb <- txtProgressBar(max = nrow(cc), style = 3)
    for(row in 1:nrow(cc)){
      i <- cc[row,1]
      j <- cc[row,2]
      
      # Update the depth of the different soil layer to match RZ
      s. <- soil
      s.$SoilDepth <- RZ[j]
      dCum <- cumsum(s.$dVec)
      layersWithinRZ <- dCum < RZ[j]
      layersWithinRZ <- c(T,layersWithinRZ[-length(layersWithinRZ)])
      s.$dVec <- s.$dVec[layersWithinRZ] # remove the layers not included
      nl <- length(s.$dVec) #new number of layers
      s.$dVec[nl] <- s.$dVec[nl]-dCum[nl]+RZ[j] # adjust the width of the last layer
      s.$Water_FC[nl] = soil$Water_FC[nl]*(s.$dVec[nl]/soil$dVec[nl]) #Adjust volume of the last layer
      # Adjust the other soil parameters to the new number of layers
      for(k in 1:length(s.)){
        if(length(s.[[k]]) > 1) s.[[k]] <- s.[[k]][1:nl]
      }
      
      V[,i,j] <- x_1sp$below$V <- root.ldrDistribution(Z50 = Z50[i,j], Z95 = RZ[j], d=s.$dVec)
      
      # Run the model
      swb <- swb(x = x_1sp, meteo = meteo, soil = s.)
      
      # Outputs
      years <- substr(as.Date(rownames(meteo)), start = 1, stop = 4)
      ma <- function(x,n=10){
        f = filter(x,rep(1/n,n), method = "convolution", sides = 2)
        f = f[!is.na(f)]
        # print(sum(is.na(f)))
        f
      }
      PsiMin[sp,i,j] <- mean(aggregate(swb$PlantPsi[op_days], 
                                       by = list(years[op_days]),
                                       FUN = function(x) min(ma(x)))$x)
      E[sp,i,j] <- mean(swb$PlantTranspiration[op_days])
      setTxtProgressBar(pb, row)
    }
    cat("\n")
  }
  
  # To see the variations of the root distribution in the different layers
  # V_melt <- melt(V)
  # for(l in 1:dim(V)[1]){
  #   p <- ggplot(V_melt[V_melt$layer == l,], aes(x = RZ, y = V1, fill = value, z = value)) + geom_raster()+geom_contour()
  #   print(p)
  # }
  
  optim <- data.frame(SP = x$cohorts$SP, psi_crit = psi_crit[1:length(x$cohorts$SP)], Z50 = NA, Z95 = NA, V1 = NA)
  for (i in 1:length(x$cohorts$SP)){
    if(!is.na(psi_crit[i])){
      psimin <- PsiMin[i,,]
      e <- E[i,,]
      # emax <- max(e)
      # e[e >= emax-0.05*emax] <- emax - 0.05*emax
      # cost <- (matrix(z, ncol = 1)%*%(1-v) + matrix(300, ncol = 1, nrow = length(z))%*%v)^(3/2)
      supinf <- matrix(0, ncol = ncol(psimin), nrow = nrow(psimin))
      supinf[psimin >= psi_crit[i]] <- 1
      subselb <- rbind(supinf[-nrow(supinf),]-supinf[-1,], rep(0, ncol(supinf))) 
      subselt <- rbind(rep(0, ncol(supinf)), supinf[-1,]-supinf[-nrow(supinf),]) 
      subsell <- cbind(rep(0, nrow(supinf)), supinf[,-1]-supinf[,-ncol(supinf)])
      subselr <- cbind(supinf[,-ncol(supinf)]-supinf[,-1], rep(0, nrow(supinf)))
      sel <- matrix(F, ncol = ncol(psimin), nrow = nrow(psimin))
      sel[subselb == 1 | subselt == 1 | subsell == 1 | subselr == 1] <- T
      if(length(e[sel])==0) {
        warning(paste("Psi value", psi_crit[i],"for species",SP[i],"not reached for any combination."))
        optim[i,] <- NA
      } else {
        point <- which(sel & e == max(e[sel]), arr.ind = T)
        optim$Z50[i] <- Z50[point[1], point[2]]
        optim$V1[i] <- V1[point[1]]
        optim$Z95[i] <- RZ[point[2]]
      }
    }
    else{
      optim[i,c("Z50", "Z95")] <- NA
    }
  }
  if(explore_out) return(list(optim=optim, explore_out=list(E = E, PsiMin = PsiMin)))
  return(optim)
}

# Function for plotting the outputs of swb.ldrOptimization
# works with the libraries ggplot2, reshape and viridis
# x is the output of the function swb.ldrOptimization with explore_out = T
# .plot.ldrOptimization <- function(x, SP = 1, raster_var = "E", contour_var = "E", special_breaks_var = "Psi",
#                              legend_pos = c(1,1), xaxis_pos = "bottom", yaxis_pos = "left", special_breaks = 0, axis_trans = "identity"){
#   
#   Psi.xyz <- melt(x$explore_out$PsiMin[SP,,])
#   E.xyz <- melt(x$explore_out$E[SP,,]*365)
#   xy <- Psi.xyz[,c("V1", "RZ")]
#   
#   # Raster layer
#   if(raster_var == "Psi"){
#     leg_title <- expression(paste(Psi[min],"(MPa)"))
#     data_raster <- Psi.xyz
#   }
#   if(raster_var == "E"){
#     leg_title <- expression(paste("E (mm ", yr^{-1}, ")"))
#     data_raster <- E.xyz
#   }
#   
#   # Contour layer
#   if(contour_var == "Psi"){
#     data_contour <- Psi.xyz
#     bw1 <- 1
#     bw2 <- 0.2
#   }
#   if(contour_var == "E"){
#     data_contour <- E.xyz
#     bw1 <- 50
#     bw2 <- 10
#   }
#   
#   # Add special break
#   if(special_breaks_var == "Psi"){
#     data_special_breaks <- Psi.xyz
#   }
#   if(special_breaks_var == "E"){
#     data_special_breaks <- E.xyz
#   }
#   
#   # Optimized parameters
#   x$optim$RZ <- x$optim$Z95
#   
#   # Plot
#   p <- ggplot(xy, aes(x = RZ, y = V1))+
#     geom_raster(data = data_raster, aes(fill = value))+
#     geom_contour(data = data_contour, aes(z = value), colour = "white", binwidth = bw1, size = 1)+
#     geom_contour(data = data_contour, aes(z = value), colour = "white", binwidth = bw2, size = 0.5)+
#     geom_contour(data = data_special_breaks, aes(z = value), colour = "red", breaks = special_breaks, size = 1)+
#     geom_point(data = x$optim[SP,], aes(x = RZ, y = V1), color = "black", fill = "red", shape = 21, size = 4, inherit.aes = F)+
#     scale_fill_viridis(name = leg_title)+
#     coord_cartesian(expand = F)+
#     ylab(expression(paste(V[1])))+
#     xlab(expression(paste(RZ, "(mm)")))+
#     theme_bw()+
#     theme(legend.position = legend_pos, legend.justification = legend_pos, legend.background = element_rect(fill = rgb(1,1,1,0.7)))+
#     scale_x_continuous(position = xaxis_pos, trans = axis_trans)+
#     scale_y_continuous(position = yaxis_pos, trans = "identity")
#   
#   return(p)
# }