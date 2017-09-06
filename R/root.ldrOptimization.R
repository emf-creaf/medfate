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

swb.ldrOptimization<-function(x, meteo, soil, psi_crit,
                              Z95min = 301, Z95max = 4000, V1min = 0.01, V1max = 0.94, resolution = 20, 
                              heatStop = NULL, transformation = "identity", explore.out = FALSE, verbose = FALSE) {
  
  

  # define the heating period
  if(!is.null(heatStop)){
    op_days <- which(as.Date(rownames(meteo)) > heatStop)
  } else {op_days <- 1:nrow(meteo)}
  
  # Define the values of Z50 and Z90 to be explored
  trans <- function(x) do.call(transformation, list(x))
  inverse_trans <- .inverse(trans, lower = 0.01, upper = 100) # inverse of the function used for the transformation
  
  if(Z95max > soil$SoilDepth){
    if(verbose) cat("\n Z95max is larger than soil depth\n")
  }
  
  Z95_trans <- seq(trans(Z95min), trans(Z95max), length.out = resolution)
  Z95 <- unlist(sapply(Z95_trans, FUN = inverse_trans))
  
  # the case where Z95 = Z1 will create problems when using the LDR model -> remove if it exists
  Z1 <- soil$dVec[1]
  if(sum(Z95 == Z1) > 0){
    if(verbose) cat("\nThe function to derive the root proportion in each soil layer is not defined for Z95 = Z1 (depth of the first soil layer)\n",
        "This value is removed\n", 
        paste("Resolution is now\n", resolution-1))
    Z95 <- Z95[-which(Z95 == Z1)]
  }
  
  if(V1max >= 0.9497887){
    if(verbose) cat("\nThe function to derive the root proportion in each soil layer is only defined for V1 c ]0,0.949[\nV1max is set to 0.94\n")
    V1max <- 0.94  
  }
  if(V1min <= 0){
    if(verbose) cat("\nThe function to derive the root proportion in each soil layer is only defined for V1 c ]0,0.949[\nV1min is set to 0.01\n")
    V1min <- 0.01  
  }
  
  V1 <- seq(V1min,V1max,length.out = length(Z95)) # the proportion of root in the first soil layer
  # Create a matrix with V1 as rows and Z95 as column, filled with logical values to indicate the parameter combinations to explore
  mExplore <- matrix(T, nrow = length(V1), ncol = length(Z95), dimnames = list(V1 = V1, Z95 = Z95))
  # mExplore[lower.tri(mExplore, diag = T)] <- F
  
  # Calculate Z50
  Z50 <- .root.ldrZ50(V = array(V1,dim = dim(mExplore)), Z = array(Z1, dim = dim(mExplore)), Z95 = t(array(Z95, dim = dim(mExplore))))
  dimnames(Z50) <- dimnames(mExplore)
  # Prepare array for V 
  V <- array(dim = c(length(soil$dVec),length(V1), length(Z95)), 
             dimnames = list(layer = 1:length(soil$dVec), V1 = V1, Z95 = Z95))
  
  # Sum LAI of all species
  x$above$LAI_live <- sum(x$above$LAI_live) 
  x$above$LAI_expanded <- sum(x$above$LAI_expanded) 
  x$above$LAI_dead <- sum(x$above$LAI_dead)
  
  # Data outputs
  E <- PsiMin <- array(dim = c(nrow(x$above), length(V1), length(Z95)), dimnames = list(SP = x$above$SP, V1 = V1, Z95 = Z95))
  
  # rename soil to avoid conflict with swb
  # s <- soil 
  
  # Start loop
  cc <- which(mExplore == T, arr.ind = T)
  for(sp in 1:nrow(x$above)){
    SP <- x$above$SP[sp]
    cat(paste("Exploring root distribution of species", SP,"\n"))
    
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
      
      # Update the depth of the different soil layer to match Z95
      s. <- soil
      s.$SoilDepth <- Z95[j]
      dCum <- cumsum(s.$dVec)
      layersWithinZ95 <- dCum < Z95[j]
      layersWithinZ95 <- c(T,layersWithinZ95[-length(layersWithinZ95)])
      s.$dVec <- s.$dVec[layersWithinZ95] # remove the layers not included
      nl <- length(s.$dVec) #new number of layers
      s.$dVec[nl] <- s.$dVec[nl]-dCum[nl]+Z95[j] # adjust the width of the last layer
      s.$Water_FC[nl] = soil$Water_FC[nl]*(s.$dVec[nl]/soil$dVec[nl]) #Adjust volume of the last layer
      # Adjust the other soil parameters to the new number of layers
      for(k in 1:length(s.)){
        if(length(s.[[k]]) > 1) s.[[k]] <- s.[[k]][1:nl]
      }
      # print(s.)
      
      # Calculate V
      # for(l in 1:length(s.$dVec)){
      #   V[l,i,j] <- root.ldrCumul(Z50 = Z50[i,j], Z95 = Z95[j], z = sum(s.$dVec[1:l]))
      #   if(l > 1){
      #     V[l,i,j] <- V[l,i,j]-sum(V[1:(l-1),i,j])
      #   }
      # }
      # 
      # # Make the sum of the root proportions in each layers be 1 (by definition of the ldr distribution the sum can be 0.95)
      # v <- V[1:length(s.$dVec),i,j]
      # v <- v/sum(v)
      x_1sp$below$V <- root.ldrDistribution(Z50 = Z50[i,j], Z95 = Z95[j], d=s.$dVec)
      
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
                                       FUN = ma)$x)
      E[sp,i,j] <- mean(swb$PlantTranspiration[op_days])
      setTxtProgressBar(pb, row)
    }
    cat("\n")
  }
  
  # To see the variations of the root distribution in the different layers
  # V_melt <- melt(V)
  # for(l in 1:dim(V)[1]){
  #   p <- ggplot(V_melt[V_melt$layer == l,], aes(x = Z95, y = V1, fill = value, z = value)) + geom_raster()+geom_contour()
  #   print(p)
  # }
  
  
  optim <- data.frame(Z50 = NA, Z95 = NA)
  # V1 <- as.numeric(dimnames(E)$V1)
  # Z95 <- as.numeric(dimnames(E)$Z95)
  # stop()
  n_sp <- dim(E)[1]
  for (i in 1:n_sp){
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
        optim[i,] <- c(Z50 = V1[point[1]]*Z95[point[2]], Z95 = Z95[point[2]])
      }
    }
    else{
      optim[i,] <- NA
    }
  }
  if(explore.out) return(list(optim=optim, explore.out=list(E = E, PsiMin = PsiMin)))
  return(optim)
}