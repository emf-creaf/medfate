# Function to calculate Z50 knowing Z95 and the cumulated root proportion V at a given depth Z
.root_ldrZ50 <- function(V,Z,Z95){
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

spwb_ldrExploration<-function(x, soil, meteo, cohorts = NULL, 
                              RZmin = 301, RZmax = 4000, V1min = 0.01, V1max = 0.94, resolution = 10, 
                              heat_stop = 0, transformation = "identity", verbose = FALSE, 
                              ...) {
  # define the days to keep in the analysis
  op_days <- (heat_stop+1):nrow(meteo)
  
  # Define the values of Z50 and Z90 to be explored
  trans <- function(x) do.call(transformation, list(x))
  inverse_trans <- .inverse(trans, lower = 0.01, upper = 100) # inverse of the function used for the transformation
  
  if(RZmax > soil$SoilDepth){
    if(verbose) cat("\n RZmax is larger than soil depth\n")
  }
  
  RZ_trans <- seq(trans(RZmin), trans(RZmax), length.out = resolution)
  RZ <- as.numeric(unlist(sapply(RZ_trans, FUN = inverse_trans)))
  
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
  Z50 <- .root_ldrZ50(V = array(V1,dim = dim(mExplore)), Z = array(Z1, dim = dim(mExplore)), Z95 = t(array(RZ, dim = dim(mExplore))))
  dimnames(Z50) <- dimnames(mExplore)
  # Prepare array for V 
  V <- array(dim = c(length(soil$dVec),length(V1), length(RZ)), 
             dimnames = list(layer = 1:length(soil$dVec), V1 = V1, RZ = RZ))
  
  # Sum LAI of all species
  x$above$LAI_live <- sum(x$above$LAI_live) 
  x$above$LAI_expanded <- sum(x$above$LAI_expanded) 
  x$above$LAI_dead <- sum(x$above$LAI_dead)
  
  # Data outputs
  if(is.null(cohorts)) cohorts = row.names(x$cohorts)
  An<-E <- PsiMin <- array(dim = c(length(cohorts), length(V1), length(RZ)), dimnames = list(cohort = cohorts, V1 = V1, RZ = RZ))
  
  # Start loop
  cc <- which(mExplore == T, arr.ind = T)
  
  # Reset input
  spwb_resetInputs(x, soil)
  
  for(ci in 1:length(cohorts)){
    coh = cohorts[ci]
    sp = which(row.names(x$cohorts)==coh)
    
    cat(paste("Exploring root distribution of cohort", coh,"(", x$cohorts$Name[sp],"):\n"))
    
    x_1sp <- x
    x_1sp$cohorts <- x$cohorts[sp,,drop = FALSE]
    x_1sp$above <- x$above[sp,,drop = FALSE]
    x_1sp$below <- x$below
    x_1sp$below$V <- x$below$V[sp,,drop = FALSE] 
    x_1sp$paramsBase <- x$paramsBase[sp,,drop = FALSE] 
    x_1sp$paramsTransp <- x$paramsTransp[sp,,drop = FALSE] 
    x_1sp$Transpiration <- x$Transpiration[sp,drop = FALSE] 
    x_1sp$Photosynthesis <- x$Photosynthesis[sp,drop = FALSE] 
    if(x_1sp$control$transpirationMode=="Granier") {
      x_1sp$PLC <- x$PLC[sp,drop = FALSE] 
    } else {
      x_1sp$below$VGrhizo_kmax <- x$below$V[sp,,drop = FALSE] 
      x_1sp$below$VCroot_kmax <- x$below$V[sp,,drop = FALSE] 
      x_1sp$paramsAnatomy <- x$paramsAnatomy[sp,,drop = FALSE] 
      x_1sp$paramsWaterStorage <- x$paramsWaterStorage[sp,,drop = FALSE] 
      x_1sp$PLCstem <- x$PLCstem[sp,drop = FALSE] 
      x_1sp$Einst <- x$Einst[sp,drop = FALSE] 
      x_1sp$psiRhizo <- x$psiRhizo[sp,,drop = FALSE] 
      x_1sp$psiRootCrown <- x$psiRootCrown[sp,drop = FALSE] 
      x_1sp$psiSympStem <- x$psiSympStem[sp,drop = FALSE] 
      x_1sp$psiStem1 <- x$psiStem1[sp,drop = FALSE] 
      x_1sp$psiStem2 <- x$psiStem2[sp,drop = FALSE] 
      x_1sp$psiSympLeaf <- x$psiSympLeaf[sp,drop = FALSE] 
      x_1sp$psiLeaf <- x$psiLeaf[sp,drop = FALSE] 
    }
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
      # s.$Water_FC[nl] = soil$Water_FC[nl]*(s.$dVec[nl]/soil$dVec[nl]) #Adjust volume of the last layer
      # Adjust the other soil parameters to the new number of layers
      s.[["sand"]] <- s.[["sand"]][1:nl]
      s.[["clay"]] <- s.[["clay"]][1:nl]
      s.[["om"]] <- s.[["om"]][1:nl]
      s.[["rfc"]] <- s.[["rfc"]][1:nl]
      s.[["macro"]] <- s.[["macro"]][1:nl]
      s.[["W"]] <- s.[["W"]][1:nl]
      s.[["Temp"]] <- s.[["Temp"]][1:nl]
      s.[["VG_alpha"]] <- s.[["VG_alpha"]][1:nl]
      s.[["VG_theta_res"]] <- s.[["VG_theta_res"]][1:nl]
      s.[["VG_theta_sat"]] <- s.[["VG_theta_sat"]][1:nl]
      s.[["usda_Type"]] <- s.[["usda_Type"]][1:nl]
      
      V[,i,j] <- 0
      x_1sp$below$V = x$below$V[sp,1:nl,drop = FALSE]
      x_1sp$below$V[1,] <- root_ldrDistribution(Z50 = Z50[i,j], Z95 = RZ[j], d=s.$dVec)
      V[1:length(x_1sp$below$V),i,j] <- x_1sp$below$V
      if(x_1sp$control$transpirationMode=="Sperry"){
        x_1sp$below$VCroot_kmax = x$below$VCroot_kmax[sp,1:nl,drop = FALSE]
        x_1sp$below$VGrhizo_kmax = x$below$VGrhizo_kmax[sp,1:nl,drop = FALSE]
        x_1sp$below$VCroot_kmax[1,] = x_1sp$paramsTransp$VCroot_kmax*root_xylemConductanceProportions(x_1sp$below$V, s.$dVec)
        x_1sp$below$VGrhizo_kmax[1,] = x_1sp$below$V*sum(x$below$VGrhizo_kmax[sp,])
      }
      
      s_res <- spwb(x = x_1sp, meteo = meteo, soil = s., ...)
      
      # Outputs
      years <- substr(as.Date(rownames(meteo)), start = 1, stop = 4)
      ma <- function(x,n=10){
        f = filter(x,rep(1/n,n), method = "convolution", sides = 2)
        f = f[!is.na(f)]
        # print(sum(is.na(f)))
        f
      }
      if(x_1sp$control$transpirationMode=="Granier") {
        PsiMin[ci,i,j] <- mean(aggregate(s_res$PlantPsi[op_days], 
                                         by = list(years[op_days]),
                                         FUN = function(x) min(ma(x)))$x)
      } else {
        PsiMin[ci,i,j] <- mean(aggregate(s_res$StemPsi[op_days], 
                                         by = list(years[op_days]),
                                         FUN = function(x) min(ma(x)))$x)
      }
      # if(verbose) print(s_res$spwbInput)
      E[ci,i,j] <- mean(s_res$PlantTranspiration[op_days], na.rm=T)
      An[ci,i,j] <- mean(s_res$PlantPhotosynthesis[op_days], na.rm=T)
      setTxtProgressBar(pb, row)
    }
    cat("\n")
  }
  res <-list(cohorts = cohorts, RZ = RZ, V1 = V1, Z50 = Z50, E = E, An = An, PsiMin = PsiMin)
  class(res)<-list("spwb_ldrExploration","list")
  return(res)
}

spwb_ldrOptimization<-function(y, psi_crit, opt_mode = 1) {
  

  E = y$E
  An = y$An
  PsiMin = y$PsiMin
  V1 = y$V1
  RZ = y$RZ
  Z50 = y$Z50
  cohorts = y$cohorts

  if(length(psi_crit)!= length(cohorts)) stop("The length of 'psi_crit' must be equal to the number of cohorts in 'y'.")
  
  optim <- data.frame(psi_crit = psi_crit, Z50 = NA, Z95 = NA, V1 = NA)
  row.names(optim) = cohorts
  for (i in 1:length(cohorts)){
    psimin <- PsiMin[i,,]
    e <- E[i,,]
    an <- An[i,,]
    if(opt_mode==1) {
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
        warning(paste("Psi value", psi_crit[i],"for cohort ",row.names(cohorts)[i],"not reached for any combination."))
        optim[i,] <- NA
      } else {
        point <- which(sel & e == max(e[sel]), arr.ind = T)
        optim$Z50[i] <- Z50[point[1], point[2]]
        optim$V1[i] <- V1[point[1]]
        optim$Z95[i] <- RZ[point[2]]
      }
    } 
    else if(opt_mode==2) {
      selPsi = (psimin > psi_crit[i]) # Select combinations with less stress than psi_crit
      if(sum(selPsi)==0) selPsi = (psimin == max(psimin)) # If none, select combination of minimum stress
      maxE = max(e[selPsi]) # Find maximum transpiration (among combinations selected)
      sel2 = selPsi & (e == maxE) # Select combination with maximum transpiration
      point = which(sel2, arr.ind = TRUE)
      optim$Z50[i] <- Z50[point[1], point[2]]
      optim$V1[i] <- V1[point[1]]
      optim$Z95[i] <- RZ[point[2]]
    } 
    else if(opt_mode==3) {
      selPsi = (psimin > psi_crit[i]) # Select combinations with less stress than psi_crit
      if(sum(selPsi)==0) selPsi = (psimin == max(psimin)) # If none, select combination of minimum stress
      maxAn = max(an[selPsi]) # Find maximum transpiration (among combinations selected)
      sel2 = selPsi & (an == maxAn) # Select combinations with maximum photosynthesis
      point = which(sel2, arr.ind = TRUE)
      optim$Z50[i] <- Z50[point[1], point[2]]
      optim$V1[i] <- V1[point[1]]
      optim$Z95[i] <- RZ[point[2]]
    }
    else if(opt_mode==4) {
      selPsi = (psimin > psi_crit[i]) # Select combinations with less stress than psi_crit
      if(sum(selPsi)==0) selPsi = (psimin == max(psimin)) # If none, select combination of minimum stress
      maxE = max(e[selPsi]) # Find maximum transpiration (among combinations selected)
      sel2 = selPsi & (e >= maxE*0.95) # Select combinations with > 95% of maximum transpiration
      points = as.data.frame(which(sel2, arr.ind = TRUE))
      minZ = min(points$RZ, na.rm=T) # Minimum rooting depth
      maxV1 = max(points$V1[points$RZ==minZ], na.rm=T) # Maximum V1
      point = c(maxV1, minZ)
      optim$Z50[i] <- Z50[point[1], point[2]]
      optim$V1[i] <- V1[point[1]]
      optim$Z95[i] <- RZ[point[2]]
    } 
    else if(opt_mode==5) {
      selPsi = (psimin > psi_crit[i]) # Select combinations with less stress than psi_crit
      if(sum(selPsi)==0) selPsi = (psimin == max(psimin)) # If none, select combination of minimum stress
      maxAn = max(an[selPsi]) # Find maximum transpiration (among combinations selected)
      sel2 = selPsi & (an >= maxAn*0.95) # Select combinations with > 95% of maximum photosynthesis
      points = as.data.frame(which(sel2, arr.ind = TRUE))
      minZ = min(points$RZ, na.rm=T) # Minimum rooting depth
      maxV1 = max(points$V1[points$RZ==minZ], na.rm=T) # Maximum V1
      point = c(maxV1, minZ)
      optim$Z50[i] <- Z50[point[1], point[2]]
      optim$V1[i] <- V1[point[1]]
      optim$Z95[i] <- RZ[point[2]]
    } 
  }
  return(optim)
}

# Function for plotting the outputs of spwb_ldrOptimization
# works with the libraries ggplot2, reshape and viridis
# x is the output of the function spwb_ldrOptimization with explore_out = T
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