spwb_ldrCalibration <- function(x, soil, meteo, calibVar, obs,
                               RZmin = 301, RZmax = 4000, V1min = 0.01,
                               V1max = 0.94, resolution = 20, heat_stop = 0,
                               transformation = "identity", verbose = FALSE) {
  
  # match calibVar argument
  calibVar <- match.arg(calibVar, c('SWC', 'Eplanttot', 'Cohorts'))
  
  ## MAE calculator
  .MAE <- function(measured, predicted) {
    # calculate MAE
    res <- mean(abs(measured - predicted), na.rm = TRUE)
    # return it
    return(res)
  }
  
  # define the days to keep the analysis
  op_days <- (heat_stop + 1):nrow(meteo)
  
  # define Z50 and Z95 values to be explored
  trans <- function(x) do.call(transformation, list(x))
  ## inverse of the function used for the transformation
  inverse_trans <- .inverse(trans, lower = 0.01, upper = 100)
  ## Inform if the RZmax is deeper than the SoilDepth
  if (RZmax > soil$SoilDepth) {
    if (verbose) {
      cat("\n RZmax is larger than soil depth\n")
    }
  }
  ## RZ vals
  RZ_trans <- seq(trans(RZmin), trans(RZmax), length.out = resolution)
  RZ <- unlist(sapply(RZ_trans, FUN = inverse_trans))
  ## the case where RZ = Z1 will create problems when using
  ## the LDR model -> remove if it exists
  Z1 <- soil$dVec[1]
  if (sum(RZ == Z1) > 0) {
    if (verbose) {
      cat("\nThe function to derive the root proportion in each soil layer is ",
          "not defined for RZ = Z1 (depth of the first soil layer)\n",
          "This value is removed\n", 
          paste("Resolution is now\n", resolution - 1))
    }
    RZ <- RZ[-which(RZ == Z1)]
  }
  ## check V1min and V1max values
  if (V1max >= 0.9497887) {
    if (verbose) {
      cat("\nThe function to derive the root proportion in each soil layer is ",
          "only defined for V1 c ]0,0.949[\nV1max is set to 0.94\n")
    }
    V1max <- 0.94  
  }
  if (V1min <= 0) {
    if (verbose) {
      cat("\nThe function to derive the root proportion in each soil layer is ",
          "only defined for V1 c ]0,0.949[\nV1min is set to 0.01\n")
    }
    V1min <- 0.001  
  }
  ## V1 vals (proportion of root in the first soil layer)
  V1 <- seq(V1min,V1max,length.out = length(RZ))
  ## We create a matrix indicating the combinations to explore, and also to help
  ## calculate Z50 vals
  mExplore <- matrix(TRUE, nrow = length(V1), ncol = length(RZ),
                     dimnames = list(V1 = V1, RZ = RZ))
  ## Get the Z50 vals
  Z50 <- .root_ldrZ50(
    V = array(V1, dim = dim(mExplore)),
    Z = array(Z1, dim = dim(mExplore)),
    Z95 = t(array(RZ, dim = dim(mExplore)))
  )
  dimnames(Z50) <- dimnames(mExplore)
  
  # Prepare V array
  V <- array(dim = c(length(soil$dVec),length(V1), length(RZ)), 
             dimnames = list(layer = 1:length(soil$dVec), V1 = V1, RZ = RZ))
  
  # store the original input to be able to access the cohorts LAI in case it
  # is needed (calibVar = 'Cohorts')
  old_x <- x
  # Sum LAI for all the species
  x[['above']][['LAI_live']] <- sum(x[['above']][['LAI_live']])
  x[['above']][['LAI_expanded']] <- sum(x[['above']][['LAI_expanded']])
  x[['above']][['LAI_dead']] <- sum(x[['above']][['LAI_dead']])
  
  # Outputs
  res_by_sp <- data.frame(
    SP = vector(),
    Z95 = vector(),
    Z50 = vector(),
    V1 = vector(),
    MAE = vector()
  )
  
  # Let's go to the looping thing
  indexes_to_explore <- which(mExplore == TRUE, arr.ind = TRUE)
  
  # Main loop by cohorts
  for (sp in 1:nrow(x$above)) {
    
    # check if shrub (not doing anything at the moment with shrubs)
    sp_code <- rownames(x[['cohorts']])[sp]
    if (strtrim(sp_code, 1) == 'S') {
      res_by_sp[sp, ] <- c(
        x[['cohorts']][['SP']][sp],
        NA, NA, NA, NA
      )
      
      next()
    }
    
    # get the sp
    SP <- x[['cohorts']][['SP']][sp]
    cat(paste("Exploring root distribution of cohort",
              row.names(x[['cohorts']])[sp],
              "(", x[['cohorts']][['Name']][sp], "):\n"))
    
    # adapting the input object
    x_1sp <- x
    x_1sp$cohorts <- x_1sp[['cohorts']][sp,]
    x_1sp$above <- x_1sp[['above']][sp,]
    x_1sp$paramsInterception <- x_1sp[['paramsInterception']][sp,] 
    x_1sp$paramsTranspiration <- x_1sp[['paramsTranspiration']][sp,] 
    x_1sp$control$verbose <- FALSE
    
    # temp output
    if (calibVar == 'SWC') {
      MAE_res <- array(
        dim = c(1, length(V1)),
        dimnames = list(MAE = 'MAE', V1 = V1)
      )
      
      Z50_vals <- array(
        dim = c(1, length(V1)),
        dimnames = list(Z50 = 'Z50', V1 = V1)
      )
      
    } else {
      MAE_res <- array(
        dim = c(1, length(V1), length(RZ)),
        dimnames = list(MAE = 'MAE', V1 = V1, RZ = RZ)
      )
    }
    
    # Inner loop to iterate by Z and V combinations using indexes_to_explore
    for (row in 1:nrow(indexes_to_explore)) {
      
      i <- indexes_to_explore[row, 1]
      j <- indexes_to_explore[row, 2]
      
      ### Id calib and valid are done by SWC
      if (calibVar == 'SWC') {
        
        # calculate Z50 assuming Z95 as soil depth and Z as depth of first layer
        Z50_val <- .root_ldrZ50(
          V = V1[i],
          Z = soil$dVec[1],
          Z95 = soil$SoilDepth
        )
        # Calculate the V and modify x accordingly
        x_1sp[['below']][['V']] <- root_ldrDistribution(
          Z50 = Z50_val,
          Z95 = soil$SoilDepth,
          d = soil$dVec
        )
        
        # Run the model
        res_model <- spwb(x = x_1sp, meteo = meteo, soil = soil)
        
        # calculate the MAE
        predicted <- res_model[['SoilWaterBalance']][['W.1']]*soil[["Theta_FC"]][[1]]
        
        # Store the MAE value
        MAE_res[,i] <- .MAE(obs, predicted)
        Z50_vals[,i] <- Z50_val
        
      } else {
        
        ### Id calib and valid are done by E
        # Here we need to update the depth of the different soil layers to match
        # RZ value
        s. <- soil
        s.$SoilDepth <- RZ[j]
        dCum <- cumsum(s.$dVec)
        layersWithinRZ <- dCum < RZ[j]
        layersWithinRZ <- c(TRUE, layersWithinRZ[-length(layersWithinRZ)])
        # maintain only the layers included
        s.$dVec <- s.$dVec[layersWithinRZ]
        nl <- length(s.$dVec)
        # modify the width of the last layer
        s.$dVec[nl] <- s.$dVec[nl] - dCum[nl] + RZ[j]
        # adjust the water content of the last layer
        s.$Water_FC[nl] <- soil$Water_FC[nl]*(s.$dVec[nl]/soil$dVec[nl])
        # adjust the other soil parameters to the new number of layers
        for (k in 1:length(s.)) {
          if (length(s.[[k]]) > 1) {
            s.[[k]] <- s.[[k]][1:nl]
          }
        }
        
        # Calculate the V and modify x accordingly
        x_1sp[['below']][['V']] <- root_ldrDistribution(
          Z50 = Z50[i,j],
          Z95 = RZ[j],
          d = s.$dVec
        )
        
        # Run the model
        res_model <- spwb(x = x_1sp, meteo = meteo, soil = s.)
        # build a data frame with transpiration values to compare to obs
        E_df <- data.frame(Eplanttot = res_model[["DailyBalance"]][["Eplanttot"]])
        E_df <- cbind(E_df, res_model[["PlantTranspiration"]])
        
        if (calibVar == 'Eplanttot') {
          # calculate the MAE
          real <- obs[['Eplanttot']]
          predicted <- E_df[['Eplanttot']]
          
          # Store the MAE value
          MAE_res[,i,j] <- .MAE(real, predicted)
        } else {
          # calculate the MAE
          real <- obs[[1 + sp]]
          # scaling the LAI to be able to compare (here the cohort analyzed is always
          # the second column of E_df, as we create an input with only one cohort)
          predicted <- E_df[[2]] * old_x[['above']][['LAI_live']][sp] / x[['above']][['LAI_live']][sp]
          
          # Store the MAE value
          MAE_res[,i,j] <- .MAE(real, predicted)
        }
      }
    }
    
    # build the res data frame
    if (calibVar == 'SWC') {
      
      mae_min_index <- which(MAE_res == min(MAE_res, na.rm = TRUE), arr.ind = TRUE)
      
      res_by_sp[sp, 'SP'] <- SP
      res_by_sp[sp, 'MAE'] <- MAE_res[mae_min_index[1], mae_min_index[2]]
      res_by_sp[sp, 'Z95'] <- soil$SoilDepth
      res_by_sp[sp, 'Z50'] <- Z50_vals[mae_min_index[1], mae_min_index[2]]
      res_by_sp[sp, 'V1'] <- V1[mae_min_index[2]]
      
    } else {
      
      mae_min_index <- which(MAE_res == min(MAE_res, na.rm = TRUE), arr.ind = TRUE)
      
      res_by_sp[sp, 'SP'] <- SP
      res_by_sp[sp, 'MAE'] <- MAE_res[mae_min_index[1], mae_min_index[2], mae_min_index[3]]
      res_by_sp[sp, 'Z95'] <- RZ[mae_min_index[3]]
      res_by_sp[sp, 'Z50'] <- Z50[mae_min_index[2], mae_min_index[3]]
      res_by_sp[sp, 'V1'] <- V1[mae_min_index[2]]
    }
  }
  
  return(res_by_sp)
}