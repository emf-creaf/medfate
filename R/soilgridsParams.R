# function to use soilgrids to estimate soil characteristics (requires GSIF)
soilgridsParams <- function(points, depths = c(300, 500, 1200), verbose = FALSE) {
  
  if(!inherits(points, "SpatialPoints")) stop("Object 'points' has to be a SpatialPoints.")
  
  coords_df = spTransform(as(points,"SpatialPoints"), sp::CRS("+proj=longlat +datum=WGS84")) # Transform to long lat
  
  npoints = length(points)
  # spatial points data frame
  # coords_df <- data.frame(lon = long, lat = lat, id = 'Site')
  # sp::coordinates(coords_df) <- ~lon+lat
  # sp::proj4string(coords_df) <- sp::CRS("+proj=longlat +datum=WGS84")
  
  # soilgrids REST API query
  query <- GSIF::REST.SoilGrids(
    attributes = c('BDTICM', 'BDRICM', 'BDRLOG', 'BLDFIE', 'CRFVOL',
                   'CLYPPT', 'SLTPPT', 'SNDPPT', 'ORCDRC')
  )
  
  if(is.null(query)) stop("Could not create GSIF query!")
  
  # data frame with retrieved values
  sg_description <- GSIF::over(query, coords_df)
  
  
  lsg = vector("list", npoints)
  names(lsg) = rownames(points@coords)
  if(verbose) pb = txtProgressBar(0, npoints, style=3)
  for(ip in 1:npoints) {
    if(verbose) setTxtProgressBar(pb, ip)
    # check if data has been retrieved
    if (sg_description[ip,1] == 'nodata') {
      warning(paste0("No data has been retrieved from SoilGrids for point ",ip,": [",paste0(round(points@coords[ip,],3), collapse=","),"]"))
    } else {
      # get default depths in mm (they come in m)
      def_depths <- abs(
        c(sg_description[ip,'depthCodesMeters.sl1'], sg_description[ip,'depthCodesMeters.sl2'],
          sg_description[ip,'depthCodesMeters.sl3'], sg_description[ip,'depthCodesMeters.sl4'],
          sg_description[ip,'depthCodesMeters.sl5'], sg_description[ip,'depthCodesMeters.sl6'],
          sg_description[ip,'depthCodesMeters.sl7']) * 1000
      )
      # get the soilGrids values for each variable
      BDTICM <- sg_description[ip,'BDTICM.BDTICM_M']*10 #From cm to mm
      BDRICM <- sg_description[ip,'BDRICM.BDRICM_M']*10 #From cm to mm
      BDRLOG <- sg_description[ip,'BDRLOG.BDRLOG_M']
      BLDFIE <- sg_description[ip, paste0('BLDFIE.M.sl', 1:7)]/1000 # from Kg/m3 to Kg/dm3
      CRFVOL <- sg_description[ip, paste0('CRFVOL.M.sl', 1:7)]
      CLYPPT <- sg_description[ip, paste0('CLYPPT.M.sl', 1:7)]
      SLTPPT <- sg_description[ip, paste0('SLTPPT.M.sl', 1:7)]
      SNDPPT <- sg_description[ip, paste0('SNDPPT.M.sl', 1:7)]
      ORCDRC <- sg_description[ip, paste0('ORCDRC.M.sl', 1:7)]/10 # from g/Kg to %
      
      # get the indexes of soilgrids layers based on user layers
      depths_indexes <- vapply(depths, function(x) {
        i = which(x - def_depths <= 0)[1]
        if(is.na(i)) i = length(def_depths)
        return(i)
      },
      numeric(1))
      
      # get the user layers number
      n_layers <- length(depths)
      
      # build the empty result object
      res <- data.frame(
        widths = c(depths[1], depths[-1] - depths[-length(depths)]),
        clay = NA,
        sand = NA,
        om = NA,
        bd= NA,
        rfc = NA
      )
      
      # iterate by desired layer
      for (layer in 1:n_layers) {
        
        # n_index (start of the soilgrids values subsetting)
        n <- ifelse(layer == 1, 1, min(depths_indexes[layer - 1] + 1, length(def_depths)))
        
        # subset the variables and def_depths values
        BLDFIE_layer <- as.numeric(BLDFIE[n:depths_indexes[layer]])
        CRFVOL_layer <- as.numeric(CRFVOL[n:depths_indexes[layer]])
        CLYPPT_layer <- as.numeric(CLYPPT[n:depths_indexes[layer]])
        SLTPPT_layer <- as.numeric(SLTPPT[n:depths_indexes[layer]])
        SNDPPT_layer <- as.numeric(SNDPPT[n:depths_indexes[layer]])
        ORCDRC_layer <- as.numeric(ORCDRC[n:depths_indexes[layer]])
        depths_layer <- as.numeric(def_depths[n:depths_indexes[layer]])
        
        # if only one value, res directly
        if (length(depths_layer) == 1) {
          res[['clay']][layer] <- CLYPPT_layer
          res[['sand']][layer] <- SNDPPT_layer
          res[['om']][layer] <- ORCDRC_layer
          res[['rfc']][layer] <- CRFVOL_layer
          res[['bd']][layer]<-BLDFIE_layer
        } else {
          # trapezoidal rule formula, vectorized
          BLDFIE_res <- (sum((depths_layer[-1] - depths_layer[-length(depths_layer)])*(BLDFIE_layer[-1] + BLDFIE_layer[-length(BLDFIE_layer)]))) / (2*(depths_layer[length(depths_layer)] - depths_layer[1]))
          CRFVOL_res <- (sum((depths_layer[-1] - depths_layer[-length(depths_layer)])*(CRFVOL_layer[-1] + CRFVOL_layer[-length(CRFVOL_layer)]))) / (2*(depths_layer[length(depths_layer)] - depths_layer[1]))
          CLYPPT_res <- (sum((depths_layer[-1] - depths_layer[-length(depths_layer)])*(CLYPPT_layer[-1] + CLYPPT_layer[-length(CLYPPT_layer)]))) / (2*(depths_layer[length(depths_layer)] - depths_layer[1]))
          SLTPPT_res <- (sum((depths_layer[-1] - depths_layer[-length(depths_layer)])*(SLTPPT_layer[-1] + SLTPPT_layer[-length(SLTPPT_layer)]))) / (2*(depths_layer[length(depths_layer)] - depths_layer[1]))
          SNDPPT_res <- (sum((depths_layer[-1] - depths_layer[-length(depths_layer)])*(SNDPPT_layer[-1] + SNDPPT_layer[-length(SNDPPT_layer)]))) / (2*(depths_layer[length(depths_layer)] - depths_layer[1]))
          ORCDRC_res <- (sum((depths_layer[-1] - depths_layer[-length(depths_layer)])*(ORCDRC_layer[-1] + ORCDRC_layer[-length(ORCDRC_layer)]))) / (2*(depths_layer[length(depths_layer)] - depths_layer[1]))
          
          # update res object
          res[['clay']][layer] <- CLYPPT_res
          res[['sand']][layer] <- SNDPPT_res
          res[['om']][layer] <- ORCDRC_res
          res[['rfc']][layer] <- CRFVOL_res
          res[['bd']][layer]<-BLDFIE_res
        }
      }
      # message("Absolute depth to bedrock : ", BDTICM, " mm\nR horizon depth (up to 2m): ", BDRICM," mm")
      lsg[[ip]] = list(
        soilparams = res,
        soilgrids_Rhorizondepth = BDRICM,
        soilgrids_Rhorizonprob = BDRLOG,
        soilgrids_absolutesoildepth = BDTICM
      )
    }
  }
  if(npoints==1) lsg = lsg[[1]]
  return(lsg)
}
