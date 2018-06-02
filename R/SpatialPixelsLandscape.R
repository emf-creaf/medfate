SpatialPixelsLandscape<-function(spxt, lct, forestlist, soillist, verbose=TRUE) {
  #check input
  if(!inherits(spxt,"SpatialPixelsTopography")) 
    stop("'spxt' has to be of class 'SpatialPixelsTopography'.")
  if(!inherits(forestlist,"list")) 
    stop("'forestlist' has to be a list of 'forest' objects.")
  if(!inherits(soillist,"list")) 
    stop("'soillist' has to be a list of 'soil' objects.")
  
  #Take grid
  grid = spxt@grid
  coords = coordinates(grid)
  elevation = spxt@data$elevation
  
  if(verbose) cat(" - Queen neighbours")
  queenNeigh = cell2nb(grid@cells.dim[1],grid@cells.dim[2], type="queen")
  class(queenNeigh)<-"list"
  
  if(verbose) cat(" - Water discharge order")
  waterOrder = order(elevation, decreasing=TRUE)
  waterQ = vector("list", length(queenNeigh))
  qfun<-function(xi, yi, zi, X, Y, Z) {
    n = length(X)
    Li = sqrt((X-xi)^2+(Y-yi)^2+(Z-zi)^2)
    dZ = zi-Z #dif. in elevation
    dZLi = dZ/Li 
    dZLi[dZ<=0] = 0 #Set to zero for neighbour cells at higher or equal elevation
    if(sum(dZLi)>0) return(dZLi/sum(dZLi))
    return(rep(0, n)) #In a flat area no discharge will be applied
  }
  for(i in 1:length(queenNeigh)) {
    wne = queenNeigh[[i]]
    waterQ[[i]] = qfun(xi = coords[i,1], yi=coords[i,2],zi = elevation[i],
                       X = coords[wne,1], Y = coords[wne,2], Z = elevation[wne])
  }  
  if(verbose) cat(" - done.\n")
  
  
  spxl = new("SpatialPixelsLandscape",
            lct = lct,
            forestlist = forestlist, 
            soillist = soillist,
            data = spxt@data,
            queenNeigh = queenNeigh,
            waterOrder = waterOrder,
            waterQ = waterQ,
            grid =grid, 
            bbox = spxt@bbox, 
            proj4string = spxt@proj4string)
  return(spxl)
}