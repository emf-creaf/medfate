SpatialPointsLandscape<-function(spt, forestlist, soillist) {
  #check input
  if(!inherits(spt,"SpatialPointsTopography")) 
    stop("'spt' has to be of class 'SpatialPointsTopography'.")
  if(!inherits(forestlist,"list")) 
    stop("'forestlist' has to be a list of 'forest' objects.")
  if(!inherits(soillist,"list")) 
    stop("'soillist' has to be a list of 'soil' objects.")
  
  spl = new("SpatialPointsLandscape", 
            forestlist = forestlist, 
            soillist = soillist,
            data = spt@data,
            coords =spt@coords, 
            bbox = spt@bbox, 
            proj4string = spt@proj4string)
  return(spl)
}