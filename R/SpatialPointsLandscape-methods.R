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
setMethod("[", signature("SpatialPointsLandscape"),definition =
            function (x, i, j, ..., drop = TRUE) 
            {
              if (!missing(j)) 
                warning("j index ignored")
              if (is.matrix(i)) 
                stop("matrix argument not supported in SpatialPointsLandscape selection")
              if (is.character(i)) 
                i <- match(i, row.names(x))
              else if (is(i, "Spatial")) 
                i = !is.na(over(x, geometry(i)))
              if (any(is.na(i))) 
                stop("NAs not permitted in row index")
              sp = as(x,"SpatialPoints")[i, , drop=drop]
              x@coords = sp@coords
              x@bbox = sp@bbox
              x@forestlist = x@forestlist[i]
              x@soillist = x@soillist[i]
              x@data = x@data[i, , ..., drop = FALSE]
              x
            }
)
setMethod("spatialForestSummary", signature("SpatialPointsLandscape"), 
          function(object, summaryFunction, ...) {
            l = object@forestlist
            if(length(l)==0) return(NULL)
            firstNoNa = which(!unlist(lapply(l,is.na)))[1]
            s = unlist(do.call(summaryFunction, args=list(object=l[[firstNoNa]],...)))
            sm = data.frame(matrix(NA, nrow=length(l), ncol=length(s)))
            colnames(sm) = names(s)
            for(i in 1:length(l)) {
              if(!is.na(l[i])) sm[i,] = unlist(do.call(summaryFunction, args=list(object=l[[i]],...)))
            }
            rownames(sm) = rownames(object@coords)
            s = sm
            return(SpatialPointsDataFrame(coords=object@coords, data = s, 
                                          proj4string=object@proj4string, 
                                          bbox = object@bbox))
          })

setMethod("spatialSoilSummary", signature("SpatialPointsLandscape"), function(object, summaryFunction, ...) {
  l = object@soillist
  if(length(l)==0) return(NULL)
  firstNoNa = which(!unlist(lapply(l,is.na)))[1]
  s = do.call(summaryFunction, args=list(object=l[[firstNoNa]],...))
  sm = data.frame(matrix(NA, nrow=length(l), ncol=length(s)))
  colnames(sm) = names(s)
  for(i in 1:length(l)) {
    if(!is.na(l[[i]])) sm[i,] = do.call(summaryFunction, args=list(object=l[[i]],...))
  }
  rownames(sm) = rownames(object@coords)
  s = sm
  return(SpatialPointsDataFrame(coords=object@coords, data = s, 
                                proj4string=object@proj4string, 
                                bbox = object@bbox))
})

.print.SpatialPointsLandscape = function(x, ..., digits = getOption("digits")) {
  cat("Object of class SpatialPointsLandscape\n")
  cat(paste("Number of points:",length(x@forestlist),"\n"))
  cc = substring(paste(as.data.frame(
    t(signif(coordinates(x), digits)))),2,999)
  df = data.frame("coordinates" = cc, x@data)
  row.names(df) = row.names(x@data)
  cat(paste("Coordinates and topography:\n"))
  print(df, ..., digits = digits)
}
setMethod("print", "SpatialPointsLandscape", function(x, ..., digits = getOption("digits")) .print.SpatialPointsLandscape(x, ..., digits))

setMethod("show", "SpatialPointsLandscape", function(object) .print.SpatialPointsLandscape(object))

.head.SpatialPointsLandscape <- function(x, n=6L, ...) {
  n <- min(n, length(x))
  ix <- sign(n)*seq(abs(n))
  x[ ix , , drop=FALSE]
}
setMethod("head", "SpatialPointsLandscape", function(x, n=6L, ...) .head.SpatialPointsLandscape(x, n, ...))

.tail.SpatialPointsLandscape <- function(x, n=6L, ...) {
  n <- min(n, length(x))
  ix <- sign(n)*rev(seq(nrow(x), by=-1L, len=abs(n)))
  x[ ix , , drop=FALSE]
}
setMethod("tail", "SpatialPointsLandscape", function(x, n=6L, ...) .tail.SpatialPointsLandscape(x, n, ...))
