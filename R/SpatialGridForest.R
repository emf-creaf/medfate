SpatialGridForest<-function(landTopo, SFItreeData, SFIshrubData, 
                            SpatialPointsIDs, SpParams, SoilParamData = NULL, 
                            lctInput , forestLCTs, 
                            shrublandLCTs= numeric(0), grasslandLCTs= numeric(0), agricultureLCTs= numeric(0), rockLCTs= numeric(0), staticLCTs= numeric(0), 
                            SFIcodes=NULL, FMcodes = NULL, control = defaultControl()){

  #Check SpatialPointsIDs
  if(!inherits(SpatialPointsIDs, "SpatialPoints")) stop("'SpatialPointsIDs' has to be a 'SpatialPoints' object")
  xyplots = SpatialPointsIDs@coords
  IDs = rownames(xyplots)
  if(is.null(IDs)) stop("SpatialPoints must have forest plot IDs in coordinate row names.")
  
  
  #Check landTopo
  if(!(inherits(landTopo, "SpatialGridTopography"))) stop("'landTopo' has to be a 'SpatialGridTopography' object")

  grid = landTopo@grid
  elevation = landTopo@elevation
  slope = landTopo@slope
  aspect = landTopo@aspect
  
  #Cell size
  patchsize = prod(grid@cellsize)  
  coords = coordinates(grid)
  nCells = nrow(coords)
  sg = SpatialGrid(grid)
  
  if(control$verbose) {
    cat(paste("Grid cells: ", nCells,", cell size: ", patchsize, " m2,  area: ", 
              areaSpatialGrid(sg)/10000," ha\n", sep=""))
  }
  
  #Check SoilParamData
  if(!is.null(SoilParamData)) if(nCells!=nrow(SoilParamData)) stop("The number of grid cells has to be the same as the number of rows in SoilParamsData")

  #Check LCTs
  allLCTs = c(forestLCTs, shrublandLCTs, grasslandLCTs, agricultureLCTs, rockLCTs, staticLCTs)
  if(sum(lctInput %in% allLCTs)!=length(lctInput)) stop("LCTs do not match. Please check.")
  if(nCells!=length(lctInput)) stop("The length of 'lctInput' has to be equal to number of grid cells")
  
  
  if(control$verbose) cat("Extracting SFI data")  
  #Find the plots that are inside the grid
  plotIndices = getGridIndex(xyplots, sg@grid, all.inside=FALSE) 
  #Discard plots that are outside the grid
  IDs = IDs[!is.na(plotIndices)] 
  xyplots = xyplots[!is.na(plotIndices),]
  plotIndices = plotIndices[!is.na(plotIndices)] 
  
  #Extract data
  x = SFItreeData[SFItreeData$ID %in% IDs, ]
  y = SFIshrubData[SFIshrubData$ID %in% IDs, ]
  if(!is.null(SFIcodes)) {
    x$Species = translateSpeciesCodes(x, SFIcodes)
    y$Species = translateSpeciesCodes(y, SFIcodes)
  }   
  lx = split(x, factor(x$ID, levels=IDs))
  ly = split(y, factor(y$ID, levels=IDs))
  setDefaults = control$setDefaults
  plotforestlist = Map(function(x,y, id) {
    extractSFIforest(x,y, id, SpParams=SpParams,patchsize=patchsize, setDefaults = setDefaults)
  }, lx, ly, IDs)  
  
  lctsPlots =lctInput[plotIndices]
  forLCTs = lctInput[lctInput %in% forestLCTs]
  if(sum(forLCTs %in% lctsPlots)!=length(forLCTs)) warning("Some forest land cover types do not have plots in them!")
  
  if(control$verbose) cat(" - Initializing soils")
  soillist = vector("list",nCells)
  for(i in 1:nCells) {
    if(!is.null(SoilParamData)) {
      soilParams = as.list(SoilParamData[i,])
      if(sum(is.na(soilParams))>0) soilParams = defaultSoilParams()
    } else {
      soilParams = defaultSoilParams()
    }
    soillist[[i]] = soil(soilParams)
  }
  
  
  if(control$verbose) cat(" - Initializing forest data")
  radius = control$initRadius  
  nInit = c(0,0,0)  
  nForest = 0
  nGrassland = 0
  nShrubland = 0
  
  findForestFor<-function(index, isShrubland=FALSE) {
    #Is index among plot indices? If yes, return corresponding plot
    pi = which(plotIndices==index)
    if(length(pi)==1) {
      forest= plotforestlist[[pi]]
      initType = 1
    }  else {      
      d = sqrt(rowSums(sweep(xyplots,2,coords[i,],"-")^2)) #Calculate distance to all plots 
      sameType = lctsPlots==lctInput[i] #Finds those plots that are of the same land cover type as the cell
      numPot = which((d<radius) & sameType) 
      if(length(numPot)>0) {
        forest = plotforestlist[[sample(numPot, 1)]]
        initType = 2
      } else {
        forest = plotforestlist[[which.min(d)]]
        initType = 3
      }
    }
    if(isShrubland) { #If shrubland remove all trees/saplings
      forest$treeData = forest$treeData[-(1:nrow(forest$treeData)),]
      forest$saplingData = forest$saplingData[-(1:nrow(forest$saplingData)),]
    }
    return(list(forest, initType))
  }
  forestlist = vector("list",nCells)    
  lct = vector("character",nCells)
  for(i in 1:nCells) {
    lcInputi  = lctInput[i]
    if(lcInputi %in% rockLCTs) {
      lct[i] = "Rock"
      forestlist[[i]] = NA
      soillist[[i]] = NA
    }
    else if(lcInputi %in% agricultureLCTs) {
      lct[i] = "Agriculture"
      forestlist[[i]] = emptyforest(patchsize=patchsize)
    }
    else if(lcInputi %in% grasslandLCTs) {
      lct[i] = "Wildland"
      forestlist[[i]] = emptyforest(patchsize=patchsize)
      nGrassland = nGrassland +1
    }
    else if(lcInputi %in% shrublandLCTs) {
      lct[i] = "Wildland"
      r = findForestFor(i, isShrubland = TRUE)
      forestlist[[i]] = r[[1]]
      nInit[r[[2]]] = nInit[r[[2]]]+1
      nShrubland = nShrubland +1
    } 
    else if(lcInputi %in% forestLCTs) {
      lct[i] = "Wildland"
      r = findForestFor(i, isShrubland=FALSE)
      forestlist[[i]] = r[[1]]
      nInit[r[[2]]] = nInit[r[[2]]]+1
      nForest = nForest +1
    } 
    else if(lcInputi %in% staticLCTs) {
      lct[i] = "Static"
      forestlist[[i]] = NA
      soillist[[i]] = NA
    }
  }
  if(control$verbose) cat(" - done\n")
  if(control$verbose) {
    cat(paste("Cells of wildland: ", sum(lct=="W"), ", agriculture: ", sum(lct=="A"),
              ", rocks: ", sum(lct=="R"), ", static: ", sum(lct=="S"),"\n", sep=""))
    cat(paste("Wildland cells with forest: ", nForest,", shrubland: ", nShrubland,
              ", grassland: ", nGrassland, "\n", sep=""))
    cat(paste("Wildland cells initialized from plot coordinates: ", nInit[1],"\n", sep=""))
    cat(paste("Wildland cells initialized from near plots of the same land cover type: ", nInit[2],"\n", sep=""))
    cat(paste("Wildland cells initialized from nearest plot: ", nInit[3],"\n", sep=""))
  }
  
  
  if(control$verbose==TRUE) cat(" - Queen neighbours")
  queenNeigh = cell2nb(grid@cells.dim[1],grid@cells.dim[2], type="queen")
  class(queenNeigh)<-"list"

  if(control$verbose==TRUE) cat(" - Water discharge order")
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
  if(control$verbose==TRUE) cat(" - done.\n")
  
  l = new("SpatialGridForest", 
           forestlist = forestlist,
           soillist = soillist,
           elevation = elevation,
           slope = slope,
           aspect = aspect,
           queenNeigh = queenNeigh,
           waterOrder = waterOrder,
           waterQ = waterQ,
           lct = lct,
           grid = grid, 
           bbox = landTopo@bbox, 
           proj4string = landTopo@proj4string)
   return(l)
}