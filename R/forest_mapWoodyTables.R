forest_mapTreeTable<-function(x, mapping, SpParams, plot.size = NULL) {
  n = nrow(x)
  treeData = data.frame(
    Species = rep(NA, n),
    N = rep(NA, n),
    Height = rep(NA, n),
    DBH = rep(NA, n),
    Z50 = rep(NA, n),
    Z95 = rep(NA, n))
  
  if("Height" %in% names(mapping)) {
    treeData$Height = x[[mapping[["Height"]]]]
  }
  if("DBH" %in% names(mapping)) {
    treeData$DBH = x[[mapping[["DBH"]]]]
  }
  if("N" %in% names(mapping)) {
    treeData$N = x[[mapping[["N"]]]]
  } else {
    treeData$N = 1
  }
  if("Z50" %in% names(mapping)) {
    treeData$Z50 = x[[mapping[["Z50"]]]]
  }
  if("Z95" %in% names(mapping)) {
    treeData$Z95 = x[[mapping[["Z95"]]]]
  }
  if("Species" %in% names(mapping)) {
    treeData$Species = x[[mapping[["Species"]]]]
  }
  if("plot.size" %in% names(mapping)) {
    plot.size = x[[mapping[["plot.size"]]]]
  } 
  if(!is.null(plot.size)) {
    treeData$N = treeData$N*(10000/plot.size)
  }
  if("Species.name" %in% names(mapping)) {
    Species.name = x[[mapping[["Species.name"]]]]
    for(i in 1:n) {
      indices = which(SpParams$Name==Species.name[i])
      if(length(indices)>0) {
        treeData$Species[i] = SpParams$SpIndex[indices]
      }
    }
  }
  return(treeData)
}

forest_mapShrubTable<-function(y, mapping, SpParams) {
  n = nrow(y)
  shrubData = data.frame(
    Species = rep(NA, n),
    Height = rep(NA, n),
    Cover = rep(NA, n),
    Z50 = rep(NA, n),
    Z95 = rep(NA, n))
  
  if("Height" %in% names(mapping)) {
    shrubData$Height = y[[mapping[["Height"]]]]
  }
  if("Cover" %in% names(mapping)) {
    shrubData$Cover = y[[mapping[["Cover"]]]]
  }
  if("Z50" %in% names(mapping)) {
    shrubData$Z50 = y[[mapping[["Z50"]]]]
  }
  if("Z95" %in% names(mapping)) {
    shrubData$Z95 = y[[mapping[["Z95"]]]]
  }
  if("Species" %in% names(mapping)) {
    shrubData$Species = y[[mapping[["Species"]]]]
  }
  if("Species.name" %in% names(mapping)) {
    Species.name = y[[mapping[["Species.name"]]]]
    for(i in 1:n) {
      indices = which(SpParams$Name==Species.name[i])
      if(length(indices)>0) {
        shrubData$Species[i] = SpParams$SpIndex[indices]
      }
    }
  }
  return(shrubData)
}

forest_mapWoodyTables<-function(x, y, mapping, SpParams, plot.size=NULL) {
  f = emptyforest()
  f$treeData = forest_mapTreeTable(x, mapping = mapping,  SpParams=SpParams, plot.size = plot.size)
  f$shrubData = forest_mapShrubTable(y, mapping = mapping,  SpParams=SpParams)
  return(f)
} 