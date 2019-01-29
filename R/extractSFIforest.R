extractSFIforest<-function(SFItreeData, SFIshrubData, ID, SpParams, 
                           SFIherbData = NULL, SFIcodes=NULL, 
                           patchsize = 10000, setDefaults=TRUE) {
  f = list()
  xid = SFItreeData[SFItreeData$ID==ID,]
  yid = SFIshrubData[SFIshrubData$ID==ID,]
  if(!is.null(SFIcodes)) {
    xid$Species = translateSpeciesCodes(xid, SFIcodes)
    yid$Species = translateSpeciesCodes(yid, SFIcodes)
  } 
  #Remove NA species
  xid = xid[!is.na(xid$Species),]
  yid = yid[!is.na(yid$Species),]
  
  f$ID = ID
  f$patchsize = patchsize
  f$treeData = data.frame(Species = xid$Species, N = round(xid$N), DBH = xid$DG, Height = xid$Ht)
  f$treeData$Z50 = rep(NA, nrow(f$treeData))
  f$treeData$Z95 = rep(NA, nrow(f$treeData))
  f$shrubData = data.frame(Species = yid$Species, Cover = as.numeric(yid$FCC), Height = yid$Ht)
  f$shrubData$Z50 =rep(NA, nrow(f$shrubData))
  f$shrubData$Z95 =rep(NA, nrow(f$shrubData))
  if(setDefaults) {
    f$treeData$Z95 = SpParams$Zmax[f$treeData$Species+1]
    f$treeData$Z95[is.na(f$treeData$Z95)] = 3000
    f$treeData$Z50 = f$treeData$Z95/4
    f$shrubData$Z95 = SpParams$Zmax[f$shrubData$Species+1]
    f$shrubData$Z95[is.na(f$shrubData$Z95)] = 1000
    f$shrubData$Z50 = f$shrubData$Z95/4
  }
  
  if(!is.null(SFIherbData)) {
    f$herbCover = SFIherbData$Cover[SFIherbData$ID==ID]
    f$herbHeight = SFIherbData$Height[SFIherbData$ID==ID]
    if(length(f$herbCover)>1) f$herbCover = mean(f$herbCover, na.rm=TRUE)
    if(length(f$herbHeight)>1) f$herbHeight = mean(f$herbHeight, na.rm=TRUE)
  } else {
    f$herbCover = NA
    f$herbHeight = NA
  }
  
  if(setDefaults) {
    usp = sort(unique(c(f$treeData$Species, f$saplingData$Species, f$shrubData$Species)))
    f$seedBank = data.frame(Species = usp, Abundance = rep(100, length(usp)))
  } else {
    f$seedBank = data.frame(Species = numeric(0), Abundance = numeric(0))
  }

  class(f)<-c("forest","list")
  return(f)
}