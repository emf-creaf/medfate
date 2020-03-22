summary.forest<-function(object, SpParams, mode = "MED", detailed = FALSE, ...) {
  summaryBasalArea<-function(x, SpParams) {
    treeSpecies = SpParams$SpIndex[SpParams$GrowthForm=="Tree"]    
    if(nrow(x$treeData)>0) {
      ft = factor(x$treeData$Species, levels=treeSpecies)
      tba = .treeBasalArea(x$treeData$N, x$treeData$DBH)
      tBAsp = tapply(tba,ft, FUN=sum)
      tBAsp[is.na(tBAsp)] = 0
      BA_trees = sum(tba)
    } else {
      tBAsp = rep(0, length(treeSpecies))
      BA_trees = 0
    }
    BAsp = tBAsp
    BA = BA_trees
    names(BAsp)<-paste("BA",treeSpecies,sep="_")
    if(detailed) return(c(BA = BA, BAsp))
    else return(c(BA = BA))
  }
  summaryNumber<-function(x, SpParams) {
    treeSpecies = SpParams$SpIndex[SpParams$GrowthForm=="Tree"]
    if(nrow(x$treeData)>0) {
      ft = factor(x$treeData$Species, levels=treeSpecies)
      Nt = sum(x$treeData$N)
      tNsp = tapply(x$treeData$N,ft, FUN=sum)
      tNsp[is.na(tNsp)] = 0
    } else {
      tNsp = rep(0, length(treeSpecies))
      Nt = 0
    }
    Nsp = tNsp
    names(Nsp)<-paste("N",treeSpecies,sep="_")
    if(detailed) return(c(N = Nt, Nsp))
    else return(c(N = Nt))
  }
  summaryLAI<-function(x, SpParams) {
    ntree = nrow(x$treeData)
    nshrub = nrow(x$treeData)
    LAIc = plant_LAI(x, SpParams, mode= mode)
    LAIt = 0
    LAIsh = 0
    if(ntree>0) LAIt = sum(LAIc[1:ntree], na.rm=T)
    if(nshrub>0) LAIsh = sum(LAIc[(1+ntree):(1+ntree+nshrub)], na.rm=T)
    LAIsp = species_LAI(x, SpParams, mode= mode)
    if(detailed) return(c(LAI= (LAIt+LAIsh), LAI_trees = LAIt,
             LAI_shrubs = LAIsh, LAIsp))
    else return(c(LAI= (LAIt+LAIsh), LAI_trees = LAIt,
                  LAI_shrubs = LAIsh))
  }
  summaryFuel<-function(x, SpParams) {
    ntree = nrow(x$treeData)
    nshrub = nrow(x$treeData)
    Fuelc = plant_fuel(x, SpParams, mode= mode)
    Fuelt = 0
    Fuelsh = 0
    if(ntree>0) Fuelt = sum(Fuelc[1:ntree], na.rm=T)
    if(nshrub>0) Fuelsh = sum(Fuelc[(1+ntree):(1+ntree+nshrub)], na.rm=T)
    Fuelsp = species_fuel(x, SpParams, mode= mode)
    if(detailed) return(c(Fuel= (Fuelt+Fuelsh), Fuel_trees = Fuelt, 
                          Fuel_shrubs = Fuelsh, Fuelsp))
    else return(c(Fuel= (Fuelt+Fuelsh), Fuel_trees = Fuelt,
                  Fuel_shrubs = Fuelsh))
  }
  s = c(summaryNumber(object, SpParams),
        summaryBasalArea(object, SpParams),
        summaryLAI(object,SpParams),
        summaryFuel(object,SpParams))
  s["Phytovolume"] = sum(plant_phytovolume(object, SpParams),na.rm=TRUE)
  class(s)<-c("summary.forest","list")
  return(s)
}
print.summary.forest<-function(x, digits=getOption("digits"),...) {
  cat("Forest summary:\n\n")
  cat(paste("  Tree density (ind/ha):", x["N"],"\n"))
  cat(paste("  BA (m2/ha):", round(x["BA"],digits),"\n"))
  cat(paste("  Shrub crown phytovolume (m3/m2):", round(x["Phytovolume"],digits),"\n"))
  cat(paste("  LAI (m2/m2) total:", round(x["LAI"], digits)," trees:", round(x["LAI_trees"], digits),
            " shrubs:", round(x["LAI_shrubs"], digits),"\n"))
  cat(paste("  Live fine fuel (kg/m2) total:", round(x["Fuel"], digits)," trees:", round(x["Fuel_trees"], digits),
            " shrubs:", round(x["Fuel_shrubs"], digits),"\n"))
}