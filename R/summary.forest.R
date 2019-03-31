summary.forest<-function(object, SpParams, detailed = FALSE, ...) {
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
    species = SpParams$SpIndex
    if(nrow(x$treeData)>0) {
      ft = factor(x$treeData$Species, levels=species)
      tLAI = .treeLAI(SP=x$treeData$Species, N=x$treeData$N, dbh=x$treeData$DBH, SpParams=SpParams)
      tLAI[is.na(tLAI)] = 0
      tLAIsp = tapply(tLAI,ft, FUN=sum)
      tLAIsp[is.na(tLAIsp)] = 0
      LAIt = sum(tLAIsp)
    } else {
      tLAIsp = rep(0, length(species))
      LAIt = 0
    }
    if(nrow(x$shrubData)>0) {
      fsh = factor(x$shrubData$Species, levels=species)
      shLAI = .shrubLAI(SP=x$shrubData$Species, Cover=x$shrubData$Cover, H=x$shrubData$Height,SpParams=SpParams)
      shLAI[is.na(shLAI)] = 0
      shLAIsp = tapply(shLAI,fsh, FUN=sum)
      shLAIsp[is.na(shLAIsp)] = 0
      LAIsh = sum(shLAI)
    } else {
      shLAIsp = rep(0, length(species))
      LAIsh = 0
    }
    LAIsp = tLAIsp+shLAIsp
    names(LAIsp)<-paste("LAI",species,sep="_")
    if(detailed) return(c(LAI= (LAIt+LAIsh), LAI_trees = LAIt,
             LAI_shrubs = LAIsh, LAIsp))
    else return(c(LAI= (LAIt+LAIsh), LAI_trees = LAIt,
                  LAI_shrubs = LAIsh))
  }
  summaryFuel<-function(x, SpParams) {
    species = SpParams$SpIndex
    if(nrow(x$treeData)>0) {
      ft = factor(x$treeData$Species, levels=species)
      tFuel = .treeFuel(SP=x$treeData$Species, N=x$treeData$N, dbh=x$treeData$DBH, SpParams=SpParams)
      tFuel[is.na(tFuel)] = 0
      tFuelsp = tapply(tFuel,ft, FUN=sum)
      tFuelsp[is.na(tFuelsp)] = 0
      Fuelt = sum(tFuelsp)
    } else {
      tFuelsp = rep(0, length(species))
      Fuelt = 0
    }
    if(nrow(x$shrubData)>0) {
      fsh = factor(x$shrubData$Species, levels=species)
      shFuel = .shrubFuel(SP=x$shrubData$Species, Cover=x$shrubData$Cover, H=x$shrubData$Height, 
                          CR = .shrubCrownRatio(x$shrubData$Species, SpParams),
                          SpParams=SpParams)
      shFuel[is.na(shFuel)] = 0
      shFuelsp = tapply(shFuel,fsh, FUN=sum)
      shFuelsp[is.na(shFuelsp)] = 0
      Fuelsh = sum(shFuel)
    } else {
      shFuelsp = rep(0, length(species))
      Fuelsh = 0
    }
    Fuelsp = tFuelsp+shFuelsp
    names(Fuelsp)<-paste("Fuel",species,sep="_")
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
  class(s)<-c("summary.forest","numeric")
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