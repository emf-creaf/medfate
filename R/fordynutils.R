.summarizeCohorts<-function(step, 
                            treeTableYear, shrubTableYear,
                            deadTreeTableYear, deadShrubTableYear,
                            cutTreeTableYear, cutShrubTableYear) {
  treecohnames = unique(c(treeTableYear$Cohort,deadTreeTableYear$Cohort))
  if(!is.null(cutTreeTableYear)) treecohnames = unique(c(treecohnames, cutTreeTableYear$Cohort))
  shrubcohnames = unique(c(shrubTableYear$Cohort,deadShrubTableYear$Cohort))
  if(!is.null(cutShrubTableYear)) shrubcohnames = unique(c(shrubcohnames, cutShrubTableYear$Cohort))
  cohnames = c(treecohnames, shrubcohnames)
  isTree = c(rep(TRUE, length(treecohnames)), rep(FALSE, length(shrubcohnames)))
  cohortSummary = data.frame("Step" = rep(step, length(cohnames)),
                             "Species" = NA,
                             "Name" = NA,
                             "Cohort" = cohnames,
                             "TreeDensityLive" = 0,
                             "TreeBasalAreaLive"= 0,
                             "ShrubCoverLive"= 0,
                             "BasalAreaDead" = 0,
                             "ShrubCoverDead" = 0,
                             "BasalAreaCut" = 0,
                             "ShrubCoverCut" = 0)
  cohortSummary$ShrubCoverLive[isTree] = NA
  cohortSummary$ShrubCoverDead[isTree] = NA
  cohortSummary$ShrubCoverCut[isTree] = NA
  cohortSummary$TreeDensityLive[!isTree] = NA
  cohortSummary$TreeBasalAreaLive[!isTree] = NA
  cohortSummary$BasalAreaDead[!isTree] = NA
  cohortSummary$BasalAreaCut[!isTree] = NA
  
  ba_live = .treeBasalArea(treeTableYear$N, treeTableYear$DBH)
  for(i in 1:nrow(treeTableYear)) {
    icoh = which(cohortSummary$Cohort==treeTableYear$Cohort[i])
    cohortSummary$Cohort[icoh] = treeTableYear$Cohort[i]
    cohortSummary$Species[icoh] = treeTableYear$Species[i]
    cohortSummary$Name[icoh] = treeTableYear$Name[i]
    cohortSummary$TreeDensityLive[icoh] = treeTableYear$N[i]
    cohortSummary$TreeBasalAreaLive[icoh] = ba_live[i]
  }
  for(i in 1:nrow(shrubTableYear)) {
    icoh = which(cohortSummary$Cohort==shrubTableYear$Cohort[i])
    cohortSummary$Cohort[icoh] = shrubTableYear$Cohort[i]
    cohortSummary$Species[icoh] = shrubTableYear$Species[i]
    cohortSummary$Name[icoh] = shrubTableYear$Name[i]
    cohortSummary$ShrubCoverLive[icoh] = shrubTableYear$Cover[i]
  }
  ba_dead = .treeBasalArea(deadTreeTableYear$N, deadTreeTableYear$DBH)
  for(i in 1:nrow(deadTreeTableYear)) {
    icoh = which(cohortSummary$Cohort==deadTreeTableYear$Cohort[i])
    cohortSummary$Cohort[icoh] = deadTreeTableYear$Cohort[i]
    cohortSummary$Species[icoh] = deadTreeTableYear$Species[i]
    cohortSummary$Name[icoh] = deadTreeTableYear$Name[i]
    cohortSummary$BasalAreaDead[icoh] = ba_dead[i]
  }
  for(i in 1:nrow(deadShrubTableYear)) {
    icoh = which(cohortSummary$Cohort==deadShrubTableYear$Cohort[i])
    cohortSummary$Cohort[icoh] = deadShrubTableYear$Cohort[i]
    cohortSummary$Species[icoh] = deadShrubTableYear$Species[i]
    cohortSummary$Name[icoh] = deadShrubTableYear$Name[i]
    cohortSummary$ShrubCoverDead[icoh] = deadShrubTableYear$Cover[i]
  }
  if(!is.null(cutTreeTableYear)) {
    ba_cut = .treeBasalArea(cutTreeTableYear$N, cutTreeTableYear$DBH)
    for(i in 1:nrow(cutTreeTableYear)) {
      icoh = which(cohortSummary$Cohort==cutTreeTableYear$Cohort[i])
      cohortSummary$Cohort[icoh] = cutTreeTableYear$Cohort[i]
      cohortSummary$Species[icoh] = cutTreeTableYear$Species[i]
      cohortSummary$Name[icoh] = cutTreeTableYear$Name[i]
      cohortSummary$BasalAreaCut[icoh] = ba_cut[i]
    }
  }
  if(!is.null(cutShrubTableYear)) {
    for(i in 1:nrow(cutShrubTableYear)) {
      icoh = which(cohortSummary$Cohort==cutShrubTableYear$Cohort[i])
      cohortSummary$Cohort[icoh] = cutShrubTableYear$Cohort[i]
      cohortSummary$Species[icoh] = cutShrubTableYear$Species[i]
      cohortSummary$Name[icoh] = cutShrubTableYear$Name[i]
      cohortSummary$ShrubCoverCut[icoh] = cutShrubTableYear$Cover[i]
    }
  }
  return(cohortSummary)
}
.summarizeSpecies<-function(step, cohSum, x, SpParams) {
  nl_sp = tapply(cohSum$TreeDensityLive, cohSum$Species, sum, na.rm=FALSE)
  bal_sp = tapply(cohSum$TreeBasalAreaLive, cohSum$Species, sum, na.rm=FALSE)
  bad_sp = tapply(cohSum$BasalAreaDead, cohSum$Species, sum, na.rm=FALSE)
  bac_sp = tapply(cohSum$BasalAreaCut, cohSum$Species, sum, na.rm=FALSE)
  shl_sp = tapply(cohSum$ShrubCoverLive, cohSum$Species, sum, na.rm=FALSE)
  shd_sp = tapply(cohSum$ShrubCoverDead, cohSum$Species, sum, na.rm=FALSE)
  shc_sp = tapply(cohSum$ShrubCoverCut, cohSum$Species, sum, na.rm=FALSE)
  spSumYear <-data.frame("Step" = rep(step, length(nl_sp)),
                         "Species" = as.integer(names(nl_sp)),
                         "Name" = species_characterParameter(as.integer(names(nl_sp)), SpParams, "Name"),
                         "NumCohorts" = as.numeric(table(cohSum$Species)),
                         "TreeDensityLive"= as.numeric(nl_sp),
                         "TreeBasalAreaLive"= as.numeric(bal_sp),
                         "ShrubCoverLive"= as.numeric(shl_sp),
                         "BasalAreaDead" = as.numeric(bad_sp),
                         "ShrubCoverDead" = as.numeric(shd_sp),
                         "BasalAreaCut" = as.numeric(bac_sp),
                         "ShrubCoverCut" = as.numeric(shc_sp))
}
.summarizeStand<-function(step, cohSum, x) {
  isTree = !is.na(x$above$DBH)
  HB = .hartBeckingIndex(x$above$N[isTree], x$above$H[isTree], x$above$DBH[isTree])
  domH = .dominantTreeHeight(x$above$N[isTree], x$above$H[isTree], x$above$DBH[isTree])
  domDBH = .dominantTreeDiameter(x$above$N[isTree], x$above$DBH[isTree])
  qmDBH = .quadraticMeanTreeDiameter(x$above$N[isTree], x$above$DBH[isTree])
  standSumYear = data.frame("Step" = step,
                            "NumTreeSpecies" = length(unique(x$above$SP[isTree])),
                            "NumTreeCohorts" = sum(isTree),
                            "NumShrubSpecies" = length(unique(x$above$SP[!isTree])),
                            "NumShrubCohorts" = sum(!isTree),
                            "TreeDensityLive"= sum(cohSum$TreeDensityLive, na.rm = TRUE),
                            "TreeBasalAreaLive"= sum(cohSum$TreeBasalAreaLive, na.rm = TRUE),
                            "DominantTreeHeight"= domH,
                            "DominantTreeDiameter"= domDBH,
                            "QuadraticMeanTreeDiameter"= qmDBH,
                            "HartBeckingIndex"= HB,
                            "ShrubCoverLive"= sum(cohSum$ShrubCoverLive, na.rm=TRUE),
                            "BasalAreaDead" = sum(cohSum$BasalAreaDead, na.rm=TRUE),
                            "ShrubCoverDead" = sum(cohSum$ShrubCoverDead, na.rm=TRUE),
                            "BasalAreaCut" = sum(cohSum$BasalAreaCut, na.rm=TRUE),
                            "ShrubCoverCut" = sum(cohSum$ShrubCoverCut, na.rm=TRUE))
}

.createTreeTable<-function(step, year, x) {
  isTree = !is.na(x$above$DBH)
  range = 1:sum(isTree)
  tt<-data.frame(Step = rep(step, length(range)), 
                 Year = rep(year, length(range)),
                 Cohort = row.names(x$cohorts)[range],
                 Species = x$above$SP[range],
                 Name = x$cohorts$Name[range],
                 N = x$above$N[range],
                 DBH = x$above$DBH[range],
                 Height = x$above$H[range],
                 Z50 = x$below$Z50[range],
                 Z95 = x$below$Z95[range])
  tt = tt[tt$N>0,, drop=FALSE]
  return(tt)
}
.createDeadTreeTable<-function(step, year, x) {
  isTree = !is.na(x$above$DBH)
  range = numeric(0)
  if(sum(isTree)>0) range = 1:sum(isTree)
  dtt<-data.frame(Step = rep(step, length(range)), 
                  Year = rep(year, length(range)),
                  Cohort = row.names(x$cohorts)[range],
                  Species = x$above$SP[range],
                  Name = x$cohorts$Name[range],
                  N = x$internalMortality$N_dead[range],
                  DBH = x$above$DBH[range],
                  Height = x$above$H[range],
                  Z50 = x$below$Z50[range],
                  Z95 = x$below$Z95[range])
  dtt = dtt[dtt$N>0,, drop = FALSE]
  return(dtt)
}
.createCutTreeTable<-function(step, year, x, N_cut) {
  range = 1:length(N_cut)
  tt<-data.frame(Step = rep(step, length(N_cut)), 
                 Year = rep(year, length(N_cut)),
                 Cohort = row.names(x$cohorts)[range],
                 Species = x$above$SP[range],
                 Name = x$cohorts$Name[range],
                 N = N_cut,
                 DBH = x$above$DBH[range],
                 Height = x$above$H[range],
                 Z50 = x$below$Z50[range],
                 Z95 = x$below$Z95[range])
  tt = tt[tt$N>0,, drop=FALSE]
  return(tt)
}
.createShrubTable<-function(step, year, x) {
  isShrub = !is.na(x$above$Cover)
  nt = sum(!isShrub)
  numCohorts = length(isShrub)
  range = numeric(0)
  if(numCohorts>nt) range = (nt+1):numCohorts
  st<-data.frame(Step = rep(step, length(range)), 
                 Year = rep(year, length(range)),
                 Cohort = row.names(x$cohorts)[range],
                 Species = x$above$SP[range],
                 Name = x$cohorts$Name[range],
                 Cover = x$above$Cover[range],
                 Height = x$above$H[range],
                 Z50 = x$below$Z50[range],
                 Z95 = x$below$Z95[range])
  st = st[st$Cover>0,, drop = FALSE]
  return(st)
}
.createDeadShrubTable<-function(step, year, x) {
  isShrub = !is.na(x$above$Cover)
  nt = sum(!isShrub)
  numCohorts = length(isShrub)
  range = numeric(0)
  if(numCohorts>nt) range = (nt+1):numCohorts
  dst<-data.frame(Step = rep(step, length(range)), 
                  Year = rep(year, length(range)),
                  Cohort = row.names(x$cohorts)[range],
                  Species = x$above$SP[range],
                  Name = x$cohorts$Name[range],
                  Cover = x$internalMortality$Cover_dead[range],
                  Height = x$above$H[range],
                  Z50 = x$below$Z50[range],
                  Z95 = x$below$Z95[range])
  dst = dst[dst$Cover>0,,drop=FALSE]
  return(dst)
}
.createCutShrubTable<-function(step, year, x, Cover_cut) {
  isShrub = !is.na(x$above$Cover)
  nt = sum(!isShrub)
  numCohorts = length(isShrub)
  range = numeric(0)
  if(numCohorts>nt) range = (nt+1):numCohorts
  st<-data.frame(Step = rep(step, length(range)), 
                 Year = rep(year, length(range)),
                 Cohort = row.names(x$cohorts)[range],
                 Species = x$above$SP[range],
                 Name = x$cohorts$Name[range],
                 Cover = Cover_cut,
                 Height = x$above$H[range],
                 Z50 = x$below$Z50[range],
                 Z95 = x$below$Z95[range])
  st = st[st$Cover>0,, drop = FALSE]
  return(st)
}