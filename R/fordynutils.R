.summarizeCohorts<-function(step, x) {
  isTree = !is.na(x$above$DBH)
  cohortSummary = data.frame("Step" = rep(step, nrow(x$cohorts)),
                             "Species" = x$cohorts$SP,
                             "Name" = x$cohorts$Name,
                             "Cohort" = row.names(x$cohorts),
                             "LeafAreaIndex" = x$above$LAI_live,
                             "TreeDensityLive" = x$above$N,
                             "TreeBasalAreaLive"= x$above$N*pi*(x$above$DBH/200)^2,
                             "ShrubCoverLive"= x$above$Cover)
  cohortSummary$TreeDensityLive[!isTree] = NA
  return(cohortSummary)
}
.summarizeSpecies<-function(step, cohSum, x) {
  lai_sp = tapply(cohSum$LeafAreaIndex, cohSum$Species, sum, na.rm=TRUE)
  nl_sp = tapply(cohSum$TreeDensityLive, cohSum$Species, sum, na.rm=TRUE)
  bal_sp = tapply(cohSum$TreeBasalAreaLive, cohSum$Species, sum, na.rm=TRUE)
  shl_sp = tapply(cohSum$ShrubCoverLive, cohSum$Species, sum, na.rm=TRUE)
  mh_sp = tapply(x$above$H, x$above$SP, max, na.rm=TRUE)
  spSumYear <-data.frame("Step" = rep(step, length(lai_sp)),
                         "Species" = as.numeric(names(lai_sp)),
                         "LeafAreaIndex" = as.numeric(lai_sp),
                         "TreeDensityLive"= as.numeric(nl_sp),
                         "TreeBasalAreaLive"= as.numeric(bal_sp),
                         "ShrubCoverLive"= as.numeric(shl_sp),
                         "MaxHeight"= as.numeric(mh_sp))
}
.summarizeStand<-function(step, cohSum, x) {
  maxH = 0
  if(nrow(x$above)>0) maxH = max(x$above$H, na.rm=TRUE)
  standSumYear = data.frame("Step" = step,
                            "LeafAreaIndex" = sum(cohSum$LeafAreaIndex, na.rm = TRUE),
                            "TreeDensityLive"= sum(cohSum$TreeDensityLive, na.rm = TRUE),
                            "TreeBasalAreaLive"= sum(cohSum$TreeBasalAreaLive, na.rm = TRUE),
                            "ShrubCoverLive"= sum(cohSum$ShrubCoverLive, na.rm=TRUE),
                            "MaxHeight"= maxH)
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
  range = 1:sum(isTree)
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