library(tidyverse)
library(openxlsx)
DB_path = "/home/miquel/OneDrive/Professional/MedfateWorks/MedfateSpeciesParametrization/"

genfun <- function(x) {
  s = strsplit(x,split = " ")
  for(i in 1:length(s)) s[i] = s[[i]][[1]]
  return(unlist(s))
}
famfun <- function(gen) {
  df = taxize::get_gbifid_(gen, messages = FALSE)[[1]]
  if(sum(df$kingdom=="Plantae")>0) return(df$family[df$kingdom=="Plantae"][1])
  return(NA)
}
TFM_fun<-function(species, values, sp_fun = "mean") {
  trait_means = tapply(values, species, sp_fun, na.rm=TRUE)
  gen = genfun(names(trait_means))
  gen_unique = unique(gen)
  fam = rep(NA, length(gen))
  pb = txtProgressBar(0, length(gen_unique), style=3)
  for(i in 1:length(gen_unique)) {
    setTxtProgressBar(pb, i)
    fam[gen==gen_unique[i]] = famfun(gen_unique[i])
  }
  trait_fam<-tapply(trait_means, fam, mean, na.rm=TRUE)
  return(trait_fam[!is.na(trait_fam)])
}

# Leaf density
TRY_LeafDensity  <-readRDS(file.path(DB_path,"TraitDatabases/TRY/TRY_traits/TRY_48.rds"))
TRY_LeafDensity<- TRY_LeafDensity[TRY_LeafDensity$ErrorRisk < 3,, drop= FALSE] 
leaf_dens_fam<-TFM_fun(TRY_LeafDensity$AccSpeciesName, TRY_LeafDensity$StdValue)


# Wood density
TRY_WoodDensity  <-readRDS(file.path(DB_path,"TraitDatabases/TRY/TRY_traits/TRY_4.rds"))
TRY_WoodDensity <- TRY_WoodDensity[TRY_WoodDensity$ErrorRisk<3, c("AccSpeciesName", "StdValue")]
names(TRY_WoodDensity) <- c("Species", "WoodDensity")
#Global wood density database
# Zanne, A.E., Lopez-Gonzalez, G.*, Coomes, D.A., Ilic, J., Jansen, S., Lewis, S.L., Miller, R.B., Swenson, N.G., Wiemann, M.C., and Chave, J. 2009. Global wood density database. 
# Dryad. Identifier: http://hdl.handle.net/10255/dryad.235. 
gwdd <- readxl::read_excel(file.path(DB_path,"TraitDatabases/GlobalWoodDensityDatabase/GlobalWoodDensityDatabase.xls"), sheet=2)
gwdd_WoodDensity<-gwdd[,c("Binomial", "Wood density (g/cm^3), oven dry mass/fresh volume")]
names(gwdd_WoodDensity)<-c("Species", "WoodDensity")
## HydraTRY
fn_hydratry <- file.path(DB_path,"TraitDatabases/HydraTRY/HydraTRY.xlsx")
hydratry <- openxlsx::read.xlsx(fn_hydratry, rowNames = FALSE) 
names(hydratry)[names(hydratry)=="WD"] = "WoodDensity"
##BROT2
fn_brot2 <- file.path(DB_path, "TraitDatabases/BROT2/BROT2_dat.csv")
brot2 <- read.csv2(fn_brot2)
brot2_WoodDensity <- brot2 %>% filter(Trait %in% c("StemDensity")) %>% 
  select(c("Taxon", "Data"))
brot2_WoodDensity$Data <- as.numeric(brot2_WoodDensity$Data)
names(brot2_WoodDensity) <- c("Species", "WoodDensity")

WoodDensity_data<-rbind(hydratry[,c("Species","WoodDensity")],
                        TRY_WoodDensity,
                        brot2_WoodDensity,
                        gwdd_WoodDensity)
wood_dens_fam<-TFM_fun(WoodDensity_data$Species, WoodDensity_data$WoodDensity)

# Leaf pressure-volume curve parameters (LeafPI0, LeafEPS and LeafAF)
fn_bartlett <- file.path(DB_path, "TraitDatabases/Bartlett_et_al_2012/Bartlett_2012_ELE_data.xlsx")
bartlett <- openxlsx::read.xlsx(fn_bartlett, rowNames = FALSE) 
names(bartlett)[names(bartlett)=="πo.(MPa)"] = "LeafPI0"
names(bartlett)[names(bartlett)=="ε.(MPa)"] = "LeafEPS"
TRY_LeafPI0  <-readRDS(file.path(DB_path, "TraitDatabases/TRY/TRY_traits/TRY_188.rds"))
TRY_LeafPI0 <- TRY_LeafPI0[TRY_LeafPI0$ErrorRisk<3, c("AccSpeciesName", "StdValue")]
names(TRY_LeafPI0) <- c("Species", "LeafPI0")
TRY_LeafPI0$LeafPI0 <- -1*TRY_LeafPI0$LeafPI0 
TRY_LeafEPS  <-readRDS(file.path(DB_path, "TraitDatabases/TRY/TRY_traits/TRY_190.rds"))
TRY_LeafEPS$StdValue <- as.numeric(TRY_LeafEPS$OrigValueStr) 
TRY_LeafEPS <- TRY_LeafEPS[, c("AccSpeciesName", "StdValue")]
names(TRY_LeafEPS) <- c("Species", "LeafEPS")
LeafPI0_data = rbind(TRY_LeafPI0[, c("Species", "LeafPI0")],
                     bartlett[,c("Species", "LeafPI0")])
LeafEPS_data = rbind(TRY_LeafEPS[, c("Species", "LeafEPS")],
                     bartlett[,c("Species", "LeafEPS")])
leaf_PI0_fam = TFM_fun(LeafPI0_data$Species, LeafPI0_data$LeafPI0)
leaf_EPS_fam = TFM_fun(LeafEPS_data$Species, LeafEPS_data$LeafEPS)
leaf_AF_fam = TFM_fun(bartlett$Species, bartlett$af)

# Maximum stomatal conductance
fn_hoshika <- file.path(DB_path, "TraitDatabases/Hoshika_et_al_2018/Hoshika_et_al_2018.xlsx")
hoshika <- openxlsx::read.xlsx(fn_hoshika) 
hoshika$Species <- str_replace(string = hoshika$Species,pattern = "_", replacement = " ")
Gswmax_fam = TFM_fun(hoshika$Species, hoshika$`gmax.(mol.m-2.s-1)`)

# Minimum stomatal conductance
fn_duursma <- file.path(DB_path, "TraitDatabases/Duursma_et_al_2018/gmindatabase.csv")
duursma <- read.csv2(fn_duursma, sep=",")
duursma$gmin <- as.numeric(duursma$gmin)/1000
Gswmin_fam = TFM_fun(duursma$species, duursma$gmin)

# N leaf
TRY_Nleaf <-readRDS(file.path(DB_path,"TraitDatabases/TRY/TRY_traits/TRY_14.rds"))
TRY_Nleaf <- TRY_Nleaf[TRY_Nleaf$ErrorRisk<3, c("AccSpeciesName", "StdValue")]
Nleaf_fam = TFM_fun(TRY_Nleaf$AccSpeciesName, TRY_Nleaf$StdValue)

# N sapwood
TRY_Nsapwood <-readRDS(file.path(DB_path,"TraitDatabases/TRY/TRY_traits/TRY_1229.rds"))
TRY_Nsapwood <- TRY_Nsapwood[TRY_Nsapwood$ErrorRisk<3, c("AccSpeciesName", "StdValue")]
Nsapwood_fam = TFM_fun(TRY_Nsapwood$AccSpeciesName, TRY_Nsapwood$StdValue)

# N fineroot
TRY_Nfineroot <-readRDS(file.path(DB_path,"TraitDatabases/TRY/TRY_traits/TRY_475.rds"))
TRY_Nfineroot <- TRY_Nfineroot[TRY_Nfineroot$ErrorRisk<3, c("AccSpeciesName", "StdValue")]
Nfineroot_fam = TFM_fun(TRY_Nfineroot$AccSpeciesName, TRY_Nfineroot$StdValue)

# Wood C
TRY_WoodC <-readRDS(file.path(DB_path,"TraitDatabases/TRY/TRY_traits/TRY_407.rds"))
TRY_WoodC <- TRY_WoodC[TRY_WoodC$ErrorRisk<3, c("AccSpeciesName", "StdValue")]
WoodC_fam = TFM_fun(TRY_WoodC$AccSpeciesName, TRY_WoodC$StdValue/1000)


# Kmax_stemxylem
fn_hydratry <- file.path(DB_path,"TraitDatabases/HydraTRY/HydraTRY.xlsx")
hydratry <- openxlsx::read.xlsx(fn_hydratry, rowNames = FALSE) 
Ks_fam = TFM_fun(hydratry$Species, hydratry$Ks)

# Psi_Critic = P50
fn_martin_stpaul <- file.path(DB_path,"TraitDatabases/MartinStPaul_et_al_2017/DataBase.xlsx")
martin_stpaul <- read.xlsx(fn_martin_stpaul, sheet="ALL")
martin_stpaul<-martin_stpaul[,c("Species", "P50", "P12", "Ptlp")]
fn_choat <-file.path(DB_path,"TraitDatabases/Choat_et_al_2012/nature11688-s2.xls")
choat <- readxl::read_excel(fn_choat)
choat <-choat[-1, c("Species", "P50", "P88")] # Remove unit row
P50_data = rbind(martin_stpaul[, c("Species", "P50")],
                 choat[, c("Species", "P50")],
                 hydratry[, c("Species", "P50")])
P50_data$P50 <- as.numeric(P50_data$P50)
TRY_VC <-readRDS(file.path(DB_path,"TraitDatabases/TRY/TRY_traits/TRY_719.rds"))
VC_str <- c("Reference P50",
            "Xylem water potential at which 50% of conductivity is lost (P50)",
            "Stem P50")
TRY_P50 <- TRY_VC[TRY_VC$DataName %in% VC_str, c("AccSpeciesName", "OrigValueStr")]
names(TRY_P50)<-c("Species", "P50")
TRY_P50$P50 <- as.numeric(TRY_P50$P50)
TRY_P50 <- TRY_P50[!is.na(TRY_P50$P50),]
TRY_P50$P50[TRY_P50$P50>0] <- -TRY_P50$P50[TRY_P50$P50>0]
P50_data = rbind(P50_data, TRY_P50)
P50_fam = TFM_fun(P50_data$Species, P50_data$P50)

# Al2As
TRY_Al2As  <-readRDS(file.path(DB_path,"TraitDatabases/TRY/TRY_traits/TRY_171.rds"))
TRY_Al2As$StdValue[TRY_Al2As$OriglName=="Huber value"]<-1/as.numeric(TRY_Al2As$OrigValueStr[TRY_Al2As$OriglName=="Huber value"])
TRY_Al2As$StdValue[TRY_Al2As$OriglName=="leaf.area.per.sapwood.area"]<-TRY_Al2As$StdValue[TRY_Al2As$OriglName=="leaf.area.per.sapwood.area"]*100
TRY_Al2As$StdValue[TRY_Al2As$OriglName=="Sapwood: leaf area ratio"]<-10000/as.numeric(TRY_Al2As$OrigValueStr[TRY_Al2As$OriglName=="Sapwood: leaf area ratio"])
TRY_Al2As$StdValue[TRY_Al2As$OriglName=="The ratio of leaf area attached per unit sapwood cross-section area (m2 cm-2)"]<-10000*as.numeric(TRY_Al2As$OrigValueStr[TRY_Al2As$OriglName=="The ratio of leaf area attached per unit sapwood cross-section area (m2 cm-2)"])
TRY_Al2As$StdValue[TRY_Al2As$OriglName=="values at base of living crown, m2/cm2"]<-10000*as.numeric(TRY_Al2As$OrigValueStr[TRY_Al2As$OriglName=="values at base of living crown, m2/cm2"])
TRY_Al2As$StdValue[TRY_Al2As$OriglName=="values at breast height, m2/cm2"]<-10000*as.numeric(TRY_Al2As$OrigValueStr[TRY_Al2As$OriglName=="values at breast height, m2/cm2"])
TRY_Al2As <- TRY_Al2As[TRY_Al2As$ErrorRisk <3, c("AccSpeciesName", "StdValue")]
names(TRY_Al2As) <- c("Species", "Al2As")
hydratry_Al2As <-hydratry[,c("Species","Hv")]
hydratry_Al2As$Hv<- 10000/hydratry_Al2As$Hv
names(hydratry_Al2As) <- c("Species", "Al2As")
Al2As_data<-rbind(TRY_Al2As,
               hydratry_Al2As)
Al2As_fam = TFM_fun(Al2As_data$Species, Al2As_data$Al2As)

# WUE
TRY_WUE  <-readRDS(file.path(DB_path,"TraitDatabases/TRY/TRY_traits/TRY_134.rds"))
TRY_WUE <- TRY_WUE[TRY_WUE$ErrorRisk<3, c("AccSpeciesName", "StdValue")]
WUE_fam = TFM_fun(TRY_WUE$AccSpeciesName, TRY_WUE$StdValue, sp_fun = "max")

#conduit2sapwood ratio
fn_morris <- file.path(DB_path, 
                       "TraitDatabases/Morris_et_al_2016/nph13737-sup-0002-tables1.xlsx")
morris <- openxlsx::read.xlsx(fn_morris, sheet=2)
morris_stem <- morris[morris$Organ=="stem",]
conduit2sapwood_fam = TFM_fun(morris_stem[["Species.name"]], 1 - (morris_stem[["Radial.and.Axial.Parenchyma.(%)"]]/100))

# Compile family mean values
fams = sort(unique(c(names(leaf_dens_fam), 
                     names(wood_dens_fam),
                     names(leaf_PI0_fam),
                     names(leaf_EPS_fam),
                     names(leaf_AF_fam),
                     names(Gswmin_fam),
                     names(Gswmax_fam),
                     names(WoodC_fam),
                     names(Nleaf_fam),
                     names(Nsapwood_fam),
                     names(Nfineroot_fam),
                     names(Ks_fam),
                     names(P50_fam),
                     names(WUE_fam),
                     names(Al2As_fam),
                     names(conduit2sapwood_fam))))
trait_family_means = data.frame(row.names=fams)
trait_family_means$LeafDensity = NA
trait_family_means[names(leaf_dens_fam), "LeafDensity"] = leaf_dens_fam
trait_family_means$WoodDensity = NA
trait_family_means[names(wood_dens_fam), "WoodDensity"] = wood_dens_fam
trait_family_means$LeafPI0 = NA
trait_family_means[names(leaf_PI0_fam), "LeafPI0"] = leaf_PI0_fam
trait_family_means$LeafEPS = NA
trait_family_means[names(leaf_EPS_fam), "LeafEPS"] = leaf_EPS_fam
trait_family_means$LeafAF = NA
trait_family_means[names(leaf_AF_fam), "LeafAF"] = leaf_AF_fam
trait_family_means$Gswmax = NA
trait_family_means[names(Gswmax_fam), "Gswmax"] = Gswmax_fam
trait_family_means$Gswmin = NA
trait_family_means[names(Gswmin_fam), "Gswmin"] = Gswmin_fam
trait_family_means$Nleaf = NA
trait_family_means[names(Nleaf_fam), "Nleaf"] = Nleaf_fam
trait_family_means$Nsapwood = NA
trait_family_means[names(Nsapwood_fam), "Nsapwood"] = Nsapwood_fam
trait_family_means$Nfineroot = NA
trait_family_means[names(Nfineroot_fam), "Nfineroot"] = Nfineroot_fam
trait_family_means$WoodC = NA
trait_family_means[names(WoodC_fam), "WoodC"] = WoodC_fam
trait_family_means$Kmax_stemxylem = NA
trait_family_means[names(Ks_fam), "Kmax_stemxylem"] = Ks_fam
trait_family_means$WUE = NA
trait_family_means[names(WUE_fam), "WUE"] = WUE_fam
trait_family_means$P50 = NA
trait_family_means[names(P50_fam), "P50"] = P50_fam
trait_family_means$Al2As = NA
trait_family_means[names(Al2As_fam), "Al2As"] = Al2As_fam
trait_family_means$conduit2sapwood = NA
trait_family_means[names(conduit2sapwood_fam), "conduit2sapwood"] = conduit2sapwood_fam
write.table(trait_family_means, "data-raw/trait_family_means.csv", sep=";")

