\encoding{UTF-8}
\name{SpParams}
\alias{SpParamsMED}
\alias{SpParamsUS}
\docType{data}
\title{
Data tables with species parameter values for Mediterranean and US
}
\description{
A data sets of species parameter values, resulting from bibliographic search, fit to empirical data or expert-based guesses.
}
\usage{
data("SpParamsMED")
data("SpParamsUS")
}
\format{
  A data frame with species as rows and the following columns (variables), grouped by categories.
  \describe{
  \bold{Plant identity}:
  \itemize{
    \item{\code{Name}: A factor with taxon names (mostly species names).}
    \item{\code{IFNcodes}: A factor with levels corresponding to species codes in the Third Spanish forest inventory (DGCN 2005).}
    \item{\code{SpIndex}: A numeric vector of the species index.}
    \item{\code{Group}: Either "Gymnosperm" or "Angiosperm".}
    \item{\code{Order}: Taxonomical order.}
    \item{\code{Family}: Taxonomical family.}
    \item{\code{GrowthForm}: Either \code{Shrub}, \code{Tree} and \code{Tree/Shrub}.}
  }
  \bold{Plant size}:
  \itemize{
    \item{\code{Hmax}: Maximum plant height (in cm).}
    \item{\code{Hmed}: Median plant height (in cm).}
    \item{\code{Z50}: Rooting depth (in mm) corresponding to 50\% of fine roots.}
    \item{\code{Z95}: Rooting depth (in mm) corresponding to 95\% of fine roots.}
  }
  \bold{Shrub allometric coefficients} (calibrated for Catalonia; see De \enc{Cáceres}{Caceres} et al. 2019):
  \itemize{
    \item{\code{a_ash}, \code{b_ash}: Allometric coefficients relating shrub height with the projected area of shrub individuals.}
    \item{\code{a_bsh}, \code{b_bsh}: Allometric coefficients relating phytovolume with dry weight of shrub individuals.}
    \item{\code{cr}: Ratio between shrub crown length and total height.}
  }
  \bold{Tree allometric coefficients} (calibrated for Catalonia):
  \itemize{
    \item{\code{a_fbt}, \code{b_fbt}, \code{c_fbt}, \code{d_fbt}: Regression coefficients used to calculate foliar biomass of an individual tree from its dbh and the cummulative basal area of larger trees.}
    \item{\code{a_cr}, \code{b_1cr}, \code{b_2cr}, \code{b_3cr}, \code{c_1cr}, \code{c_2cr}: Regression coefficients used to calculate crown ratio of trees.}
    \item{\code{a_cw}, \code{b_cw}: Regression coefficients used to calculated crown width of trees.}
  }
  \bold{Phenology}:
  \itemize{  
    \item{\code{PhenologyType}: Leaf phenology type, either "oneflush-evergreen" (new leaves develop in spring-summer), "progressive-evergreen" (new leaves develop during any season), "winter-deciduous" (leaf senescence in autumn, new leaves in spring-summer) or "winter-semideciduous" (same as before, but abscission of senescent leaves occurs when new leaves are produced).}
    \item{\code{LeafDuration}: Leaf duration (in years).}
    \item{\code{Sgdd}: Degree days for leaf budburst (in Celsius).}
    \item{\code{Tbgdd}: Base temperature for the calculation of degree days to leaf budburst.}
    \item{\code{Ssen}: Degree days corresponding to leaf senescence.}
    \item{\code{Phsen}: Photoperiod corresponding to start counting senescence degree-days.}
    \item{\code{Tbsen}: Base temperature for the calculation of degree days to leaf senescence.}
  }
  \bold{Anatomy}:
  \itemize{  
    \item{\code{SLA}: Specific leaf area (mm2/mg = m2/kg).}
    \item{\code{LeafDensity}: Density of leaf tissue (dry weight over volume).}
    \item{\code{WoodDensity}: Wood density (at 0 percent humidity!).}
    \item{\code{FineRootDensity}: Density of fine root tissue (dry weight over volume).}
    \item{\code{conduit2sapwood}: Proportion of sapwood corresponding to conducive elements (vessels or tracheids) as opposed to parenchymatic tissue.}
    \item{\code{r635}: Ratio between the weight of leaves plus branches and the weight of leaves alone for branches of 6.35 mm.}
    \item{\code{pDead}: Proportion of total fine fuels that are dead}
    \item{\code{Al2As}: Leaf area to sapwood area ratio (in m2·m-2).}
    \item{\code{LeafWidth}: Leaf width (in cm).}
    \item{\code{SRL}: Specific root length (cm/g).}
    \item{\code{RLD}: Fine root length density (i.e. density of root length per soil volume) (cm·cm-3).}
  }
  \bold{Tissue moisture}:
  \itemize{  
    \item{\code{maxFMC}: Maximum fuel moisture (in percent of dry weight)}
    \item{\code{minFMC}: Minimum fuel moisture (in percent of dry weight)}
    \item{\code{LeafPI0}: Osmotic potential at full turgor of leaves (MPa).}
    \item{\code{LeafEPS}: Modulus of elasticity (capacity of the cell wall to resist changes in volume in response to changes in turgor) of leaves (MPa).}
    \item{\code{LeafAF}: Apoplastic fraction (proportion of water outside the living cells) in leaves.}
    \item{\code{StemPI0}: Osmotic potential at full turgor of symplastic xylem tissue (MPa).}
    \item{\code{StemEPS}: Modulus of elasticity (capacity of the cell wall to resist changes in volume in response to changes in turgor) of symplastic xylem tissue (Mpa).}
    \item{\code{StemAF}: Apoplastic fraction (proportion of water outside the living cells) in stem xylem.}
    }
  \bold{Flammability}:
  \itemize{  
    \item{\code{LigninPercent}: Percent of lignin+cutin over dry weight in leaves.}
    \item{\code{ParticleDensity}: Particle density (kg/m3).}
    \item{\code{LeafLitterFuelType}: Fuel type for leaf litter, with levels \code{Broadleaved}, \code{LongLinear}, \code{Scale} and \code{ShortLinear}.}
    \item{\code{Flammability}: Flammability modifier (either 1 or 2 for normal or high, respectively).}
    \item{\code{SAV}: Surface-area-to-volume ratio of the small fuel (1h) fraction (leaves and branches < 6.35mm) (m2/m3).}
    \item{\code{HeatContent}: High fuel heat content (kJ/kg).}
    }
  \bold{Extinction}:
  \itemize{  
    \item{\code{gammaSWR}: Reflectance (albedo) for SWR (\code{gammaPAR = 0.8*gammaSWR}).}
    \item{\code{alphaSWR}: Absorbance coefficient for SWR (\code{alphaPAR = alphaSWR * 1.35}).}
    \item{\code{kPAR}: Extinction coeficient for PAR (\code{kSWR = kPAR/1.35}).}
    \item{\code{g}: Canopy water storage capacity per LAI unit (in mm/LAI).}
    }
  \bold{Transpiration/Photosynthesis}:
  \itemize{  
     \item{\code{Tmax_LAI}: Coefficient relating LAI with the ratio of maximum transpiration over potential evapotranspiration.}
    \item{\code{Tmax_LAIsq}: Coefficient relating squared LAI with the ratio of maximum transpiration over potential evapotranspiration.}
    \item{\code{Psi_Extract}: Water potential corresponding to 50\% stomatal conductance (in MPa).}
    \item{\code{Psi_Critic}: Water potential corresponding to 50\% of stem cavitation (in MPa).}
    \item{\code{WUE}: Water use efficiency for carbon assimilation (gross photosynthesis) (g C /mm water).}
    \item{\code{pRootDisc}: Relative root conductance leading to hydraulic disconnection from a soil layer.}
    \item{\code{Gswmin}: Minimum stomatal conductance to water vapor per leaf area unit (in mol·s-1·m-2).}
    \item{\code{Gswmax}: Maximum stomatal conductance to water vapor per leaf area unit (in mol·s-1·m-2).}
    \item{\code{VCleaf_kmax}: Leaf hydraulic conductance (in mmol H2O·s-1·m-2·MPa-1).}
    \item{\code{VCleaf_c}, \code{VCleaf_d}: Parameters of the leaf vulnerability curve (\code{VCleaf_d} in MPa).}
    \item{\code{Kmax_stemxylem}: Sapwood-specific hydraulic conductivity of stem xylem (in kg H2O·s-1·m-1·MPa-1).}
    \item{\code{VCstem_c}, \code{VCstem_d}: Parameters of the stem xylem vulnerability curve (\code{VCstem_d} in MPa).}
    \item{\code{Kmax_rootxylem}: Sapwood-specific hydraulic conductivity of root xylem (in kg H2O·s-1·m-1·MPa-1).}
    \item{\code{VCroot_c}, \code{VCroot_d}: Parameters of the root xylem vulnerability curve (\code{VCroot_d} in MPa).}
    \item{\code{Narea}: Nitrogen contrentration per leaf area (in g·m-2).}
    \item{\code{Vmax298}: Maximum Rubisco carboxilation rate at 25ºC (in micromol CO2·s-1·m-2).}
    \item{\code{Jmax298}: Maximum electron transport rate at 25ºC (in micromol electrons·s-1·m-2).}
  }
  \bold{Growth}:
  \itemize{  
    \item{\code{WoodC}: Wood carbon content per dry weight (g C /g dry).}
    \item{\code{RERleaf}: Maintenance respiration rate for leaves (g gluc · g dry-1 · day-1).}
    \item{\code{RERsapwood}: Maintenance respiration rate for living cells of sapwood (g gluc · g dry-1 · day-1).}
    \item{\code{RERfineroot}: Maintenance respiration rate for fine roots (g gluc · g dry-1 · day-1).}
    \item{\code{RGRleafmax}: Maximum leaf area growth rate relative to sapwood area (in m2·cm-2· day-1).}
    \item{\code{RGRsapwoodmax}: Maximum sapwood area relative growth rate (in cm2·cm-2·day-1).}
    \item{\code{RGRfinerootmax}: Maximum fine root relative growth rate (in g dry·g dry-1· day-1).}
    \item{\code{fHDmin}: Minimum value of the height-to-diameter ratio (dimensionless).}
    \item{\code{fHDmax}: Maximum value of the height-to-diameter ratio (dimensionless).}
    }
  \bold{Recruitment}:
  \itemize{  
    \item{\code{SeedProductionHeight}: Minimum height (cm) for seed production.}
    \item{\code{MinTempRecr}: Minimum average temperature (Celsius) of the coldest month for successful recruitment.}
    \item{\code{MinMoistureRecr}: Minimum value of the moisture index (annual precipitation over annual PET) for successful recruitment.}
    \item{\code{MinFPARRecr}: Minimum percentage of PAR at the ground level for successful recruitment.}
    \item{\code{RecrTreeDBH}: Recruitment DBH (cm) for trees.}
    \item{\code{RecrTreeHeight}: Recruitment height (cm) for trees.}
    \item{\code{RecrTreeDensity}: Recruitment density (ind·ha-1) for trees.}
    \item{\code{RecrShrubCover}: Recruitment cover (\%) for shrubs.}
    \item{\code{RecrShrubHeight}: Recruitment height (cm) for shrubs.}
    \item{\code{RecrZ50}: Soil depth (mm) corresponding to 50\% of fine roots for recruitment.}
    \item{\code{RecrZ95}: Soil depth (mm) corresponding to 95\% of fine roots for recruitment.}
  }
  }
}
\details{
See details of parameterization in De Caceres et al. (2015) and De Caceres et al. (2021).
}
\source{
De \enc{Cáceres}{Caceres} M, \enc{Martínez}{Martinez}-Vilalta J, Coll L, Llorens P, Casals P, Poyatos R, Pausas JG, Brotons L. (2015) Coupling a water balance model with forest inventory data to predict drought stress: the role of forest structural changes vs. climate changes. Agricultural and Forest Meteorology (doi:10.1016/j.agrformet.2015.06.012).

De \enc{Cáceres}{Caceres} M, Mencuccini M, Martin-StPaul N, Limousin JM, Coll L, Poyatos R, Cabon A, Granda V, Forner A, Valladares F, \enc{Martínez}{Martinez}-Vilalta J (2021) Unravelling the effect of species mixing on water use and drought stress in holm oak forests: a modelling approach. Agricultural and Forest Meteorology 296 (doi:10.1016/j.agrformet.2020.108233).

De \enc{Cáceres}{Caceres} M, Casals P, Gabriel E, Castro X (2019) Scaling-up individual-level allometric equations to predict stand-level fuel loading in Mediterranean shrublands. Ann For Sci 76:87 (doi:10.1007/s13595-019-0873-4).


DGCN (2005). Tercer Inventario Forestal Nacional (1997-2007): Catalunya. Dirección General de Conservación de la Naturaleza, Ministerio de Medio Ambiente, Madrid.

}
\seealso{ \code{\link{spwb}}}
\examples{
data(SpParamsMED)
}
\keyword{datasets}
