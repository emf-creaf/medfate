\encoding{UTF-8}
\name{carbon}
\alias{carbon_leafStarchCapacity}
\alias{carbon_leafStructuralBiomass}
\alias{carbon_sapwoodStarchCapacity}
\alias{carbon_sapwoodStructuralBiomass}
\alias{carbon_sapwoodStructuralLivingBiomass}
\alias{carbon_relativeSapViscosity}
\alias{carbon_sugarConcentration}
\alias{carbon_osmoticWaterPotential}
\alias{carbon_sugarStarchDynamicsLeaf}
\alias{carbon_sugarStarchDynamicsStem}
\alias{carbon_carbonCompartments}
\title{
Carbon-related functions
}
\description{
Set of functions used in the calculation of carbon balance.
}
\usage{
carbon_leafStarchCapacity(LAI, N, SLA, leafDensity)
carbon_leafStructuralBiomass(LAI, N, SLA)
carbon_sapwoodStarchCapacity(SA, H, L, V, woodDensity, conduit2sapwood)
carbon_sapwoodStructuralBiomass(SA, H, L, V, woodDensity)
carbon_sapwoodStructuralLivingBiomass(SA, H, L, V, woodDensity, conduit2sapwood)
carbon_sugarConcentration(osmoticWP, temp, nonSugarConc)
carbon_osmoticWaterPotential(sugarConc, temp, nonSugarConc)
carbon_relativeSapViscosity(sugarConc, temp)
carbon_sugarStarchDynamicsLeaf(sugarConc, starchConc, eqSugarConc)
carbon_sugarStarchDynamicsStem(sugarConc, starchConc, eqSugarConc)
carbon_carbonCompartments(x, biomassUnits = "g_m2")
}
\arguments{
  \item{LAI}{Leaf area index.}
  \item{N}{Density (ind·ha-1).}
  \item{SLA}{Specific leaf area (mm2/mg = m2/kg).}
  \item{leafDensity}{ Density of leaf tissue (dry weight over volume).}
  \item{SA}{Sapwood area (cm2).}
  \item{H}{Plant height (cm).}
  \item{L}{Coarse root length (mm) for each soil layer.}
  \item{V}{Proportion of fine roots in each soil layer.}
  \item{woodDensity}{Wood density (dry weight over volume).}
  \item{conduit2sapwood}{Proportion of sapwood corresponding to conducive elements (vessels or tracheids) as opposed to parenchymatic tissue.}
  \item{osmoticWP}{Osmotic water potential (MPa).}
  \item{temp}{Temperature (degrees Celsius).}
  \item{nonSugarConc}{Concentration of inorganic solutes (mol/l).}
  \item{sugarConc}{Concentration of soluble sugars (mol/l).}
  \item{starchConc}{Concentration of starch (mol/l)}
  \item{eqSugarConc}{Equilibrium concentration of soluble sugars (mol/l).}
  \item{x}{An object of class \code{\link{growthInput}}.}
  \item{biomassUnits}{A string for output biomass units, either "g_ind" (g per individual) or "g_m2" (g per square meter).}
}
\value{
Values returned for each function are:
\itemize{
\item{\code{carbon_leafStarchCapacity}: Capacity of storing starch in the leaf compartment (mol gluc/ind.).}
\item{\code{carbon_leafStructuralBiomass}: Leaf structural biomass (g dry/ind.)}
\item{\code{carbon_sapwoodStarchCapacity}: Capacity of storing starch in the sapwood compartment (mol gluc/ind.).}
\item{\code{carbon_sapwoodStructuralBiomass}: Sapwood structural biomass (g dry/ind.)}
\item{\code{carbon_sapwoodStructuralLivingBiomass}: Living sapwood (parenchyma) structural biomass (g dry/ind.)}
\item{\code{carbon_sugarConcentration}: Sugar concentration (mol gluc/l)}
\item{\code{carbon_osmoticWaterPotential}: Osmotic component of water potential (MPa)}
\item{\code{carbon_relativeSapViscosity}: Relative viscosity of sapwood with respect to pure water.}
\item{\code{carbon_sugarStarchDynamicsLeaf}: Rate of conversion from sugar to starch in leaf (mol gluc/l/s).}
\item{\code{carbon_sugarStarchDynamicsStem}: Rate of conversion from sugar to starch in leaf (mol gluc/l/s).}
\item{\code{carbon_carbonCompartments}: A data frame with the size of compartments for each plant cohort, in the specified units.}

}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{
\code{\link{growth}}
}
