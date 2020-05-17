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
\title{
Carbon-related functions
}
\description{
Set of functions used in the calculation of carbon balance.
}
\usage{
carbon_leafStarchCapacity(LAI, N, SLA, leafDensity)
carbon_leafStructuralBiomass(LAI, N, SLA)
carbon_sapwoodStarchCapacity(SA, H, Z, woodDensity, vessel2sapwood)
carbon_sapwoodStructuralBiomass(SA, H, Z, woodDensity)
carbon_sapwoodStructuralLivingBiomass(SA, H, Z, woodDensity, vessel2sapwood)
carbon_sugarConcentration(osmoticWP, temp, nonSugarConc)
carbon_osmoticWaterPotential(sugarConc, temp, nonSugarConc)
carbon_relativeSapViscosity(sugarConc, temp)
carbon_sugarStarchDynamicsLeaf(sugarConc, starchConc, eqSugarConc)
}
\arguments{
  \item{LAI}{Leaf area index.}
  \item{N}{Density (ind·ha-1).}
  \item{SLA}{Specific leaf area (mm2/mg = m2/kg).}
  \item{leafDensity}{ Density of leaf tissue (dry weight over volume).}
  \item{SA}{Sapwood area (cm2).}
  \item{H}{Plant height (cm).}
  \item{Z}{Rooting depth (mm).}
  \item{woodDensity}{Wood density (dry weight over volume).}
  \item{vessel2sapwood}{Proportion of sapwood corresponding to vessels.}
  \item{osmoticWP}{Osmotic water potential (MPa).}
  \item{temp}{Temperature (degrees Celsius).}
  \item{nonSugarConc}{Concentration of inorganic solutes (mol/l).}
  \item{sugarConc}{Concentration of soluble sugars (mol/l).}
  \item{starchConc}{Concentration of starch (mol/l)}
  \item{eqSugarConc}{Equilibrium concentration of soluble sugars (mol/l).}
}
\value{
Values returned for each function are:
\itemize{
\item{\code{carbon_leafStarchCapacity}: Capacity of storing starch in the leaf compartment (mol glucose).}
\item{\code{carbon_sapwoodStarchCapacity}: Capacity of storing starch in the sapwood compartment (mol glucose).}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{
\code{\link{growth}}
}