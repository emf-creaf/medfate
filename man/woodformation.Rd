\encoding{UTF-8}
\name{Wood formation}
\alias{woodformation_initRing}
\alias{woodformation_growRing}
\alias{woodformation_relativeExpansionRate}

\title{Wood formation}
\description{
Functions to initialize and expand a ring of tracheids to simulate secondary growth.
}
\usage{
woodformation_initRing()
woodformation_growRing(ring, psi, Tc,
                       Nc=8.85, phi0=0.13, pi0=-0.8, CRD0=8.3,
                       Y_P=0.05, Y_T=8.0, h=0.043*1.8, s=1.8)
woodformation_relativeExpansionRate(psi, Tc, pi, phi, Y_P, Y_T)
}
\arguments{
  \item{ring}{An object of class \code{\link{ring}} returned by function \code{woodformation_initRing}.}
  \item{psi}{Water potential (in MPa).}
  \item{Tc}{Temperature in Celsius.}
  \item{Nc}{Number of active cells in the cambium.}
  \item{phi0}{Initial value of cell extensibility (in MPa-1 day-1)}
  \item{pi0}{Initial value of cell osmotic potential (in MPa)}
  \item{CRD0}{Initial value of cell radial diameter}
  \item{Y_P}{Turgor pressure yield threshold (in MPa)}
  \item{Y_T}{Temperature yield threshold (in Celsius)}
  \item{h}{Cell wall hardening coefficient (in day-1)}
  \item{s}{Cell wall softening coefficient (unitless)}
  \item{pi}{Osmotic potential (in MPa)}
  \item{phi}{Cell extensibility (in MPa-1 day-1)}
}
\value{
Function \code{woodformation_initRing()} returns a list of class 'ring', that is a list containing a data frame \code{cells} and two vectors: \code{P} and \code{SA}. Dataframe \code{cells} contains the columns "formation_date", "phi", "pi" and "CRD" and as many rows as dates processed. Vectors \code{P} and \code{SA} contain, respectively, the number of cells produced and the sapwood area corresponding to the ring of cells (assuming a tangencial radius of 20 micrometers). Function \code{woodformation_growRing()} modifies the input 'ring' object according to the environmental conditions given as input. Function \code{woodformation_relativeExpansionRate()} returns a numeric scalar with the relative expansion rate.
}
\note{
Code modified from package xylomod by Antoine Cabon, available at GitHub
}
\author{
Antoine Cabon, CTFC

Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\references{
Cabon A, \enc{Fernández-de-Uña}{Fernandez-de-Una} L, Gea-Izquierdo G, Meinzer FC, Woodruff DR, \enc{Martínez-Vilalta}{Martinez-Vilalta} J, De \enc{Cáceres}{Caceres} M. 2020a. Water potential control of turgor-driven tracheid enlargement in Scots pine at its xeric distribution edge. New Phytologist 225: 209–221.

Cabon A, Peters RL, Fonti P, \enc{Martínez-Vilalta}{Martinez-Vilalta}  J, De \enc{Cáceres}{Caceres} M. 2020b. Temperature and water potential co-limit stem cambial activity along a steep elevational gradient. New Phytologist: nph.16456.
}
\seealso{
 \code{\link{growth}}
}
