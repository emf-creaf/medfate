\encoding{UTF-8}
\name{soil}
\alias{soil}
\alias{print.soil}
\alias{redefineSoilLayers}

\title{Soil initialization}
\description{
Initializes soil parameters and state variables for its use in simulations.
}
\usage{
soil(SoilParams, VG_PTF = "Toth", W = as.numeric(c(1)), SWE = 0)
\method{print}{soil}(x, model="SX", ...)
redefineSoilLayers(SoilParams, widths = c(300, 700, 1000, 2000))
}
\arguments{
  \item{SoilParams}{A data frame of soil parameters (see an example in \code{\link{defaultSoilParams}}).}
  \item{VG_PTF}{Pedotransfer functions to obtain parameters for the van Genuchten-Mualem equations. Either \code{"Carsel"} (Carsel and Parrish 1988) or \code{"Toth"} (Toth et al. 2015).}
  \item{W}{A numerical vector with the initial relative water content of each soil layer.}
  \item{SWE}{Initial snow water equivalent of the snow pack on the soil surface (mm).}
  \item{x}{An object of class \code{soil}.}
  \item{model}{Either 'SX' or 'VG' for Saxton or Van Genuchten pedotransfer models.}
  \item{widths}{A numeric vector indicating the desired layer widths, in \emph{mm}.}
  \item{...}{Additional parameters to \code{print}.}
}
\details{
Function \code{print} prompts a description of soil characteristics and state variables (water content and temperature) according to a water retention curve (either Saxton's or Van Genuchten's). Volume at field capacity is calculated assuming a soil water potential equal to -0.033 MPa. Parameter \code{Temp} is initialized as missing for all soil layers. Function \code{redefineSoilLayers} allows redefining soil layer widths of an input data frame of soil parameters.
}
\value{
Function \code{soil} returns a list of class \code{soil} with the following elements:
\itemize{
\item{\code{SoilDepth}: Soil depth (in mm).}
\item{\code{W}: State variable with relative water content of each layer (in as proportion relative to FC).}
\item{\code{Temp}: State variable with temperature (in ºC) of each layer.}
\item{\code{Ksoil}: Kappa parameter for infiltration.}
\item{\code{Gsoil}: Gamma parameter for infiltration.}
\item{\code{dVec}: Width of soil layers (in mm).}
\item{\code{sand}: Sand percentage for each layer (in percent volume).}
\item{\code{clay}: Clay percentage for each layer (in percent volume).}
\item{\code{om}: Organic matter percentage for each layer (in percent volume).}
\item{\code{VG_alpha}, \code{VG_n}, \code{VG_theta_res}, \code{VG_theta_sat}: Parameters for van Genuchten's pedotransfer functions, for each layer, corresponding to the USDA texture type.}
\item{\code{Ksat}: Saturated soil conductivity for each layer (estimated using function \code{\link{soil_saturatedConductivitySX}}.}
\item{\code{macro}: Macroporosity for each layer (estimated using Stolf et al. 2011).}
\item{\code{rfc}: Percentage of rock fragment content for each layer.}
\item{\code{Kdrain}: Saturated vertical hydraulic conductivity (mm/day) (i.e. how easy is deep drainage towards groundwater). Function \code{soil} estimates it as a function of soil saturated hydraulic conductivity, but should be parametrized as a function of bedrock material. }
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\references{
Carsel, R.F., and Parrish, R.S. 1988. Developing joint probability distributions of soil water retention characteristics. Water Resources Research 24: 755–769.

\enc{Tóth}{Toth}, B., Weynants, M., Nemes, A., \enc{Makó}{Mako}, A., Bilas, G., and \enc{Tóth}{Toth}, G. 2015. New generation of hydraulic pedotransfer functions for Europe. European Journal of Soil Science 66: 226–238.

Stolf, R., Thurler, A., Oliveira, O., Bacchi, S., Reichardt, K., 2011. Method to estimate soil macroporosity and microporosity based on sand content and bulk density. Rev. Bras. Ciencias do Solo 35, 447–459.
}
\seealso{
  \code{\link{soil_psi2thetaSX}}, \code{\link{soil_psi2thetaVG}}, \code{\link{spwb}}, \code{\link{defaultSoilParams}}
}
\examples{
# Initializes soil
s = soil(defaultSoilParams())

# Prints soil characteristics according to Saxton's water retention curve
print(s, model="SX")

# Prints soil characteristics according to Van Genuchten's water retention curve
print(s, model="VG")

}
