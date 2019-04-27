\encoding{UTF-8}
\name{soil}
\alias{soil}
\alias{print.soil}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Soil initialization}
\description{
Initializes soil parameters and state variables for its use in simulations.
}
\usage{
soil(SoilParams, VG_PTF = "Toth", W = as.numeric(c(1)), SWE = 0)
\method{print}{soil}(x, model="SX", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SoilParams}{A list of soil parameters (see \code{\link{defaultSoilParams}}).}
  \item{VG_PTF}{Pedotransfer functions to obtain parameters for the van Genuchten-Mualem equations. Either \code{"Carsel"} (Carsel & Parrish 1988) or \code{"Toth"} (Toth et al. 2015).}
  \item{W}{A numerical vector with the initial relative water content of each soil layer.}
  \item{SWE}{Initial snow water equivalent of the snow pack on the soil surface (mm).}
  \item{x}{An object of class \code{soil}.}
  \item{model}{Either 'SX' or 'VG' for Saxton or Van Genuchten pedotransfer models.}
  \item{...}{Additional parameters to \code{print}.}
}
\details{
Function \code{print} prompts a description of soil characteristics and state variables (water content and temperature) according to a water retention curve (either Saxton's or Van Genuchten's). Volume at field capacity is calculated assuming a soil water potential equal to -0.033 MPa. Parameter \code{Temp} is initialized as missing for all soil layers.
}
\value{
An list of class \code{soil} with the following elements:
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
\item{\code{usda_Type}: USDA texture type.}
\item{\code{VG_alpha}, \code{VG_n}, \code{VG_theta_res}, \code{VG_theta_sat}: Parameters for van Genuchten's pedotransfer functions, for each layer, corresponding to the USDA texture type.}
\item{\code{macro}: Macroporosity for each layer (estimated using Stolf et al. 2011).}
\item{\code{rfc}: Percentage of rock fragment content for each layer.}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\references{
Carsel, R.F., & Parrish, R.S. 1988. Developing joint probability distributions of soil water retention characteristics. Water Resources Research 24: 755–769.

\enc{Tóth}{Toth}, B., Weynants, M., Nemes, A., \enc{Makó}{Mako}, A., Bilas, G., & \enc{Tóth}{Toth}, G. 2015. New generation of hydraulic pedotransfer functions for Europe. European Journal of Soil Science 66: 226–238.

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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.