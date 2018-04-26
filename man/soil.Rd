\encoding{UTF-8}
\name{soil}
\alias{soil}
\alias{print.soil}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Soil initialization}
\description{
Initializes a soil object for its use in simulation functions.
}
\usage{
soil(SoilParams, W = as.numeric(c(1)))
\method{print}{soil}(x, model="SX", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SoilParams}{A list of soil parameters (see \code{\link{defaultSoilParams}}).}
  \item{W}{A numerical vector with the initial relative water content of each soil layer.}
  \item{x}{An object of class \code{soil}.}
  \item{model}{Either 'SX' or 'VG' for Saxton or Van Genuchten pedotransfer models.}
  \item{...}{Additional parameters to \code{print}.}
}
\details{
Function \code{print} prompts a description of soil characteristics and status (water content) according to a pedotransfer model (either Saxton's or Van Genuchten's). Volume at field capacity is calculated assuming a soil water potential equal to -0.033 MPa. Parameter \code{Temp} is initialized as missing for all layers.
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
\item{\code{macro}: Macroporosity for each layer.}
\item{\code{rfc}: Percentage of rock fragment content for each layer.}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
  \code{\link{soil.psi2thetaSX}}, \code{\link{soil.psi2thetaVG}}, \code{\link{swb}}, \code{\link{defaultSoilParams}}
}
\examples{
# Initializes soil
s = soil(defaultSoilParams())

# Prints soil characteristics according to Saxton's pedotransfer functions
print(s, model="SX")

# Prints soil characteristics according to Van Genuchten's pedotransfer functions
print(s, model="VG")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.