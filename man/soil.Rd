\encoding{UTF-8}
\name{soil}
\alias{soil}
\alias{print.soil}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Soil definition}
\description{
Defines a soil object for its use in simulation functions.
}
\usage{
soil(SoilParams, W = as.numeric(c(1)))
\method{print}{soil}(x,  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SoilParams}{A list of soil parameters (see \code{\link{defaultSoilParams}}).}
  \item{W}{A numerical vector with the initial relative water content of each soil layer.}
  \item{x}{An object of class \code{soil}.}
  \item{...}{Additional parameters to \code{print}.}
}
\details{
Function \code{print} prompts a description of soil characteristics and status (water content). Volume at field capacity is calculated from assuming a soil water potential equal to -0.033 MPa.
}
\value{
An list of class \code{soil} with the following elements:
\itemize{
\item{\code{SoilDepth}: Soil depth (in mm).}
\item{\code{W}: State variable with relative water content of each layer (in as percent relative to FC).}
\item{\code{psi}: State variable with water potential (in MPa) of each layer.}
\item{\code{Temp}: State variable with temperature (in ºC) of each layer.}
\item{\code{Ksoil}: Kappa parameter for infiltration.}
\item{\code{Gsoil}: Gamma parameter for infiltration.}
\item{\code{dVec}: Width of soil layers (in mm).}
\item{\code{sand}: Sand percentage for each layer (in percent volume).}
\item{\code{clay}: Clay percentage for each layer (in percent volume).}
\item{\code{om}: Organic matter percentage for each layer (in percent volume).}
\item{\code{usda_Type}: USDA texture type.}
\item{\code{VG_alpha}, \code{VG_n}: Parameters for van Genuchten's conductance function, for each layer, corresponding to the USDA texture type.}
\item{\code{macro}: Macroporosity for each layer.}
\item{\code{rfc}: Percentage of rock fragment content for each layer.}
\item{\code{Theta_FC}: Soil moisture at field capacity (in percent volume).}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
  \code{\link{soil.psi2theta}}, \code{\link{swb}}, \code{\link{defaultSoilParams}}
}
\examples{
s = soil(defaultSoilParams())
s
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.