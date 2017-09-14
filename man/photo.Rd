\encoding{UTF-8}
\name{photo}
\alias{photo.GammaTemp}
\alias{photo.KmTemp}
\alias{photo.VmaxTemp}
\alias{photo.JmaxTemp}
\alias{photo.electronLimitedPhotosynthesis}
\alias{photo.rubiscoLimitedPhotosynthesis}
\alias{photo.photosynthesis}
\alias{photo.photosynthesisFunction}
\alias{photo.profitMaximization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Photosynthesis submodel functions
}
\description{
Set of functions used in the calculation of photosynthesis.
}
\usage{
photo.GammaTemp(leaf_temp)
photo.KmTemp(leaf_temp, Oi = 209)
photo.VmaxTemp(Vmax298, leaf_temp)
photo.JmaxTemp(Jmax298, leaf_temp)
photo.electronLimitedPhotosynthesis(Q, Ci, GT, Jmax)
photo.rubiscoLimitedPhotosynthesis(Ci, GT, Km, Vmax)
photo.photosynthesis(Q, Catm, Gc, leaf_temp, Vmax298, Jmax298, verbose = FALSE)
photo.leafPhotosynthesisFunction(supplyFunction, Catm, Patm, Tair, vpa, u, 
                             absRad, Q, Vmax298, Jmax298, verbose = FALSE)
photo.multilayerPhotosynthesisFunction(supplyFunction, Catm, Patm, Tair, vpa, 
                                      SLarea, SHarea, u, 
                                      absRadSL, absRadSH, QSL, QSH, 
                                      Vmax298, Jmax298, Gwmin, Gwmax, verbose = FALSE)
photo.profitMaximization(supplyFunction, photosynthesisFunction, type = 1, kstemmax = NA)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{leaf_temp}{Leaf temperature (in ºC).}
  \item{Oi}{Oxigen concentration (mmol*mol-1).}
  \item{Vmax298}{Maximum Rubisco carboxylation rate per leaf area at 298ºK (i.e. 25 ºC) (micromol*s-1*m-2).}
  \item{Jmax298}{Maximum electron transport rate per leaf area at 298ºK (i.e. 25 ºC) (micromol*s-1*m-2).}
  \item{Q}{Active photon flux density (micromol * s-1 * m-2).}
  \item{Ci}{CO2 internal concentration (micromol * mol-1).}
  \item{GT}{CO2 saturation point corrected by temperature (micromol * mol-1).}
  \item{Jmax}{Maximum electron transport rate per leaf area (micromol*s-1*m-2).}
  \item{Km}{Km = Kc*(1.0+(Oi/Ko)) - Michaelis-Menten term corrected by temperature (in micromol * mol-1).}
  \item{Vmax}{Maximum Rubisco carboxylation rate per leaf area (micromol*s-1*m-2).}
  \item{Catm}{CO2 air concentration (micromol * mol-1).}
  \item{Gc}{CO2 leaf (stomatal) conductance (mol * s-1 * m-2).}
  \item{supplyFunction}{Water supply function (see \code{\link{hydraulics.supplyFunctionNetwork}}).}
  \item{Patm}{Atmospheric air pressure (in kPa).}
  \item{Tair}{Air temperature (in ºC).}
  \item{vpa}{Vapour pressure deficit (in kPa).}
  \item{u}{Wind speed above the leaf boundary layer (in m/s).}
  \item{absRad}{Absorbed long- and short-wave radiation (in W*m^-2).}
  \item{verbose}{Boolean flag to indicate console output.}
  \item{photosynthesisFunction}{Function returned by \code{photo.photosynthesisFunction()}.}
  \item{type}{Optimization type (1,2 or 3).}
  \item{kstemmax}{Maximum stem conductance (flow per pressure drop).}
%  \item{Gwmin, Gwmax}{Minimum and maximum leaf water conductance (i.e. cuticular and maximum stomatal conductance) (mol * s-1 * m-2).}
}
\details{
Details of the photosynthesis submodel are given in a vignette.
}
\value{
Values returned for each function are:
\itemize{
\item{\code{photo.GammaTemp}: CO2 compensation concentration (micromol * mol-1).}
\item{\code{photo.KmTemp}: Michaelis-Menten coefficients of Rubisco for Carbon (micromol * mol-1) and Oxigen (mmol * mol-1).}
\item{\code{photo.VmaxTemp}: Temperature correction of Vmax298.}
\item{\code{photo.JmaxTemp}: Temperature correction of Jmax298.}
\item{\code{photo.electronLimitedPhotosynthesis}: Electron-limited photosynthesis (micromol*s-1*m-2) following Farquhar et al. (1980).}
\item{\code{photo.rubiscoLimitedPhotosynthesis}: Rubisco-limited photosynthesis (micromol*s-1*m-2) following Farquhar et al. (1980).}
\item{\code{photo.photosynthesis}: Calculates gross photosynthesis (micromol*s-1*m-2) following (Farquhar et al. (1980) and Collatz et al (1991).}
\item{\code{photo.photosynthesisFunction}: Returns a list with the following elements:
\itemize{
  \item{\code{LeafTemperature}: Leaf temperature (ºC).}
  \item{\code{LeafVPD}: Leaf vapor pressure deficit (kPa).}
  \item{\code{WaterVaporConductance}: Leaf vapor conductance (mol * s-1 * m-2).}
  \item{\code{Photosynthesis}: Gross photosynthesis (micromol*s-1*m-2).}
  }
}
\item{\code{photo.profitMaximization}:Returns a list with the following elements:
\itemize{
  \item{\code{Cost}: Cost function [0-1].}
  \item{\code{Gain}: Gain function [0-1].}
  \item{\code{Profit}: Profit function [0-1].}
  \item{\code{iMaxProfit}: Index corresponding to maximum profit (starting from 0).}
}
}
}
}
\references{
Bernacchi, C. J., E. L. Singsaas, C. Pimentel, A. R. Portis, and S. P. Long. 2001. Improved temperature response functions for models of Rubisco-limited photosynthesis. Plant, Cell and Environment 24:253–259.

Collatz, G. J., J. T. Ball, C. Grivet, and J. A. Berry. 1991. Physiological and environmental regulation of stomatal conductance, photosynthesis and transpiration: a model that includes a laminar boundary layer. Agricultural and Forest Meteorology 54:107–136.

Farquhar, G. D., S. von Caemmerer, and J. A. Berry. 1980. A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species. Planta 149:78–90.

Leuning, R. 2002. Temperature dependence of two parameters in a photosynthesis model. Plant, Cell and Environment 25:1205–1210.

Sperry, J. S., M. D. Venturas, W. R. L. Anderegg, M. Mencuccini, D. S. Mackay, Y. Wang, and D. M. Love. 2016. Predicting stomatal responses to the environment from the optimization of photosynthetic gain and hydraulic cost. Plant Cell and Environment.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
\code{\link{hydraulics.supplyFunctionNetwork}}, \code{\link{biophysics.leafTemperature}}, \code{\link{swb}}
}