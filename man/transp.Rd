\encoding{UTF-8}
\name{transp}
\alias{transp.profitMaximization}
\alias{transp.dayCanopyTranspiration}
\alias{transp.dayCanopyTranspirationPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transpiration submodel functions
}
\description{
Set of high-level functions used in the calculation of stomatal conductance and transpiration. Function \code{transp.profitMaximization} calculates gain and cost functions, as well as profit maximization from supply and photosynthesis input functions. Function \code{transp.dayCanopyTranspiration} calculates profit maximization for all time steps in a day, starting from stand description, soil description and meteorological input. Function \code{transp.dayCanopyTranspirationPlot} produces a plot with the cohort supply functions against water potential and a plot with the cohort photosynthesis functions against water potential, both with the maximum profit values indicated.
}
\usage{
transp.profitMaximization(supplyFunction, photosynthesisFunction, type, Gwmax, kstemmax = NA)
transp.dayCanopyTranspiration(x, soil, meteo, day,
                             latitude, elevation, slope = 0, aspect = 0)
transp.dayCanopyTranspirationPlot(x, soil, meteo, day, timestep, 
                                  latitude, elevation, slope=0, aspect=0)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{supplyFunction}{Water supply function (see \code{\link{hydraulics.supplyFunctionNetwork}}).}
  \item{photosynthesisFunction}{Function returned by \code{photo.photosynthesisFunction()}.}
  \item{type}{Optimization type (1,2 or 3).}
  \item{Gwmax}{Maximum stomatal conductance (mol·m-2·s-1).}
  \item{kstemmax}{Maximum stem conductance (flow per pressure drop).}
  \item{x}{An object of class \code{\link{swbInput}} built using the 'Complex' transpiration mode.}
  \item{soil}{An object of class \code{\link{soil}}.}
  \item{meteo}{A data frame with daily meteorological data series:
   \itemize{
     \item{\code{DOY}: Day of the year (Julian day).}
     \item{\code{Precipitation}: Precipitation (in mm).}
     \item{\code{MeanTemperature}: Mean temperature (in degrees Celsius).}
     \item{\code{MinTemperature}: Minimum temperature (in degrees Celsius).}
     \item{\code{MaxTemperature}: Maximum temperature (in degrees Celsius).}
     \item{\code{MinRelativeHumidity}: Minimum relative humidity (in percent).}
     \item{\code{MaxRelativeHumidity}: Maximum relative humidity (in percent).}
     \item{\code{Radiation}: Solar radiation (in MJ/m2/day).}
     \item{\code{WindSpeed}: Wind speed (in m/s). If not available, this column can be left with \code{NA} values.}
    }
  }  
  \item{day}{An integer to identify a day within \code{meteo}.}
  \item{timestep}{An integer between 1 and \code{ndailysteps} specified in \code{x} (see \code{\link{defaultControl}}).}
  \item{latitude}{Latitude (in degrees).}
  \item{elevation, slope, aspect}{Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North).}
}
\details{
Details of the transpiration submodel are given in a vignette.
}

\value{
Values returned for each function are:
\itemize{
\item{\code{transp.profitMaximization}: A list with the following elements:
\itemize{
  \item{\code{Cost}: Cost function [0-1].}
  \item{\code{Gain}: Gain function [0-1].}
  \item{\code{Profit}: Profit function [0-1].}
  \item{\code{iMaxProfit}: Index corresponding to maximum profit (starting from 0).}
}
}
\item{\code{transp.dayCanopyTranspiration}: A list whose elements correspond to cohorts. Each element has the following subelements:
\itemize{
  \item{\code{supply}: Supply function (see \code{\link{hydraulics.supplyFunctionNetwork}}).}
  \item{\code{photo}: Cohort photosynthesis function for each time step (see \code{\link{photo.sunshadePhotosynthesisFunction}} or \code{\link{photo.multilayerPhotosynthesisFunction}}).}
  \item{\code{PM}: Profit maximization for each time step (see \code{transp.profitMaximization}).}
}
}
}
}
\references{
Sperry, J. S., M. D. Venturas, W. R. L. Anderegg, M. Mencuccini, D. S. Mackay, Y. Wang, and D. M. Love. 2016. Predicting stomatal responses to the environment from the optimization of photosynthetic gain and hydraulic cost. Plant Cell and Environment.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
\code{\link{hydraulics.supplyFunctionNetwork}}, \code{\link{biophysics.leafTemperature}}, \code{\link{photo.photosynthesis}}, \code{\link{swb}}
}
\examples{
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Initialize soil with default soil params (2 layers)
examplesoil2 = soil(defaultSoilParams(2))

#Initialize control parameters
control = defaultControl()
#Switch to 'Complex' transpiration mode
control$transpirationMode="Complex"

#Initialize input
x2 = forest2swbInput(exampleforest,examplesoil2, SpParamsMED, control)

transp.dayCanopyTranspirationPlot(x2, examplesoil2, examplemeteo, day=100, timestep = 12,
                                  latitude = 41.82592, elevation = 100)
}