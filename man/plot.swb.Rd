\encoding{UTF-8}
\name{plot.swb}
\alias{plot.swb}
\alias{plot.growth}
\alias{summary.swb}
\alias{summary.growth}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Displays simulation results}
\description{
Function \code{plot} plots the results of the soil water balance model (see \code{\link{swb}}) or the forest growth model (see \code{\link{growth}}), whereas function \code{summary} summarizes the model's output in different temporal steps (i.e. weekly, annual, ...).
}
\usage{
\method{plot}{swb}(x, type="PET_Precipitation", yearAxis=FALSE, 
    xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL,...)
\method{plot}{growth}(x, type="PET_Precipitation", yearAxis=FALSE, 
    xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL,...)
\method{summary}{swb}(object, freq="years", output="DailyBalance", FUN=sum, ...)
\method{summary}{growth}(object, freq="years", output="DailyBalance", FUN=sum, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, object}{An object of class \code{swb}}
  \item{type}{The information to be plotted:
  \itemize{
    \item{\code{"PET_NetPrec"}:}{ Potential evapotranspiration and Net Precipitation.}
    \item{\code{"Export"}:}{ Water exported through deep drainage and surface runoff.}
    \item{\code{"ET"}:}{ Plant transpiration and soil evaporation.}
    \item{\code{"Psi"}:}{ Soil water potential.}
    \item{\code{"Theta"}:}{ Soil relative water content.}
    \item{\code{"Vol"}:}{ Soil water volumetric content.}
    \item{\code{"PlantStress"}:}{ Plant cohort average daily drought stress.}
    \item{\code{"PlantPsi"}:}{ Plant cohort water potential.}
    \item{\code{"AirTemperature"}:}{ Minimum/maximum/mean daily temperatures above canopy.}
    \item{\code{"CanopyTemperature"}:}{ Minimum/maximum/mean daily temperatures inside canopy.}
    \item{\code{"SoilTemperature"}:}{ Minimum/maximum/mean daily temperatures inside the first soil layer.}
   }
  }
  \item{yearAxis}{A boolean to indicate whether the units of the x-axis are years (by default they are dates).}
  \item{xlim}{Range of values for x.}
  \item{ylim}{Range of values for y.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{freq}{Frequency of summary statistics (see \code{\link{cut.Date}}).}
  \item{output}{The data table to be summarized. Accepted values are \code{"DailyBalance"}, \code{"PlantStress"}, \code{"PlantTranspiration"} and \code{"SoilWaterBalance"}.}
  \item{FUN}{The function to summarize results (e.g., \code{sum}, \code{mean}, ...)}
  \item{...}{Additional parameters for functions \code{plot} or \code{summary}.}
}
\references{
De \enc{Cáceres}{Caceres} M, \enc{Martínez}{Martinez}-Vilalta J, Coll L, Llorens P, Casals P, Poyatos R, Pausas JG, Brotons L. (submitted) Coupling a water balance model with forest inventory data to evaluate plant drought stress at the regional level. Agricultural and Forest Meteorology.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
 \code{\link{swb}}, \code{\link{swbpoints}}
}
\examples{
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Initialize soil with default soil params (2 layers)
examplesoil = soil(defaultSoilParams(2))

#Initialize control parameters
control = defaultControl()

#Initialize input
x = forest2swbInput(exampleforest,examplesoil, SpParamsMED, control)

#Call simulation function
S1<-swb(x, examplesoil, examplemeteo)
        
#Plot results
plot(S1)

#Monthly summary (averages) of soil water balance
summary(S1, freq="months",FUN=mean, output="SoilWaterBalance")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.