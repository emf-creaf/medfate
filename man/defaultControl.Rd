\encoding{UTF-8}
\name{defaultControl}
\alias{defaultControl}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Default control parameters for models}
\description{
Creates a list with global default parameters for simulation models.
}
\usage{
defaultControl() 
}
\details{
The function returns a list with default parameters. Users can change those defaults that need to be set to other values and use the list as input for model functions. The relevant parameters are different for each model function.
}
\value{
A list, with the following options:
\itemize{
\item{\code{verbose (=TRUE)}: Boolean flag to indicate console output during calculations.}
\item{\code{subdailyResults (=FALSE)}: Boolean flag to force subdaily results to be stored (as a list called 'subdaily' of \code{\link{spwb_day}} objects, one by simulated date) in calls to \code{\link{spwb}}.}
\item{\code{soilFunctions (="SX")}: Soil water retention curve and conductivity functions, either 'SX' (for Saxton) or 'VG' (for Van Genuchten).}
\item{\code{snowpack (=TRUE)}: Boolean flag to indicate the simulation of snow accumulation and melting.}
\item{\code{drainage (=TRUE)}: Boolean flag to indicate the simulation of deep drainage.}
\item{\code{leafPhenology (=TRUE)}: Boolean flag to indicate the simulation of leaf phenology for winter-deciduous species.}
\item{\code{transpirationMode (="Granier")}: Transpiration model (either 'Granier' or 'Sperry'). See \code{\link{spwbInput}}.}
\item{\code{hydraulicCostFunction (= 1)}: Variant of the hydraulic cost function used in the stomatal regulation model of Sperry & Love (2016). Values accepted are 1 (original cost function based on the derivative of supply function), 2 (leaf vulnerability curve).}
\item{\code{ndailysteps (= 24)}: Number of steps into which each day is divided for determination of stomatal conductance, transpiration and photosynthesis (24 equals 1-hour intervals).}
\item{\code{verticalLayerSize (= 100)}: The size of vertical layers (in cm) for photosynthesis calculation.}
\item{\code{nStemSegments (= 1)}: Number of segments within the stem.}
\item{\code{capacitance (=FALSE)}: Whether capacitance is considered in simulations.}
\item{\code{cavitationRefill (= TRUE)}: Whether refilling of embolized conduits is activated.}
\item{\code{klat (= 0.1)}: Symplastic-apoplastic lateral conductance.}
\item{\code{taper (= TRUE)}: Whether taper of xylem conduits is accounted for when calculating aboveground stem conductance from xylem conductivity.}
\item{\code{numericParams}: A list with the following elements:
\itemize{
\item{\code{maxNsteps (= 400)}: Maximum number of steps in supply function.}
\item{\code{ntrial (= 200)}: Number of iteration trials when finding root of equation system.}
\item{\code{psiTol (= 0.0001)}: Tolerance value for water potential.}
\item{\code{ETol (= 0.0001)}: Tolerance value for flow.}
}
}
\item{\code{thermalCapacityLAI (=1000000)}: Thermal canopy capacitance per LAI unit.}
\item{\code{defaultWindSpeed (= 2.5)}: Default wind speed value (in m/s) to be used when missing from data. }
\item{\code{Catm (=386)}: Atmospheric CO2 concentration (in micromol·mol-1 = ppm).}
\item{\code{fracRootResistance (=0.40)}: Fraction of plant total resistance (leaf+stem+root) that corresponds to root system. If \code{fracRootResistance = NA} then root conductance is calculated from root xylem conductivity and root distribution.}
\item{\code{averageFracRhizosphereResistance (=0.15)}: Fraction to total continuum (leaf+stem+root+rhizosphere) resistance that corresponds to rhizosphere (averaged across soil water potential values).}
\item{\code{storagePool (="none")}: Whether carbon storage pools are considered (either "none", "one" or "two").}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{
 \code{\link{spwbInput}}, \code{\link{spwb}}
}