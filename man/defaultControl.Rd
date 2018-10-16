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
\item{\code{soilFunctions (="SX")}: Soil water retention curve and conductivity functions, either 'SX' (for Saxton) or 'VG' (for Van Genuchten).}
\item{\code{snowpack (=TRUE)}: Boolean flag to indicate the simulation of snow accumulation and melting.}
\item{\code{drainage (=TRUE)}: Boolean flag to indicate the simulation of deep drainage (not used in \code{spwbgrid}).}
\item{\code{transpirationMode (="Simple")}: Transpiration model. See \code{\link{spwbInput}}.}
\item{\code{hydraulicCostFunction (= 1)}: Variant of the hydraulic cost function used in the stomatal regulation model of Sperry & Love (2016). Values accepted are 1 (original cost function based on the derivative of supply function), 2 (leaf vulnerability curve).}
\item{\code{ndailysteps (= 24)}: Number of steps into which each day is divided for determination of stomatal conductance, transpiration and photosynthesis (24 equals 1-hour intervals).}
\item{\code{nsubsteps (= 100)}: Number of temporal substeps for the calculation of capacitance effects (\code{nsubsteps = 100} means substeps of 36 seconds if \code{ndailysteps = 24}).}
\item{\code{canopyMode (="sunshade")}: Indicates how crowns should be described to calculate photosynthesis. Accepted values are \code{"sunshade"} (distinguishes photosynthesis in sun leaves from shade leaves) and \code{"multilayer"} (distinguishes photosynthesis of sun leaves and shade leaves in each canopy layer).}
\item{\code{verticalLayerSize (= 100)}: The size of vertical layers (in cm) for photosynthesis calculation.}
\item{\code{nStemSegments (= 1)}: Number of segments within the stem.}
\item{\code{ksymver (=0)}: Whole-stem symplastic vertical conductance.}
\item{\code{klat (=0)}: Symplastic-apoplastic lateral conductance.}
\item{\code{cavitationRefill (= TRUE)}: Whether refilling of embolized conduits is activated.}
\item{\code{taper (= TRUE)}: Whether taper of xylem conduits is accounted for when calculating aboveground stem conductance from xylem conductivity.}
\item{\code{numericParams}: A list with the following elements:
\itemize{
\item{\code{maxNsteps (= 400)}: Maximum number of steps in supply function.}
\item{\code{psiStep (= - 0.0001)}: Step in water potential (has a large influence in computational time).}
\item{\code{psiMax (= - 10.0)}: Minimum water potential allowed.}
\item{\code{ntrial (= 200)}: Number of iteration trials when finding root of equation system.}
\item{\code{psiTol (= 0.0001)}: Tolerance value for water potential.}
\item{\code{ETol (= 0.0001)}: Tolerance value for flow.}
}
}
\item{\code{thermalCapacityLAI (=1000000)}: Thermal canopy capacitance per LAI unit.}
\item{\code{defaultWindSpeed (=5)}: Default wind speed value (in m/s) to be used when missing from data. }
\item{\code{Catm (=386)}: Atmospheric CO2 concentration (in micromol·mol-1 = ppm).}
\item{\code{averageFracRhizosphereResistance (=0.15)}: Fraction to total continuum (stem+root+rhizosphere) resistance that corresponds to rhizosphere (averaged across soil water potential values).}
\item{\code{storagePool (="none")}: Whether carbon storage pools are considered (either "none", "one" or "two").}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
 \code{\link{spwbInput}}, \code{\link{spwb}}, \code{\link{spwbpoints}}
}