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
\item{\code{transpirationMode (="Simple")}: Transpiration model. See \code{\link{swbInput}}.}
\item{\code{hydraulicCostFunction (= 2)}: Variant of the hydraulic cost function used in the stomatal regulation model of Sperry & Love (2016). Values accepted are 1 (original cost function based on the derivative of supply function), 2 (xylem vulnerability curve), 3 (xylem vulnerability curve, including previous cavitation).}
\item{\code{canopyMode (="sunshade")}: Indicates how crowns should be described to calculate photosynthesis. Accepted values are \code{"sunshade"} (distinguishes photosynthesis in sun leaves from shade leaves) and \code{"multilayer"} (distinguishes photosynthesis of sun leaves and shade leaves in each canopy layer).}
\item{\code{cavitationRefill (= TRUE)}: Whether refilling of embolized conduits is activated.}
\item{\code{taper (=FALSE)}: Whether taper of xylem conduits is accounted for when calculating aboveground stem conductance from xylem conductivity.}
\item{\code{numericParams}: A list with the following elements:
\itemize{
\item{\code{maxNsteps (= 400)}: Maximum number of steps in supply function.}
\item{\code{psiStep (= - 0.001)}: Step in water potential.}
\item{\code{psiMax (= - 10.0)}: Minimum water potential allowed.}
\item{\code{ntrial (= 20)}: Number of iteration trials when finding root of equation system.}
\item{\code{psiTol (= 0.0001)}: Tolerance value for water potential.}
\item{\code{ETol (= 0.0001)}: Tolerance value for flow.}
}
}
\item{\code{fracTotalTreeResistance (=0.625)}: Fraction of the minimum tree (stem+root) resistance that corresponds to stem resistance.}
\item{\code{averageFracRhizosphereResistance (=0.05)}: Fraction to total continuum (stem+root+rhizosphere) resistance that corresponds to rhizosphere (averaged across soil water potential values).}
\item{\code{storagePool (="none")}: Whether carbon storage pools are considered (either "none", "one" or "two").}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
 \code{\link{swbInput}}, \code{\link{swb}}, \code{\link{swbpoints}}
}