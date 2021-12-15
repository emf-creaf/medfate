\encoding{UTF-8}
\name{summary.spwb}
\alias{summary.pwb}
\alias{summary.spwb}
\alias{summary.growth}

\title{Summarises simulation results}
\description{
Function \code{summary} summarizes the model's output in different temporal steps (i.e. weekly, annual, ...).
}
\usage{
\method{summary}{pwb}(object, freq="years", output="WaterBalance", FUN=sum, bySpecies = FALSE, ...)
\method{summary}{spwb}(object, freq="years", output="WaterBalance", FUN=sum, bySpecies = FALSE, ...)
\method{summary}{growth}(object, freq="years", output="WaterBalance", FUN=sum, bySpecies = FALSE, ...)
}
\arguments{
  \item{object}{An object of class \code{spwb}, \code{pwb}, \code{growth} or \code{fordyn}.}
  \item{freq}{Frequency of summary statistics (see \code{\link{cut.Date}}).}
  \item{output}{The data table to be summarized. Accepted values are \code{"DailyBalance"}, \code{"PlantStress"}, \code{"PlantPsi"}, \code{"PlantTranspiration"}, \code{"PlantPhotosynthesis"}  and \code{"SoilWaterBalance"}, \code{"Temperature"} and \code{"EnergyBalance"}.}
  \item{FUN}{The function to summarize results (e.g., \code{sum}, \code{mean}, ...)}
  \item{bySpecies}{Allows aggregating output by species before calculating summaries (only has an effect with some values of \code{output}). Aggregation can involve a sum (as for plant lai or transpiration) or a LAI-weighted mean (as for plant stress or plant water potential).}
  \item{...}{Additional parameters for function \code{summary}.}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{
 \code{\link{spwb}}
}
\examples{
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforestMED)

#Default species parameterization
data(SpParamsMED)

#Initialize soil with default soil params (2 layers)
examplesoil = soil(defaultSoilParams(2))

#Initialize control parameters
control = defaultControl("Granier")

#Initialize input
x = forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)

#Call simulation function
S1<-spwb(x, examplemeteo, latitude = 41.82592, elevation = 100)
        
#Monthly summary (averages) of soil status
summary(S1, freq="months",FUN=mean, output="Soil")

}
