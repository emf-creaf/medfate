\encoding{UTF-8}
\name{exampleobs}
\Rdversion{1.1}
\alias{exampleobs}
\docType{data}
\title{Example observed data}
\description{
  Example (fake) data set of variables measured in a plot. 
}
\usage{data(exampleobs)}
\format{
 A data frame containing daily 'observed' values for year 2001.
   \describe{
\item{\code{SWC}}{Soil moisture content (in m3/m3).}
\item{\code{ETR}}{Total evapotranspiration (mm).}
\item{\code{E_T1_148}}{Transpiration of Pinus halepensis cohort 'T1_148' (L/m2 of leaf area).}
\item{\code{E_T2_168}}{Transpiration of Quercus ilex cohort 'T2_168' (L/m2 of leaf area).}
\item{\code{FMC_T1_148}}{Fuel moisture content of Pinus halepensis cohort 'T1_148'  (in percent).}
\item{\code{FMC_T2_168}}{Fuel moisture content of Quercus ilex cohort 'T2_168' (in percent).}
\item{\code{BAI_T1_148}}{Basal area increment for Pinus halepensis cohort 'T1_148'  (in cm2).}
\item{\code{BAI_T2_168}}{Basal area increment for Quercus ilex cohort 'T2_168' (in cm2).}
\item{\code{DI_T1_148}}{Diameter increment for Pinus halepensis cohort 'T1_148'  (in cm).}
\item{\code{DI_T2_168}}{Diameter increment for Quercus ilex cohort 'T2_168' (in cm).}
  }
}
\source{
This data set was actually created by running a simulation and adding some gaussian error to the outputs.
}
\seealso{ \code{\link{evaluation}}}
\examples{
data(exampleobs)
}
\keyword{datasets}
