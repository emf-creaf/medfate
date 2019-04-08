\encoding{UTF-8}
\name{hydrology_verticalInputs}
\alias{hydrology_verticalInputs}
\title{Soil water input processes}
\description{Simulates soil water input processes, including canopy water interception, snow accumulation/melt, soil infiltration and percolation.}
\usage{
hydrology_verticalInputs(soil, soilFunctions, prec, er, tday, rad,
                         elevation, Cm, LgroundPAR, LgroundSWR, runon = 0,
                         snowpack = TRUE, drainage = TRUE, modifySoil = TRUE)
}
\arguments{
  \item{soil}{A list containing the description of the soil (see \code{\link{soil}}).}
  \item{soilFunctions}{Soil water retention curve and conductivity functions, either 'SX' (for Saxton) or 'VG' (for Van Genuchten).}
  \item{prec}{Precipitation for a given day (mm)}
  \item{er}{The ratio of evaporation rate to rainfall rate.}
  \item{tday}{Average day temperature (ºC).}
  \item{rad}{Solar radiation (in MJ/m2/day).}
  \item{elevation}{Altitude above sea level (m).}
  \item{Cm}{Canopy water storage capacity.}
  \item{LgroundPAR}{Proportion of photosynthetically-acvive radiation (PAR) reaching the ground.}
  \item{LgroundSWR}{Proportion of short-wave radiation (SWR) reaching the ground.}
  \item{runon}{Surface water amount running on the target area from upslope (in mm).}
  \item{snowpack}{Boolean flag to indicate the simulation of snow accumulation and melting.}
  \item{drainage}{Boolean flag to indicate the simulation of deep drainage.}
  \item{modifySoil}{Boolean flag to indicate that the input \code{soil} object should be modified during the simulation.}
}
\details{
The function simulates different vertical hydrological processes, which are described separately in other functions. If \code{modifySoil = TRUE} the function will modify the \code{soil} object (including both soil moisture and the snowpack on its surface) as a result of simulating hydrological processes.
}
\value{
A named vector with the following elements, all in mm:
\item{Rain}{Precipitation as rainfall.}
\item{Snow}{Precipitation as snow.}
\item{Interception}{Rainfall water intercepted by the canopy and evaporated.}
\item{Throughfall}{Rainfall reaching the ground.}
\item{Snowmelt}{Snow melted during the day, and added to the water infiltrated.}
\item{Runon}{Surface water amount running on the target area from upslope.}
\item{Infiltration}{Water infiltrated into the soil (i.e. throughfall + runon + snowmelt - runoff).}
\item{Runoff}{Surface water leaving the target area.}
\item{DeepDrainage}{Water leaving the target soil towards the water table.}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{
\code{\link{spwb_day}}, \code{\link{hydrology_rainInterception}}, \code{\link{hydrology_soilEvaporation}}
}