\encoding{UTF-8}
\name{spwb}
\alias{spwb}
\alias{spwb.resetInputs}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Soil-plant water balance}
\description{
Function \code{spwb()} is a water balance model that determines changes in soil moisture, soil water potentials and plant transpiration and drought stress at daily steps for a given forest stand during a period specified in the input climatic data. Additionally, the function also calculates plant net photosynthesis. Transpiration and photosynthesis processes are conducted with different level of detail depending on the transpiration mode.
}
\usage{
spwb(x, soil, meteo,  latitude = NA, elevation = NA, slope = NA, aspect = NA)
spwb.resetInputs(x, soil, from = NULL, day = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{\link{spwbInput}}.}
  \item{soil}{A list containing the description of the soil (see \code{\link{soil}}).}
  \item{meteo}{A data frame with daily meteorological data series. When using the 'Simple' transpiration mode the following columns are required:
   \itemize{
     \item{\code{DOY}: Day of the year.}
     \item{\code{Precipitation}: Precipitation (in mm).}
     \item{\code{MeanTemperature}: Mean temperature (in degrees Celsius).}
     \item{\code{PET}: Potential evapotranspiration (in mm).}
    }
  When using the 'Complex' transpiration mode the following columns are required:
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
  \item{latitude}{Latitude (in degrees). Required when using the 'Complex' transpiration mode.}
  \item{elevation, slope, aspect}{Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North). Required when using the 'Complex' transpiration mode.}
  \item{from}{An object of class \code{spwb} storing the results of a previous simulation, including values of state variables. If \code{from = NULL}, state variables are set to their defaults (i.e. soil moisture set to field capacity and cumulative growth degree days set to zero).}
  \item{day}{An integer with the day from which state variable values stored in \code{from} should be copied. If missing, values are copied from the first day of stored values.}
}
\details{
  Detailed model description is available in the vignettes section. The model using 'Simple' transpiration mode is described in De Caceres et al. (2015). Simulations using the 'Complex' transpiration mode are computationally much more expensive. Function \code{spwb.resetInputs()} allows resetting state variables in \code{x} and \code{soil} to their defaults, or to copy values of state variables from a previous \code{spwb()} simulation results stored in \code{from}.
}
\value{
  Function \code{spwb} returns a list of class 'spwb' with the following elements (both transpiration modes):
  \itemize{
  \item{\code{"control"}: Control parameters, copied from \code{\link{spwbInput}}.}
  \item{\code{"cohorts"}: A data frame with cohort information, copied from \code{\link{spwbInput}}.}
  \item{\code{"NumSoilLayers"}: Number of soil layers.}
  \item{\code{"DailyBalance"}: A data frame where different variables (in columns) are given for each simulated day (in rows):}
  \itemize{
  \item{\code{"GDD"}: Cumulative growth degree days.}
  \item{\code{"LAIcell"}: The LAI of the stand (accounting for leaf phenology) (in m2/m2).}
  \item{\code{"Cm"}: The water retention capacity of the stand (in mm) (accounting for leaf phenology).}
  \item{\code{"Lground"}: The proportion of light that reaches the ground (accounting for leaf phenology).}
  \item{\code{"PET"}: Potential evapotranspiration (in mm).}
  \item{\code{"Rainfall"}: Input precipitation (in mm).}
  \item{\code{"NetPrec"}: Net precipitation, after accounting for interception (in mm).}
  \item{\code{"Infiltration"}: The amount of water infiltrating into the soil (in mm).}
  \item{\code{"Runoff"}: The amount of water exported via surface runoff (in mm).}
  \item{\code{"DeepDrainage"}: The amount of water exported via deep drainage (in mm).}
  \item{\code{"Etot"}: Evapotranspiration (in mm).}
  \item{\code{"Esoil"}: Bare soil evaporation (in mm).}
  \item{\code{"Eplanttot"}: Plant transpiration (considering all soil layers) (in mm).}
  \item{\code{"Eplant.1"}, \code{...}, \code{"Eplant.k"}: Plant transpiration from each soil layer (in mm).}
  }
  \item{\code{"SoilWaterBalance"}: A data frame where different variables (in columns) are given for each simulated day (in rows):}
  \itemize{
  \item{\code{"W.1"}, \code{...}, \code{"W.k"}: Relative soil moisture content (relative to field capacity) in each soil layer.}
  \item{\code{"ML.1"}, \code{...}, \code{"ML.k"}: Soil water volume in each soil layer (in L/m2).}
  \item{\code{"MLTot"}: Total soil water volume (in L/m2).}
  \item{\code{"psi.1"}, \code{...}, \code{"psi.k"}: Soil water potential in each soil layer (in MPa).}
  }
  \item{\code{"PlantLAI"}: A data frame with the daily leaf area index for each plant cohort. Days are in rows and plant cohorts are in columns.}
  \item{\code{"PlantTranspiration"}: A data frame with the amount of daily transpiration (in mm) for each plant cohort. Days are in rows and plant cohorts are in columns.}
  \item{\code{"PlantPhotosynthesis"}: A data frame with the amount of daily net photosynthesis (in g C·m-2) for each plant cohort. Days are in rows and plant cohorts are in columns.}
 \item{\code{"PlantPsi"}: A data frame with the average daily water potential of each plant (in MPa). Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
  \item{\code{"PlantStress"}: A data frame with the amount of daily stress suffered by each plant cohort (relative whole-plant conductance). Days are in rows and plant cohorts are in columns.}
  }
  If \code{transpirationMode=="Complex"} the list also includes the following elements:
 \itemize{
   \item{\code{"PlantAbsorbedSWR"}: A data frame with the daily SWR absorbed by each plant cohort. Days are in rows and plant cohorts are in columns.}
  \item{\code{"PlantAbsorbedLWR"}: A data frame with the daily LWR absorbed by each plant cohort. Days are in rows and plant cohorts are in columns.}
  \item{\code{"EnergyBalance"}: A data frame with the daily values of energy balance components for the soil and the canopy.}
 \item{\code{"Temperature"}: A data frame with the daily values of minimum/mean/maximum temperatures for the atmosphere (input), canopy and soil.}
 }
}
\references{
De \enc{Cáceres}{Caceres} M, \enc{Martínez}{Martinez}-Vilalta J, Coll L, Llorens P, Casals P, Poyatos R, Pausas JG, Brotons L. (2015) Coupling a water balance model with forest inventory data to predict drought stress: the role of forest structural changes vs. climate changes. Agricultural and Forest Meteorology 213: 77-90 (doi:10.1016/j.agrformet.2015.06.012).
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\note{
State variables stored in objects \code{x} and \code{soil} are modified during the simulation. Use function \code{spwb.resetInputs()} to reset state variables to defaults. Daily transpiration and photosynthesis values are stored in columns \code{Transpiration} and \code{Photosynthesis} of object \code{x}. Water content relative to field capacity (vector \code{W}) of \code{soil} is also modified.  
}
\seealso{
\code{\link{spwbInput}}, \code{\link{spwb.day}}, \code{\link{plot.spwb}},  \code{\link{spwbpoints}}, \code{\link{spwbgrid}}, 
\code{\link{spwb.ldrOptimization}}, \code{\link{forest}}
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
x1 = forest2spwbInput(exampleforest,examplesoil, SpParamsMED, control)

#Call simulation function
S1<-spwb(x1, examplesoil, examplemeteo)
        
#Plot results
plot(S1)

#Monthly summary (averages) of soil water balance
summary(S1, freq="months",FUN=mean, output="SoilWaterBalance")



#Initialize soil with default soil params (2 layers)
examplesoil2 = soil(defaultSoilParams(2))

#Switch to 'Complex' transpiration mode
control$transpirationMode="Complex"

#Initialize input
x2 = forest2spwbInput(exampleforest,examplesoil2, SpParamsMED, control)

#Call simulation function (5 days)
S2<-spwb(x2, examplesoil2, examplemeteo[100:105,], latitude = 41.82592, elevation = 100)

}