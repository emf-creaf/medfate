\encoding{UTF-8}
\name{spwb}
\alias{spwb}
\alias{pwb}
\alias{spwb_resetInputs}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Soil-plant water balance}
\description{
Function \code{spwb()} is a water balance model that determines changes in soil moisture, soil water potentials and plant transpiration and drought stress at daily steps for a given forest stand during a period specified in the input climatic data. Additionally, the function also calculates plant net photosynthesis. Function \code{pwb()} performs plant water balance only (i.e. soil moisture dynamics is an input) at daily steps for a given forest stand during a period specified in the input climatic data. On both simulation functions transpiration and photosynthesis processes are conducted with different level of detail depending on the transpiration mode.
}
\usage{
spwb(x, soil, meteo,  latitude = NA, elevation = NA, slope = NA, aspect = NA)
pwb(x, soil, meteo, W, latitude = NA, elevation = NA, 
    canopyEvaporation = numeric(0), 
    snowMelt = numeric(0), 
    soilEvaporation = numeric(0))
spwb_resetInputs(x, soil, from = NULL, day = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{\link{spwbInput}}.}
  \item{soil}{A list containing the description of the soil (see \code{\link{soil}}).}
  \item{meteo}{A data frame with daily meteorological data series. Row names of the data frame should correspond to date strings with format "yyyy-mm-dd" (see \code{\link{Date}}). When using the 'Simple' transpiration mode the following columns are required:
   \itemize{
     \item{\code{MeanTemperature}: Mean temperature (in degrees Celsius).}
     \item{\code{Precipitation}: Precipitation (in mm).}
     \item{\code{Radiation}: Solar radiation (in MJ/m2/day), required only if \code{snowpack = TRUE}.}
     \item{\code{PET}: Potential evapotranspiration (in mm).}
    }
  When using the 'Complex' transpiration mode the following columns are required:
   \itemize{
     \item{\code{MeanTemperature}: Mean temperature (in degrees Celsius).}
     \item{\code{MinTemperature}: Minimum temperature (in degrees Celsius).}
     \item{\code{MaxTemperature}: Maximum temperature (in degrees Celsius).}
     \item{\code{MinRelativeHumidity}: Minimum relative humidity (in percent).}
     \item{\code{MaxRelativeHumidity}: Maximum relative humidity (in percent).}
     \item{\code{Precipitation}: Precipitation (in mm).}
     \item{\code{Radiation}: Solar radiation (in MJ/m2/day).}
     \item{\code{WindSpeed}: Wind speed (in m/s). If not available, this column can be left with \code{NA} values.}
    }
  }
  \item{W}{A matrix with the same number of rows as \code{meteo} and as many columns as soil layers, containing the soil moisture of each layer as proportion of field capacity.}
  \item{latitude}{Latitude (in degrees). Required when using the 'Complex' transpiration mode.}
  \item{elevation, slope, aspect}{Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North). Required when using the 'Complex' transpiration mode. Elevation is also required for 'Simple' if snowpack dynamics are simulated.}
  \item{from}{An object of class \code{spwb} storing the results of a previous simulation, including values of state variables. If \code{from = NULL}, state variables are set to their defaults (i.e. soil moisture set to field capacity and cumulative growth degree days set to zero).}
  \item{day}{An integer with the day from which state variable values stored in \code{from} should be copied. If missing, values are copied from the first day of stored values.}
}
\details{
  Detailed model description is available in the vignettes section. The model using 'Simple' transpiration mode is illustrated by function \code{\link{transp_transpirationGranier}} and described in De Caceres et al. (2015). Simulations using the 'Complex' transpiration mode are computationally much more expensive and are illustrated by function \code{\link{transp_transpirationSperry}}. Function \code{spwb_resetInputs()} allows resetting state variables in \code{x} and \code{soil} to their defaults, or to copy values of state variables from a previous \code{spwb()} simulation results stored in \code{from}.
}
\value{
  Function \code{spwb} returns a list of class 'spwb' whereas Function \code{pwb} returns a list of class 'pwb'. There are many elements in common in these lists, so they are listed here together:
  \itemize{
     \item{\code{"spwbInput"}: A copy of the object \code{x} of class \code{\link{spwbInput}} given as input.}
     \item{\code{"soilInput"}: A copy of the object \code{soil} of class \code{\link{soil}} given as input.}
     \item{\code{"WaterBalance"}: A data frame where different variables (in columns) are given for each simulated day (in rows):}
       \itemize{
          \item{\code{"GDD"}: Cumulative growth degree days.}
          \item{\code{"LAIcell"}: The LAI of the stand (accounting for leaf phenology) (in m2/m2).}
          \item{\code{"Cm"}: The water retention capacity of the stand (in mm) (accounting for leaf phenology).}
          \item{\code{"Lground"}: The proportion of light that reaches the ground (accounting for leaf phenology).}
          \item{\code{"PET"}: Potential evapotranspiration (in mm).}
          \item{\code{"Precipitation"}: Input precipitation (in mm).}
          \item{\code{"Rain"}: Precipitation as rain (in mm).}
          \item{\code{"Snow"}: Precipitation as snow (in mm).}
          \item{\code{"NetRain"}: Net rain, after accounting for interception (in mm).}
          \item{\code{"Infiltration"}: The amount of water infiltrating into the soil (in mm).}
          \item{\code{"Runoff"}: The amount of water exported via surface runoff (in mm).}
          \item{\code{"DeepDrainage"}: The amount of water exported via deep drainage (in mm).}
          \item{\code{"Evapotranspiration"}: Evapotranspiration (in mm).}
          \item{\code{"SoilEvaporation"}: Bare soil evaporation (in mm).}
          \item{\code{"PlantExtraction"}: Amount of water extracted from soil by plants (in mm) (can only be diferent from transpiration for \code{transpirationMode = "Complex"} when capacitance is considered).}
          \item{\code{"Transpiration"}: Plant transpiration (considering all soil layers) (in mm).}
          \item{\code{"HydraulicRedistribution"}: Water redistributed among soil layers, transported through the plant hydraulic network (only for \code{transpirationMode = "Complex"}).}
          
      }
    \item{\code{"Soil"}: A data frame where different variables (in columns) are given for each simulated day (in rows):}
        \itemize{
          \item{\code{"W.1"}, \code{...}, \code{"W.k"}: Relative soil moisture content (relative to field capacity) in each soil layer.}
          \item{\code{"ML.1"}, \code{...}, \code{"ML.k"}: Soil water volume in each soil layer (in L/m2).}
          \item{\code{"MLTot"}: Total soil water volume (in L/m2).}
          \item{\code{"SWE"}: Snow water equivalent (mm) of the snow pack.}
          \item{\code{"PlantExt.1"}, \code{...}, \code{"PlantExt.k"}: Plant extraction from each soil layer (in mm).}
          \item{\code{"HydraulicInput.1"}, \code{...}, \code{"HydraulicInput.k"}: Water that entered the layer coming from other layers and transported via the plant hydraulic network (in mm) (only for \code{transpirationMode = "Complex"}).}
          \item{\code{"psi.1"}, \code{...}, \code{"psi.k"}: Soil water potential in each soil layer (in MPa).}
        }
     \item{\code{"PlantLAI"}: A data frame with the daily leaf area index for each plant cohort. Days are in rows and plant cohorts are in columns.}
     \item{\code{"PlantTranspiration"}: A data frame with the amount of daily transpiration (in mm) for each plant cohort. Days are in rows and plant cohorts are in columns.}
     \item{\code{"PlantPhotosynthesis"}: A data frame with the amount of daily net photosynthesis (in g C·m-2) for each plant cohort. Days are in rows and plant cohorts are in columns.}
     \item{\code{"PlantPsi"}: A data frame with the average daily water potential of each plant (in MPa). Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
     \item{\code{"PlantStress"}: A data frame with the amount of daily stress [0-1] suffered by each plant cohort (relative whole-plant conductance). Days are in rows and plant cohorts are in columns.}
     \item{\code{"subdaily"}: A list of objects of class \code{\link{spwb_day}}, one per day simulated (only if required in \code{control} parameters, see \code{\link{defaultControl}}).}
  }
  If \code{transpirationMode=="Complex"} the list also includes the following elements:
 \itemize{
   \item{\code{"PlantAbsorbedSWR"}: A data frame with the daily SWR absorbed by each plant cohort. Days are in rows and plant cohorts are in columns.}
  \item{\code{"PlantAbsorbedLWR"}: A data frame with the daily LWR absorbed by each plant cohort. Days are in rows and plant cohorts are in columns.}
  \item{\code{"dEdP"}: A data frame with mean daily values of soil-plant conductance (derivative of the supply function) for each plant cohort. Days are in rows and plant cohorts in columns.}
  \item{\code{"EnergyBalance"}: A data frame with the daily values of energy balance components for the soil and the canopy.}
 \item{\code{"Temperature"}: A data frame with the daily values of minimum/mean/maximum temperatures for the atmosphere (input), canopy and soil.}
 \item{\code{"LeafPsi"}: A data frame with the minimum daily leaf water potential of each plant (in MPa). Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
 \item{\code{"LeafRWC"}: A data frame with the average daily leaf relative water content of each plant (in percent). Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
 \item{\code{"StemPsi"}: A data frame with the minimum daily stem water potential of each plant (in MPa). Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
 \item{\code{"StemPLC"}: A data frame with the average daily proportion of stem conductance loss of each plant ([0-1]). Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
 \item{\code{"StemRWC"}: A data frame with the average daily stem relative water content of each plant (in percent). Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
 \item{\code{"RootPsi"}: A data frame with the minimum daily root water potential of each plant (in MPa). Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
 \item{\code{"RhizoPsi"}: A list of data frames (one per plant cohort) with the minimum daily root water potential of each plant (in MPa). In each data frame, days are in rows and layers are in columns.}
 }
}
\references{
De \enc{Cáceres}{Caceres} M, \enc{Martínez}{Martinez}-Vilalta J, Coll L, Llorens P, Casals P, Poyatos R, Pausas JG, Brotons L. (2015) Coupling a water balance model with forest inventory data to predict drought stress: the role of forest structural changes vs. climate changes. Agricultural and Forest Meteorology 213: 77-90 (doi:10.1016/j.agrformet.2015.06.012).
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\note{
State variables stored in objects \code{x} and \code{soil} are modified during the simulation. Use function \code{spwb_resetInputs()} to reset state variables to defaults. Daily transpiration and photosynthesis values are stored in columns \code{Transpiration} and \code{Photosynthesis} of object \code{x}. Water content relative to field capacity (vector \code{W}) of \code{soil} is also modified.  
}
\seealso{
\code{\link{spwbInput}}, \code{\link{spwb_day}}, \code{\link{plot.spwb}},  
\code{\link{spwb_ldrOptimization}}, \code{\link{forest}}
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
S1<-spwb(x1, examplesoil, examplemeteo, elevation = 100)
        
#Plot results
plot(S1)

#Monthly summary (averages) of soil water balance
summary(S1, freq="months",FUN=mean, output="Soil")



#Initialize soil with default soil params (2 layers)
examplesoil2 = soil(defaultSoilParams(2))

#Switch to 'Complex' transpiration mode
control$transpirationMode="Complex"

#Initialize input
x2 = forest2spwbInput(exampleforest,examplesoil2, SpParamsMED, control)

#Call simulation function (5 days)
S2<-spwb(x2, examplesoil2, examplemeteo[100:105,], latitude = 41.82592, elevation = 100)

}