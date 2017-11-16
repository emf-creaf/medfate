\encoding{UTF-8}
\name{swb}
\alias{swb}
\alias{swb.day}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Soil water balance model}
\description{
Function \code{swb} is a soil water balance model that determines changes in soil moisture, soil water potentials and plant transpiration and drought stress at daily steps for a given forest stand during a period specified in the input climatic data. Additionally, the function also calculates plant net photosynthesis. Transpiration and photosynthesis processes are conducted with different level of detail depending on the transpiration mode. Function \code{swb.day} performs soil water balance for a single day.
}
\usage{
swb(x, soil, meteo,  latitude = NA, elevation = NA, slope = NA, aspect = NA)
swb.day(x, soil, date, doy, tmin, tmax, rhmin, rhmax, rad, wind, 
        latitude, elevation, slope, aspect, rain, runon = 0.0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{\link{swbInput}}.}
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
  \item{date}{Date as string "yyyy-mm-dd".}
  \item{doy}{Day of the year.}
  \item{tmin, tmax}{Minimum and maximum temperature (in degrees Celsius).}
  \item{rhmin, rhmax}{Minimum and maximum relative humidity (in percent).}
  \item{rad}{Solar radiation (in MJ/m2/day).}
  \item{wind}{Wind speed (in m/s).}
  \item{rain}{Precipitation (in mm).}
  \item{runon}{Surface water amount running on the target area from upslope (in mm).}
}
\details{
  Detailed model description is available in the vignettes section. The model using 'Simple' transpiration mode is described in De Caceres et al. (2015). Simulations using the 'Complex' transpiration mode are computationally much more expensive. 
}
\value{
  Function \code{swb} returns a list of class 'swb' with the following elements:
  \itemize{
  \item{\code{"control"}: Control parameters, copied from \code{\link{swbInput}}.}
  \item{\code{"cohorts"}: A data frame with cohort information, copied from \code{\link{swbInput}}.}
  \item{\code{"NumSoilLayers"}: Number of soil layers.}
  \item{\code{"DailyBalance"}: A data frame where different variables (in columns) are given for each simulated day (in rows):}
  \itemize{
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
  \item{\code{"PlantTranspiration"}: A data frame with the amount of daily transpiration (in mm) for each plant cohort. Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
  \item{\code{"PlantPhotosynthesis"}: A data frame with the amount of daily net photosynthesis (in g C·m-2) for each plant cohort. Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
 \item{\code{"PlantPsi"}: A data frame with the average daily water potential of each plant (in MPa). Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
  \item{\code{"PlantStress"}: A data frame with the amount of daily stress suffered by each plant cohort (relative whole-plant conductance). Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
  }
}
\references{
De \enc{Cáceres}{Caceres} M, \enc{Martínez}{Martinez}-Vilalta J, Coll L, Llorens P, Casals P, Poyatos R, Pausas JG, Brotons L. (2015) Coupling a water balance model with forest inventory data to predict drought stress: the role of forest structural changes vs. climate changes. Agricultural and Forest Meteorology (doi:10.1016/j.agrformet.2015.06.012).
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\note{
Objects \code{x} and \code{soil} are modified during the simulation. Daily transpiration and photosynthesis values are stored in columns \code{Transpiration} and \code{Photosynthesis} of object \code{x}. Water content relative to field capacity (vector \code{W}) of \code{soil} is also modified.  
}
\seealso{
\code{\link{swbInput}}, \code{\link{plot.swb}},  \code{\link{swbpoints}}, \code{\link{swbgrid}}, 
\code{\link{swb.ldrOptimization}}, \code{\link{forest}}
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
x1 = forest2swbInput(exampleforest,examplesoil, SpParamsMED, control)

#Call simulation function
S1<-swb(x1, examplesoil, examplemeteo)
        
#Plot results
plot(S1)

#Monthly summary (averages) of soil water balance
summary(S1, freq="months",FUN=mean, output="SoilWaterBalance")



#Initialize soil with default soil params (2 layers)
examplesoil2 = soil(defaultSoilParams(2))

#Switch to 'Complex' transpiration mode
control$transpirationMode="Complex"

#Initialize input
x2 = forest2swbInput(exampleforest,examplesoil2, SpParamsMED, control)

#Call simulation function (5 days)
S2<-swb(x2, examplesoil2, examplemeteo[100:105,], latitude = 41.82592, elevation = 100)


#Simulate one day only
examplesoil2 = soil(defaultSoilParams(2))
x2 = forest2swbInput(exampleforest,examplesoil2, SpParamsMED, control)
d = 1
sd1<-swb.day(x2, examplesoil2, rownames(examplemeteo)[d], d , 
             examplemeteo$MinTemperature[d], examplemeteo$MaxTemperature[d], 
             examplemeteo$MinRelativeHumidity[d], examplemeteo$MaxRelativeHumidity[d], 
             examplemeteo$Radiation[d], examplemeteo$WindSpeed[d], 
             latitude = 41.82592, elevation = 100, 
             slope= 0, aspect = 0, rain = examplemeteo$Precipitation[d])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.