\encoding{UTF-8}
\name{spwb_day}
\alias{spwb_day}
\title{Soil-plant water balance for a single day}
\description{
Function \code{spwb_day} performs water balance for a single day.
}
\usage{
spwb_day(x, soil, date, tmin, tmax, rhmin, rhmax, rad, wind, 
        latitude, elevation, slope, aspect, prec, runon = 0.0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{\link{spwbInput}}.}
  \item{soil}{A list containing the description of the soil (see \code{\link{soil}}).}
  \item{date}{Date as string "yyyy-mm-dd".}
  \item{tmin, tmax}{Minimum and maximum temperature (in degrees Celsius).}
  \item{rhmin, rhmax}{Minimum and maximum relative humidity (in percent).}
  \item{rad}{Solar radiation (in MJ/m2/day).}
  \item{wind}{Wind speed (in m/s).}
  \item{prec}{Precipitation (in mm).}
  \item{latitude}{Latitude (in degrees). Required when using the 'Sperry' transpiration mode.}
  \item{elevation, slope, aspect}{Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North). Required when using the 'Sperry' transpiration mode.}
  \item{runon}{Surface water amount running on the target area from upslope (in mm).}
}
\details{
  Detailed model description is available in the vignettes section. The model using 'Granier' transpiration mode is described in De Caceres et al. (2015). Simulations using the 'Sperry' transpiration mode are computationally much more expensive. 
}
\value{
An object (a list) of class \code{spwb_day} with the following elements:

\itemize{
  \item{\code{"cohorts"}: A data frame with cohort information, copied from \code{\link{spwbInput}}.}
  \item{\code{"WaterBalance"}: A vector of water balance components (rain, snow, net rain, infiltration, ...) for the simulated day, equivalent to one row of 'WaterBalance' object given in \code{\link{spwb}}.}
  \item{\code{"Soil"}: A data frame with results for each soil layer:
    \itemize{
      \item{\code{"SoilEvaporation"}: Water evaporated from the soil surface (in mm).}
      \item{\code{"HydraulicInput"}: Water entering each soil layer from other layers, transported via plant hydraulic network (in mm) (only for \code{transpirationMode = "Sperry"}).}
      \item{\code{"HydraulicOutput"}: Water leaving each soil layer (going to other layers or the transpiration stream) (in mm) (only for \code{transpirationMode = "Sperry"}).}
      \item{\code{"PlantExtraction"}: Water extracted by plants from each soil layer (in mm).}
      \item{\code{"psi"}: Soil water potential (in MPa).}
    }
  }
  \item{\code{"EnergyBalance"}: When using the 'Sperry' transpiration mode, the model performs energy balance of the stand and 'EnergyBalance' is a list with the following:
    \itemize{
      \item{\code{"Temperature"}: A data frame with the temperature of the atmosphere ('Tatm'), canopy ('Tcan') and soil ('Tsoil.1', 'Tsoil.2', ...) for each time step.}
      \item{\code{"CanopyEnergyBalance"}: A data frame with the components of the canopy energy balance (in W/m2) for each time step.}
      \item{\code{"SoilEnergyBalance"}: A data frame with the components of the soil energy balance (in W/m2) for each time step.}
    }  
  }
  \item{\code{"Plants"}: A data frame of results for each plant cohort. When using the 'Granier' transpiration mode this includes:
    \itemize{
      \item{\code{"LAI"}: Leaf area index of the plant cohort.}
      \item{\code{"Transpiration"}: Transpirated water (in mm) corresponding to each cohort.}
      \item{\code{"psi"}: Water potential (in MPa) of the plant cohort (average over soil layers).}
      \item{\code{"DDS"}: Daily drought stress [0-1] (relative whole-plant conductance).}
    }
  When using the 'Sperry' transpiration mode the data frame \code{"Plants"} also includes columns:
    \itemize{
      \item{\code{"Extraction"}: Water extracted from the soil (in mm) for each cohort.}
      \item{\code{"RootPsi"}: Minimum water potential (in MPa) at the root collar.}
      \item{\code{"StemPsi"}: Minimum water potential (in MPa) at the stem.}
      \item{\code{"LeafPsi"}: Minimum water potential (in MPa) at the leaf.}
      \item{\code{"StemPLC"}: Proportion of conductance loss in stem.}
      \item{\code{"StemRWC"}: Relative water content of symplastic stem tissue.}
      \item{\code{"LeafRWC"}: Relative water content of symplastic leaf tissue.}
      \item{\code{"dEdP"}: Overall soil-plant conductance (derivative of the supply function).}
    }
  }
  \item{\code{"RhizoPsi"}: Minimum water potential (in MPa) inside roots, after crossing rhizosphere, per cohort and soil layer.}
  \item{\code{"PlantsInst"}: A list with instantaneous (per time step) results for each plant cohort:
    \itemize{
      \item{\code{"LAIsunlit"}: Leaf area index of sunlit leaves of the plant cohort.}
      \item{\code{"LAIshade"}: Leaf area index of shade leaves of the plant cohort.}
      \item{\code{"AbsRad"}: A list with four data frames containing the instantaneous absorved radiation for each plant cohort during each time step. The data frames are combinations of short-wave radiation (SWR) vs long-wave radiation (LWR) and sunlit leaves ('SL') vs. shade leaves ('SH').}
      \item{\code{"E"}: A data frame containing the cumulative transpiration (mm) for each plant cohort during each time step. }
      \item{\code{"An"}: A data frame containing the cumulative net photosynthesis (gC/m2) for each plant cohort during each time step. }
      \item{\code{"GWsunlit"}: A data frame containing instantaneous stomatal conductance (in mmol/m2/s) for sunlit leaves each plant cohort during each time step. }
      \item{\code{"GWshade"}: A data frame containing instantaneous stomatal conductance (in mmol/m2/s) for shade leaves each plant cohort during each time step. }
      \item{\code{"VPDsunlit"}: A data frame containing vapour pressure deficit (in kPa) for sunlit leaves each plant cohort during each time step. }
      \item{\code{"VPDshade"}: A data frame containing vapour pressure deficit (in kPa) for shade leaves each plant cohort during each time step. }
      \item{\code{"Tempsunlit"}: A data frame containing temperature (in degrees Celsius) for sunlit leaves each plant cohort during each time step. }
      \item{\code{"Tempshade"}: A data frame containing temperature (in degrees Celsius) for shade leaves each plant cohort during each time step. }
      \item{\code{"PsiRoot"}: A data frame containing root crown water potential (in MPa) for each plant cohort during each time step.}
      \item{\code{"PsiPlant"}: A data frame containing leaf water potential (in MPa) for each plant cohort during each time step. }
      \item{\code{"PLCstem"}: A data frame containing the proportion loss of conductance [0-1] for each plant cohort during each time step. }
      \item{\code{"RWCstem"}: A data frame containing the (average) relative water content of symplastic stem tissue [0-1] for each plant cohort during each time step. }
      \item{\code{"RWCstem"}: A data frame containing the relative water content of symplastic leaf tissue [0-1] for each plant cohort during each time step. }
    }
  
  }
  
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
\code{\link{spwbInput}}, \code{\link{spwb}},  
\code{\link{spwb_ldrOptimization}}, \code{\link{forest}}
}
\examples{
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Initialize control parameters
control = defaultControl()

#Initialize soil with default soil params (2 layers)
examplesoil = soil(defaultSoilParams(2), W=c(0.5,0.5))


#Simulate one day only
x1 = forest2spwbInput(exampleforest,examplesoil, SpParamsMED, control)
d = 100
sd1<-spwb_day(x1, examplesoil, rownames(examplemeteo)[d],  
             examplemeteo$MinTemperature[d], examplemeteo$MaxTemperature[d], 
             examplemeteo$MinRelativeHumidity[d], examplemeteo$MaxRelativeHumidity[d], 
             examplemeteo$Radiation[d], examplemeteo$WindSpeed[d], 
             latitude = 41.82592, elevation = 100, slope=0, aspect=0,
             prec = examplemeteo$Precipitation[d])

#Initialize soil with default soil params (2 layers)
examplesoil2 = soil(defaultSoilParams(2))

#Switch to 'Sperry' transpiration mode
control$transpirationMode="Sperry"


#Simulate one day only
x2 = forest2spwbInput(exampleforest,examplesoil2, SpParamsMED, control)
d = 100
sd2<-spwb_day(x2, examplesoil2, rownames(examplemeteo)[d],
             examplemeteo$MinTemperature[d], examplemeteo$MaxTemperature[d], 
             examplemeteo$MinRelativeHumidity[d], examplemeteo$MaxRelativeHumidity[d], 
             examplemeteo$Radiation[d], examplemeteo$WindSpeed[d], 
             latitude = 41.82592, elevation = 100, slope=0, aspect=0,
             prec = examplemeteo$Precipitation[d])

#Plot plant transpiration (see function 'plot.swb.day()')
plot(sd2)

}