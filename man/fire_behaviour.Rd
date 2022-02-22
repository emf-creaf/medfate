\encoding{UTF-8}
\name{fire_behaviour}
\alias{fire_Rothermel}
\alias{fire_FCCS}
\title{
Fire behaviour functions
}
\description{
Function \code{fire_FCCS()} implements a modification of the fire behavior models described for the Fuel Characteristics Classification System (FCCS) in Prichard et al. (2013). Function \code{fire_Rothermel()} implements Rothermel's (1972) fire behaviour model (modified from package 'Rothermel' (Giorgio Vacchiano, Davide Ascoli)).


}
\usage{
fire_FCCS(FCCSpropsSI, MliveSI = as.numeric(c(90, 90, 60)), 
          MdeadSI = as.numeric(c(6, 6, 6, 6, 6)), 
          slope = 0, windSpeedSI = 11)
fire_Rothermel(modeltype, wSI, sSI, delta, mx_dead,
               hSI, mSI, u, windDir, slope, aspect)
                
}
\arguments{
  \item{FCCSpropsSI}{A data frame describing the properties of five fuel strata (canopy, shrub, herbs, dead woody and litter) returned by \code{\link{fuel_FCCS}}.}
  \item{MliveSI}{Moisture of live fuels (in percent of dry weight) for canopy, shrub, and herb strata. Live moisture values are drawn from column \code{ActFCM} in \code{FCCSpropsSI} if available (see \code{\link{fuel_FCCS}}). Otherwise, moisture values supplied for \code{MliveSI} are used.}
  \item{MdeadSI}{Moisture of dead fuels (in percent of dry weight) for canopy, shrub, herb, woody and litter strata.}
  \item{slope}{Slope (in degrees).}
  \item{windSpeedSI}{Wind speed (in m/s) at 20 ft (6 m) over vegetation (default 11 m/s = 40 km/h)}
  \item{modeltype}{'S'(tatic) or 'D'(ynamic)}
  \item{wSI}{A vector of fuel load (t/ha) for five fuel classes.}
  \item{sSI}{A vector of surface-to-volume ratio (m2/m3) for five fuel classes.}
  \item{delta}{A value of fuel bed depth (cm).}
  \item{mx_dead}{A value of dead fuel moisture of extinction (percent).}
  \item{hSI}{A vector of heat content (kJ/kg) for five fuel classes.}
  \item{mSI}{A vector of percent moisture on a dry weight basis (percent) for five fuel classes.}
 \item{u}{A value of windspeed (m/s) at midflame height.}
 \item{windDir}{Wind direction (in degrees from north). North means blowing from north to south.}
 \item{aspect}{Aspect (in degrees from north).}
}
\details{
Default moisture, slope and windspeed values are benchmark conditions used to calculate fire potentials (Sandberg et al. 2007) and map vulnerability to fire.
}
\value{
Both functions return list with fire behavior variables. In the case of \code{fire_FCCS}, the function returns the variables in three blocks (lists \code{SurfaceFire}, \code{CrownFire} and \code{FirePotentials}), and the values are:
\itemize{
\item{\code{SurfaceFire$`midflame_WindSpeed [m/s]`}: Midflame wind speed in the surface fire.}
\item{\code{SurfaceFire$phi_wind}: Spread rate modifier due to wind.}
\item{\code{SurfaceFire$phi_slope}: Spread rate modifier due to slope.}
\item{\code{SurfaceFire$`I_R_surf [kJ/m2/min]`}: Intensity of the surface fire reaction.}
\item{\code{SurfaceFire$`I_R_litter [kJ/m2/min]`}: Intensity of the litter fire reaction.}
\item{\code{SurfaceFire$`q_surf [kJ/m2]`}: Heat sink of the surface fire.}
\item{\code{SurfaceFire$`q_litter [kJ/m2]`}: Heat sink of the litter fire.}
\item{\code{SurfaceFire$xi_surf}: Propagating flux ratio of the surface fire.}
\item{\code{SurfaceFire$xi_litter}: Propagating flux ratio of the litter fire.}
\item{\code{SurfaceFire$`ROS_surf [m/min]`}: Spread rate of the surface fire(without accounting for faster spread in the litter layer).}
\item{\code{SurfaceFire$`ROS_litter [m/min]`}: Spread rate of the litter fire.}
\item{\code{SurfaceFire$`ROS_windslopecap [m/min]`}: Maximum surface fire spread rate according to wind speed.}
\item{\code{SurfaceFire$`ROS [m/min]`}: Final spread rate of the surface fire.}
\item{\code{SurfaceFire$`I_b [kW/m]`}: Fireline intensity of the surface fire.}
\item{\code{SurfaceFire$`FL [m]`}: Flame length of the surface fire.}
\item{\code{CrownFire$`I_R_canopy [kJ/m2/min]`}: Intensity of the canopy fire reaction.}
\item{\code{CrownFire$`I_R_crown [kJ/m2/min]`}: Intensity of the crown fire reaction (adding surface and canopy reactions).}
\item{\code{CrownFire$`q_canopy [kJ/m2]`}: Heat sink of the canopy fire.}
\item{\code{CrownFire$`q_crown [kJ/m2]`}: Heat sink of the crown fire (adding surface and canopy heat sinks).}
\item{\code{CrownFire$xi_surf}: Propagating flux ratio of the crown fire.}
\item{\code{CrownFire$`canopy_WindSpeed [m/s]`}: Wind speed in the canopy fire (canopy top wind speed).}
\item{\code{CrownFire$WAF}: Wind speed adjustment factor for crown fires.}
\item{\code{CrownFire$`ROS [m/min]`}: Spread rate of the crown fire.}
\item{\code{CrownFire$Ic_ratio}: Crown initiation ratio.}
\item{\code{CrownFire$`I_b [kW/m]`}: Fireline intensity of the crown fire.}
\item{\code{CrownFire$`FL [m]`}: Flame length of the crown fire.}
\item{\code{FirePotentials$RP}: Surface fire reaction potential ([0-9]).}
\item{\code{FirePotentials$SP}: Surface fire spread rate potential ([0-9]).}
\item{\code{FirePotentials$FP}: Surface fire flame length potential ([0-9]).}
\item{\code{FirePotentials$SFP}: Surface fire potential ([0-9]).
}
\item{\code{FirePotentials$IC}: Crown initiation potential ([0-9]).}
\item{\code{FirePotentials$TC}: Crown-to-crown transmission potential ([0-9]).}
\item{\code{FirePotentials$RC}: Crown fire spread rate potential ([0-9]).}
\item{\code{FirePotentials$CFC}: Crown fire potential ([0-9]).}
}
}
\references{
Albini, F. A. (1976). Computer-based models of wildland fire behavior: A users' manual. Ogden, UT: US Department of Agriculture, Forest Service, Intermountain Forest and Range Experiment Station.

Rothermel, R. C. 1972. A mathematical model for predicting fire spread in wildland fuels. USDA Forest Service Research Paper INT USA.

Prichard, S. J., D. V Sandberg, R. D. Ottmar, E. Eberhardt, A. Andreu, P. Eagle, and K. Swedin. 2013. Classification System Version 3.0: Technical Documentation.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\note{
Default moisture, slope and windspeed values are benchmark conditions used to calculate fire potentials (Sandberg et al. 2007) and map vulnerability to fire.
}
\seealso{
\code{\link{fuel_FCCS}}
}
\examples{
#Load example plot plant data
data(exampleforestMED)

#Default species parameterization
data(SpParamsMED)

#Calculate fuel properties according to FCCS
fccs = fuel_FCCS(exampleforestMED, SpParamsMED)

#Calculate fire behavior according to FCCS
fire_FCCS(fccs)

#Load fuel model parameter data
data(SFM_metric)

#Fuel stratification (returns heights in cm)
fs = fuel_stratification(exampleforestMED, SpParamsMED)

#Correct windspeed (transform heights to m)
u = 11 #m/s
umf = u*fuel_windAdjustmentFactor(fs$surfaceLayerTopHeight/100, fs$canopyBaseHeight/100, 
                                  fs$canopyTopHeight/100, 60)

#Call Rothermel function using fuel model 'A6'
fire_Rothermel(modeltype="D", wSI = as.numeric(SFM_metric["A6",2:6]), 
                    sSI = as.numeric(SFM_metric["A6",7:11]), 
                    delta = as.numeric(SFM_metric["A6",12]),
                    mx_dead = as.numeric(SFM_metric["A6",13]),
                    hSI = as.numeric(SFM_metric["A6",14:18]),
                    mSI = c(10,10,10,30,60),
                    u=umf, windDir=0, slope=0, aspect=0)
}
