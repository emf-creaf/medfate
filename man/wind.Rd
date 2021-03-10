\encoding{UTF-8}
\name{Wind models}
\alias{wind_canopyTurbulenceModel}
\alias{wind_canopyTurbulence}
\title{
Models for canopy turbulence
}
\description{
Models for canopy turbulence by Katul et al (2004).
}
\usage{
wind_canopyTurbulenceModel(zm, Cx, hm, d0, z0, model = "k-epsilon")
wind_canopyTurbulence(zmid, LAD, canopyHeight, 
                      u, windMeasurementHeight = 200, model = "k-epsilon")
}
\arguments{
  \item{zm}{A numeric vector with height values (m).}
  \item{Cx}{Effective drag = Cd x leaf area density.}
  \item{hm}{Canopy height (m).}
  \item{d0}{Zero displacement height (m).}
  \item{z0}{Momentum roughness height (m).}
  \item{zmid}{A numeric vector of mid-point heights (in cm) for canopy layers.}
  \item{LAD}{A numeric vector of leaf area density values (m3/m2).}
  \item{canopyHeight}{Canopy height (in cm).}
  \item{u}{Measured wind speed (m/s).}
  \item{windMeasurementHeight}{Height of wind speed measurement with respect to canopy height (cm).}
  \item{model}{Closure model.}
}
\value{ Function \code{wind_canopyTurbulenceModel} returns a data frame of vertical profiles for variables:
\itemize{
\item{\code{z1}: Height values.}
\item{\code{U1}: U/u*, where U is mean velocity and u* is friction velocity.}
\item{\code{dU1}: dUdz/u*, where dUdz is mean velocity gradient and u* is friction velocity.}
\item{\code{epsilon1}: epsilon/(u^3/h) where epsilon is the turbulent kinetic dissipation rate, u* is friction velocity and h is canopy height.}
\item{\code{k1}: k/(u*^2), where k is the turbulent kinetic energy and u* is friction velocity.}
\item{\code{uw1}: <uw>/(u*^2), where <uw> is the Reynolds stress and u* is friction velocity.}
\item{\code{Lmix1}: Mixing length.}
}
Function \code{wind_canopyTurbulence} returns a data frame of vertical profiles for transformed variables:
\itemize{
\item{\code{zmid}: Input mid-point heights (in cm) for canopy layers.}
\item{\code{u}: Wind speed (m/s).}
\item{\code{du}: Mean velocity gradient (1/s).}
\item{\code{epsilon}: Turbulent kinetic dissipation rate.}
\item{\code{k}: Turbulent kinetic energy.}
\item{\code{uw}: Reynolds stress.}
}
}
\details{
 Implementation in Rcpp of the K-epsilon canopy turbulence models by Katul et al (2004) originally in Matlab code (https://nicholas.duke.edu/people/faculty/katul/k_epsilon_model.htm).
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\references{
Katul GG, Mahrt L, Poggi D, Sanz C (2004) One- and two-equation models for canopy turbulence. Boundary-Layer Meteorol 113:81–109. https://doi.org/10.1023/B:BOUN.0000037333.48760.e5
}
\seealso{
\code{\link{vprofile_windExtinction}}
}
\examples{
#Default species parameterization
data(SpParamsMED)

#Load example plot plant data
data(exampleforestMED)

#Canopy height (in m)
h= max(exampleforestMED$treeData$Height/100) 
d0 = 0.67*h
z0 = 0.08*h

#Height values (cm)
z = seq(50,1000, by=50)
zm = z/100 # (in m)

# Leaf area density
lad = vprofile_leafAreaDensity(exampleforestMED, SpParamsMED, draw = FALSE, z = z)

# Effective drag
Cd = 0.2
Cx = Cd*lad

# canopy turbulence model
wind_canopyTurbulence(zm, Cx,h,d0,z0)
}
