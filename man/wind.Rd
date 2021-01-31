\encoding{UTF-8}
\name{Wind models}
\alias{wind_canopyTurbulence}
\title{
Models for canopy turbulence
}
\description{
Models for canopy turbulence by Katul et al (2004).
}
\usage{
wind_canopyTurbulence(z, Cx, h, d0, z0, model = "k-epsilon")
}
\arguments{
  \item{z}{A numeric vector with height values.}
  \item{Cx}{Effective drag = Cd x leaf area density.}
  \item{h}{Canopy height (m).}
  \item{d0}{Zero displacement height (m).}
  \item{z0}{Momentum roughness height (m).}
}
\value{ A data frame of vertical profiles for variables:
\itemize{
\item{\code{z1}: Height values.}
\item{\code{U1}: U/u*, where U is mean velocity and u* is friction velocity.}
\item{\code{epsilon1}: epsilon/(u*3/h) where epsilon is the turbulent kinetic dissipation rate, u* is friction velocity and h is canopy height.}
\item{\code{k1}: k/(u*^2), where k is the turbulent kinetic energy and u* is friction velocity.}
\item{\code{uw1}: <uw>/(u*^2), where <uw> is the Reynolds stress and u* is friction velocity.}
\item{\code{Lmix1}: Mixing length.}
}
}
\details{
 Implementation in Rcpp of the models by Katul et al (2004) originally in Matlab code.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
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
