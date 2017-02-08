\name{SFM_metric}
\alias{SFM_metric}
\docType{data}
\title{
Standard fuel models (Albini 1976, Scott & Burgan 2005)
}
\description{
Standard fuel models converted to metric system. Copied from package 'Rothermel' (Giorgio Vacchiano, Davide Ascoli).
}
\usage{data("SFM_metric")}
\format{
A data frame including standard fuel models as in Albini (1976) and Scott and Burgan (2005), to be used as input of \code{\link{fire.Rothermel}} function. All values converted to metric format.
  \describe{
    \item{\code{Fuel_Model_Type}}{A factor with levels \code{D} (for dynamic) or \code{S} (for static).}
    \item{\code{Load_1h}}{Loading of 1h fuel class [t/ha].}
    \item{\code{Load_10h}}{Loading of 10h fuel class [t/ha].}
    \item{\code{Load_100h}}{Loading of 100h fuel class [t/ha]}
    \item{\code{Load_Live_Herb}}{Loading of herbaceous fuels [t/ha]}
    \item{\code{Load_Live_Woody}}{Loading of woody fuels [t/ha]}
    \item{\samp{SA/V_1h}}{Surface area to volume ratio of 1h fuel class [m2/m3]}
    \item{\samp{SA/V_10h}}{Surface area to volume ratio of 10h fuel class [m2/m3]}
    \item{\samp{SA/V_100h}}{Surface area to volume ratio of 100h fuel class [m2/m3]}
    \item{\samp{SA/V_Live_Herb}}{Surface area to volume ratio of herbaceous fuels [m2/m3]}
    \item{\samp{SA/V_Live_Woody}}{Surface area to volume ratio of woody fuels [m2/m3]}
    \item{\code{Fuel_Bed_Depth}}{Fuel bed depth [cm]}
    \item{\code{Mx_dead}}{Dead fuel moisture of extinction [percent]}
    \item{\code{Heat_1h}}{Heat content of 1h fuel class [kJ/kg]}
    \item{\code{Heat_10h}}{Heat content of 10h fuel class [kJ/kg]}
    \item{\code{Heat_100h}}{Heat content of 100h fuel class [kJ/kg]}
    \item{\code{Heat_Live_Herb}}{Heat content of herbaceous fuels [kJ/kg]}
    \item{\code{Heat_Live_Woody}}{Heat content of woody fuels [kJ/kg]}
  }
}
\source{
Albini, F. A. (1976). Computer-based models of wildland fire behavior: A users' manual. Ogden, UT: US Department of Agriculture, Forest Service, Intermountain Forest and Range Experiment Station.

Scott, J., & Burgan, R. E. (2005). A new set of standard fire behavior fuel models for use with Rothermel's surface fire spread model. Gen. Tech. Rep. RMRSGTR-153. Fort Collins, CO: US Department of Agriculture, Forest Service, Rocky Mountain Research Station.
}
\seealso{
 \code{\link{fire.Rothermel}}
}
\examples{
data(SFM_metric)
}
\keyword{datasets}
