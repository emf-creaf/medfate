\encoding{UTF-8}
\name{spwb_ldrCalibration}
\alias{spwb_ldrCalibration}

\title{Calibration of root distribution}

\description{The function \code{spwb_ldrCalibration} calibrates the species root
distribution within \code{\link{spwb}}, given the arguments \code{x}, \code{meteo},
\code{psi_crit}, \code{obs} and \code{calibVar}. This calibration
is based on reference measured values. These reference measured values can be
\code{Soil water content}, \code{Total tranpiration} or \code{Transpiration by
cohort}. Return the calibrated root distribution for each tree species (no
shrub calibration is done), expressed as parameters of the function
\code{\link{root_ldrDistribution}}.}

\usage{
spwb_ldrCalibration(x, meteo, calibVar, obs,
                   RZmin = 301, RZmax = 4000,
                   V1min = 0.01, V1max = 0.94, resolution = 20, heat_stop = 0, 
                   transformation = "identity", verbose = FALSE)
}

\arguments{
  \item{x}{An object of class \code{\link{spwbInput}}.}
  \item{meteo}{A data frame with daily meteorological data series. When using the 'Granier' transpiration mode the following columns are required:
   \itemize{
     \item{\code{DOY}: Day of the year (Julian day).}
     \item{\code{Precipitation}: Precipitation (in mm).}
     \item{\code{MeanTemperature}: Mean temperature (in degrees Celsius).}
     \item{\code{PET}: Potential evapotranspiration (in mm).}
    }
  }
  \item{calibVar}{A character string indicating the calibration variable to be used.
                  It can be one of the following: \code{SWC}, \code{Eplanttot}
                  or \code{Cohorts}.}
  \item{obs}{Measured calibration variable. Depending on the value of \code{calibVar}
             it can be a numeric vector with the measured SWC values (if
             \code{calibVar = "SWC"}), or a data frame with the first column
             containing the measured total transpiration (named \code{Eplanttot}) and
             the following columns containing the cohorts transpiration.}
  \item{RZmin}{The minimum value of RZ (the rooting depth) to be explored (in mm)}
  \item{RZmax}{The maximum value of RZ (the rooting depth) to be explored (in mm)}
  \item{V1min}{The minimum value of V1 (the root proportion in the first soil layer) to be explored}
  \item{V1max}{The maximum value of V1 (the root proportion in the first soil layer) to be explored}
  \item{resolution}{An integer defining the number of values to obtain by discretization of the root parameters RZ and V1. The number of parameter combinations and therefore the computation cost increases increase with the square of resolution}
  \item{transformation}{Function to modify the size of Z intervals to be explored (by default, bins are equal).}
  \item{heat_stop}{An integer defining the number of days during to discard from the calculation of the optimal root distribution. Usefull if the soil water content initialization is not certain}
  \item{verbose}{A logical value. Print the internal messages of the function?}
}

\details{
This function performs three different kinds of calibration, selecting those
root distribution parameters that minimize the MAE between the predicted values
and the measured values provided in \code{obs} argument. If \code{calibVar = "SWC"} different
V1 values are tested running \code{\link{spwb}} maintaining the total soil depth
provided in \code{x} and assuming that value is also the depth containing 95 percent of the
roots. If \code{calibVar = "Eplanttot"} or \code{calibVar = 'Cohorts'} different 
combinations of RZ and V1 values are tested for each tree cohort and the root
paramters are selected based on the MAE between the total transpiration or the
cohort transpiration.
}

\value{
THe function returns a data frame containing the species index used in medfate,
calibrated values for Z50, Z95 and V1 and the MAE value for that combination.
}
\author{
\enc{Víctor}{Victor} Granda, CREAF

Antoine Cabon, CTFC-CREAF

Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{
 \code{\link{spwb_ldrOptimization}} for when no measured data is available,
 \code{\link{spwb}}, \code{\link{soil}}, \code{\link{root_ldrDistribution}}
}
