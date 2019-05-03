\name{spwb_validation}
\alias{spwb_validation}
\title{
Validation of soil-plant water balance simulations
}
\description{
Compares model predictions against observed values
}
\usage{
spwb_validation(x, measuredData, type = "SWC", cohort = NULL, draw = TRUE)
}
\arguments{
  \item{x}{An object of class \code{\link{spwb}}.}
  \item{measuredData}{A data frame with observed/measured values (see details).}
  \item{type}{A string with the kind of model output to be evaluated:
    \itemize{
      \item{\code{"SWC"}: Soil moisture content.}
      \item{\code{"ETR"}: Total evapotranspiration.}
      \item{\code{"E"}: Plant transpiration.}
    }
  }
  \item{cohort}{A string of the cohort to be compared (e.g. "T1_68").}
  \item{draw}{A boolean flag to indicate that plots should be produced.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Evaluation statistics (per cohort, when multiple cohorts evaluated):
\itemize{
  \item{\code{MAE}: Mean absolute error.}
  \item{\code{Bias}: Mean deviation (positive values correspond to model overestimations).}
  \item{\code{R2}: Squared linear correlation coefficient.}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{
\code{\link{spwb}}
}
