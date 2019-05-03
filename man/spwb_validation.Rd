\encoding{UTF-8}
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
  \item{x}{An object of class \code{\link{spwb}} or \code{\link{pwb}}.}
  \item{measuredData}{A data frame with observed/measured values. Dates should be in row names, whereas columns should be named according to the type of output to be evaluated (see details).}
  \item{type}{A string with the kind of model output to be evaluated. Accepted values are \code{"SWC"} (soil moisture content), \code{"ETR"} (total evapotranspiration), \code{"E"} (plant transpiration) and \code{"WP"} (plant water potentials).}
  \item{cohort}{A string of the cohort to be compared (e.g. "T1_68"). If \code{NULL} all possible cohorts are compared.}
  \item{draw}{A boolean flag to indicate that plots should be produced.}
}
\details{
   Users should provide the appropriate columns in \code{measuredData}, depending on the type of output to be evaluated:
   \itemize{
     \item{\code{"SWC"}: A column named \code{"SWC"} should be present, containing soil moisture content in percent volume.}
     \item{\code{"ETR"}: A column named \code{"ETR"} should be present, containing stand's evapotranspiration in mm/day. Observed values will be compared against modelled evapotranspiration (i.e. sum of transpiration, soil evaporation and interception loss) as well as against the sum of transpiration and soil evaporation only.}
     \item{\code{"E"}: For each plant cohort whose transpiration is to be evaluated, a column starting with \code{"E_"} and continuing with a cohort name (e.g. \code{"E_T1_68"}) with transpiration in mm/day.}
     \item{\code{"WP"}: For each plant cohort whose transpiration is to be evaluated, two columns, one starting with \code{"PD_"} (for pre-dawn) and the other with \code{"MD_"} (for midday), and continuing with a cohort name (e.g. \code{"PD_T1_68"}). They should contain leaf water potential values in MPa. These are compared against sunlit water potentials.}
   }
   Additional columns may exist with the standard error of measured quantities. These should be named as the referred quantity, followed by \code{"_err"} (e.g. \code{"PD_T1_68_err"}), and are used to draw confidence intervals around observations.
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
