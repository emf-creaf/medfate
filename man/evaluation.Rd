\encoding{UTF-8}
\name{evaluation}
\alias{evaluation}
\alias{evaluation_table}
\alias{evaluation_stats}
\alias{evaluation_plot}
\alias{evaluation_metric}
\title{
Evaluation of simulations results
}
\description{
Functions to compare model predictions against observed values.
}
\usage{
evaluation_table(out, measuredData, type = "SWC", cohort = NULL, 
                 temporalResolution = "day", SpParams = NULL)
evaluation_stats(out, measuredData, type="SWC", cohort = NULL, 
                 temporalResolution = "day", SpParams = NULL)
evaluation_plot(out, measuredData, type="SWC", cohort = NULL, 
                temporalResolution = "day", SpParams = NULL, 
                plotType = "dynamics")
evaluation_metric(out, measuredData, type="SWC", cohort=NULL, 
                  temporalResolution = "day", SpParams = NULL,
                  metric = "loglikelihood")
}
\arguments{
  \item{out}{An object of class \code{\link{spwb}}, \code{\link{growth}} or \code{\link{pwb}}.}
  \item{measuredData}{A data frame with observed/measured values. Dates should be in row names, whereas columns should be named according to the type of output to be evaluated (see details).}
  \item{type}{A string with the kind of model output to be evaluated. Accepted values are \code{"SWC"} (soil moisture content), \code{"REW"} relative extractable water, \code{"ETR"} (total evapotranspiration), \code{"E"} (transpiration per leaf area), \code{"FMC"} (fuel moisture content) and \code{"WP"} (plant water potentials).}
  \item{cohort}{A string of the cohort to be compared (e.g. "T1_68"). If \code{NULL} results for the first cohort will be evaluated.}
  \item{temporalResolution}{A string to indicate the temporal resolution of the model evaluation, which can be "day", "week", "month" or "year". Observed and modelled values are aggregated temporally (using either means or sums) before comparison.}
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}}), only needed if \code{type = "FMC"}.}
  \item{plotType}{Plot type to draw, either \code{"dynamics"} or \code{"scatter"}.}
  \item{metric}{An evaluation metric:
    \itemize{
      \item{\code{"MAE"}: Mean absolute error.}
      \item{\code{"r"}: Pearson's linear correlation coefficient.}
      \item{\code{"NSE"}: Nash-Sutcliffe model efficiency coefficient.}
      \item{\code{"NSEabs"}: Modified Nash-Sutcliffe model efficiency coefficient (L1 norm) (Legates & McCabe 1999).}
      \item{\code{"loglikelihood"}: Logarithm of the likelihood of observing the data given the model predictions, assuming independent Gaussian errors.}
    }
  }
}
\details{
   Users should provide the appropriate columns in \code{measuredData}, depending on the type of output to be evaluated:
   \itemize{
     \item{\code{"SWC" or "REW"}: A column named \code{"SWC"} should be present, containing soil moisture content in percent volume.}
     \item{\code{"ETR"}: A column named \code{"ETR"} should be present, containing stand's evapotranspiration in mm/day. Observed values will be compared against modelled evapotranspiration (i.e. sum of transpiration, soil evaporation and interception loss) as well as against the sum of transpiration and soil evaporation only.}
     \item{\code{"E"}: For each plant cohort whose transpiration is to be evaluated, a column starting with \code{"E_"} and continuing with a cohort name (e.g. \code{"E_T1_68"}) with transpiration in mm/day on a leaf area basis.}
     \item{\code{"FMC"}: For each plant cohort whose transpiration is to be evaluated, a column starting with \code{"FCM_"} and continuing with a cohort name (e.g. \code{"FMC_T1_68"}) with fuel moisture content as percent of dry weight.}
     \item{\code{"WP"}: For each plant cohort whose transpiration is to be evaluated, two columns, one starting with \code{"PD_"} (for pre-dawn) and the other with \code{"MD_"} (for midday), and continuing with a cohort name (e.g. \code{"PD_T1_68"}). They should contain leaf water potential values in MPa. These are compared against sunlit water potentials.}
   }
   Additional columns may exist with the standard error of measured quantities. These should be named as the referred quantity, followed by \code{"_err"} (e.g. \code{"PD_T1_68_err"}), and are used to draw confidence intervals around observations.
}
\value{
\itemize{
  \item{Function \code{evaluation_table} returns a data frame with observed and predicted values.}
  \item{
    Function \code{evaluation_stats} returns evaluation statistics (a vector or a data frame depending on \code{type}):
    \itemize{
      \item{\code{MAE}: Mean absolute error.}
      \item{\code{Bias}: Mean deviation (positive values correspond to model overestimations).}
      \item{\code{r}: Pearson's linear correlation coefficient.}
      \item{\code{NSE}: Nash-Sutcliffe model efficiency coefficient.}
      \item{\code{NSEabs}: Modified Nash-Sutcliffe model efficiency coefficient (L1 norm) (Legates & McCabe 1999).}
    }
  }
  \item{Function \code{evaluation_plot} returns a ggplot object.}
  \item{Function \code{evaluation_metric} returns a scalar with the desired metric.}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\references{
Legates, D.R., McCabe, G.J., 1999. Evaluating the use of “goodness-of-fit” measures in hydrologic and hydroclimatic model validation. Water Resour. Res. 35, 233–241. 
}
\seealso{
\code{\link{spwb}}, \code{\link{growth}}, \code{\link{optimization}}, \code{\link{exampleobs}}
}
\examples{
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforestMED)

#Default species parameterization
data(SpParamsMED)

#Initialize soil with default soil params (4 layers)
examplesoil1 = soil(defaultSoilParams(4))

#Initialize control parameters
control = defaultControl()

#Initialize input
x1 = forest2spwbInput(exampleforestMED,examplesoil1, SpParamsMED, control)

#Call simulation function
S1<-spwb(x1, examplesoil1, examplemeteo, latitude = 41.82592, elevation = 100)

#Load observed data (in this case the same simulation results with some added error)  
data(exampleobs)

#Evaluation statistics for soil water content
evaluation_stats(S1, exampleobs)

#NSE only
evaluation_metric(S1, exampleobs, metric="NSE")


#Comparison of temporal dynamics
evaluation_plot(S1, exampleobs)

#Loglikelihood value
evaluation_metric(S1, exampleobs)
}
