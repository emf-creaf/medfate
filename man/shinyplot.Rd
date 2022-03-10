\name{shinyplot}
\alias{shinyplot}
\alias{shinyplot.spwb}
\alias{shinyplot.pwb}
\alias{shinyplot.growth}
\alias{shinyplot.fordyn}
\alias{shinyplot.spwb_day}
\alias{shinyplot.pwb_day}
\alias{shinyplot.growth_day}
\title{
Shiny app with interactive plots
}
\description{
Creates a shiny app with interactive plots for simulation results and evaluation
}
\usage{
\method{shinyplot}{spwb}(x, measuredData = NULL, ...)
\method{shinyplot}{pwb}(x, measuredData = NULL, ...)
\method{shinyplot}{growth}(x, measuredData = NULL, ...)
\method{shinyplot}{fordyn}(x, measuredData = NULL, ...)
\method{shinyplot}{spwb_day}(x, ...)
\method{shinyplot}{pwb_day}(x, ...)
\method{shinyplot}{growth_day}(x, ...)
}
\arguments{
  \item{x}{An object of the right class.}
  \item{measuredData}{A data frame with observed/measured values (see \code{\link{evaluation_plot}}).}
  \item{...}{Additional parameters.}
}
\details{Only run this function in interactive mode. When \code{measuredData} is not \code{NULL}, an additional panel is shown for evaluation plots.}
\value{An object that represents the shiny app}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{
\code{\link{plot.spwb}}, \code{\link{evaluation_plot}}
}
\examples{
## Only run this example in interactive R sessions
\dontrun{
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforestMED)

#Default species parameterization
data(SpParamsMED)

#Initialize soil with default soil params (4 layers)
examplesoil = soil(defaultSoilParams(4))

#Initialize control parameters
control = defaultControl("Granier")

#Initialize input
x1 = forest2spwbInput(exampleforestMED, examplesoil, SpParamsMED, control)

#Call simulation function
S1<-spwb(x1, examplemeteo, latitude = 41.82592, elevation = 100)

#Load observed data (in this case the same simulation results with some added error)  
data(exampleobs)

#Call interactive plot shiny app
shinyplot(S1, exampleobs)
}
}
