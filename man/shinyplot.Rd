\name{shinyplot}
\alias{shinyplot}
\title{
Shiny app with interactive plots
}
\description{
Creates a shiny app with interactive plots for simulation results and evaluation
}
\usage{
shinyplot(out, measuredData = NULL, SpParams = NULL)
}
\arguments{
  \item{out}{An object of class \code{\link{spwb}} or \code{\link{growth}}}
  \item{measuredData}{A data frame with observed/measured values (see \code{\link{evaluation_plot}}).}
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}}).}
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
if (interactive()) {
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforestMED)

#Default species parameterization
data(SpParamsMED)

#Initialize soil with default soil params (4 layers)
examplesoil1 = soil(defaultSoilParams(4))

#Initialize control parameters
control = defaultControl("Granier")

#Initialize input
x1 = forest2spwbInput(exampleforestMED,examplesoil1, SpParamsMED, control)

#Call simulation function
S1<-spwb(x1, examplemeteo, latitude = 41.82592, elevation = 100)

#Load observed data (in this case the same simulation results with some added error)  
data(exampleobs)

#Call interactive plot shiny app
shinyplot(S1, exampleobs, SpParamsMED)
}
}
