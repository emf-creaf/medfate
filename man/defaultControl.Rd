\encoding{UTF-8}
\name{defaultControl}
\alias{defaultControl}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Default control parameters for models}
\description{
Creates a list with global default parameters for simulation models.
}
\usage{
defaultControl() 
}
\details{
The function returns a list with default parameters. Users can change those defaults that need to be set to other values and use the list as input for model functions. The relevant parameters are different for each model function.
}
\value{
A list, with the following options:
\itemize{
\item{\code{verbose (=TRUE)}: Boolean flag to indicate console output during calculations.}
\item{\code{transpirationMode (="Simple")}: Transpiration model. See \code{\link{swbInput}}.}
\item{\code{cavitationRefill (= TRUE)}: Whether refilling of embolized conduits is activated.}
\item{\code{fracTotalTreeResistance (=0.625)}: Fraction of the minimum tree (stem+root) resistance that corresponds to stem resistance.}
\item{\code{averageFracRhizosphereResistance (=0.05)}: Fraction to total continuum (stem+root+rhizosphere) resistance that corresponds to rhizosphere (averaged across soil water potential values).}
\item{\code{storagePool (="none")}: Whether carbon storage pools are considered (either "none", "one" or "two").}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
 \code{\link{swbInput}}, \code{\link{swb}}, \code{\link{swbpoints}}
}