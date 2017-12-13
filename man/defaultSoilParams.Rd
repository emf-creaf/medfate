\encoding{UTF-8}
\name{defaultSoilParams}
\alias{defaultSoilParams}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Default soil parameters}
\description{
Creates a list with default soil parameters for model functions.
}
\usage{
defaultSoilParams(n = 4) 
}
\arguments{
 \item{n}{An integer with the number of soil layers.}
}
\details{
The function returns a list with default parameters. Users can change those that need to be set to other values and use the list as input for function \code{\link{soil}}.
}
\value{
A list with the following elements (and default values):
\itemize{
\item{\code{widths (= c(300,700,1000,2000)}: Width of soil layers (in mm).}
\item{\code{clay (= 25)}: Clay percentage for each layer.}
\item{\code{sand (= 25)}: Sand percentage for each layer.}
\item{\code{om (= NA)}: Organic matter percentage for each layer.}
\item{\code{macro (= 0.1)}: Macroporosity for each layer.}
\item{\code{rfc (= c(20,40,60,85))}: Percentage of rock fragment content for each layer.}
\item{\code{Gsoil (= 0.5)}: Gamma parameter for infiltration.}
\item{\code{Ksoil (= 0.05)}: Kappa parameter for infiltration.}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
 \code{\link{soil}}, \code{\link{soilgridsParams}}, \code{\link{defaultControl}}, \code{\link{SpParamsMED}}
}
\examples{
defaultSoilParams(4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.