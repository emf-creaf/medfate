\encoding{UTF-8}
\name{defaultSoilParams}
\alias{defaultSoilParams}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Default soil parameters}
\description{
Creates a data frame with default soil physical description for model functions.
}
\usage{
defaultSoilParams(n = 4) 
}
\arguments{
 \item{n}{An integer with the number of soil layers.}
}
\details{
The function returns a data frame with default physical soil description, with soil layers in rows. Users can change those that need to be set to other values and use the list as input for function \code{\link{soil}}.
}
\value{
A data frame with layers in rows and the following columns (and default values):
\itemize{
\item{\code{widths (= c(300,700,1000,2000)}: Width of soil layers (in mm).}
\item{\code{clay (= 25)}: Clay percentage for each layer (in \%).}
\item{\code{sand (= 25)}: Sand percentage for each layer (in \%).}
\item{\code{om (= NA)}: Organic matter percentage for each layer (in \%).}
\item{\code{bd (= 1.5)}: Bulk density for each layer (in g/cm3).}
\item{\code{rfc (= c(20,40,60,85))}: Percentage of rock fragment content for each layer.}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{
 \code{\link{soil}}, \code{\link{soilgridsParams}}, \code{\link{defaultControl}}, \code{\link{SpParamsMED}}
}
\examples{
defaultSoilParams(4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.