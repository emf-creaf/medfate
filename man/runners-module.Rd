\name{runners}
\alias{runners}
\title{
	Rcpp module: Running simulations
}
\description{
	This Rcpp module allows running simulations using C++ background objects. They are intended for non-standard users.
}
\details{
	The module contains the following items: 

	classes: \describe{
        \item{single_runner}{An S4 class to run spwb_day() or growth_day() simulations on a single site.}
        \item{multiple_runner}{An S4 class to run spwb_day() or growth_day() simulations on multiple sites, if required using parallelization.}
        \item{watershed_runner}{An S4 class to run one day of watershed simulations.}
		}
}
\examples{
show( runners )
}
\keyword{datasets}
