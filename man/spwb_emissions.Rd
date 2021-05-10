\name{spwb_emissions}
\alias{spwb_emissions}
\title{
Estimates BVOC emissions
}
\description{
Allows estimating biogenic volatile organic compound (BVOC) emissions from \code{\link{spwb}} objects using the MEGAN (2.1) model.
}
\usage{
spwb_emissions(x, verbose = FALSE)
}
\arguments{
  \item{x}{An object of class \code{spwb}.}
  \item{verbose}{A boolean flag to indicate extra console output.}
}
\details{
TO BE DONE
}
\value{
A matrix of BVOC emissions by day and compound class.
}
\references{
Guenther et al. 2012
}
\author{
Shengli Huang
}

\seealso{
\code{\link{spwb}}
}

