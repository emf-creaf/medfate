\encoding{UTF-8}
\name{modifySpParams}
\alias{modifySpParams}
\alias{modifyCohortParams}
\title{
Modify species/cohort parameters
}
\description{
Routines to modify input species parameters or model input objects
}
\usage{
modifySpParams(SpParams, customParams, subsetSpecies = TRUE)
modifyCohortParams(x, customParams)
}
\arguments{
  \item{x}{A model input object of class \code{\link{spwbInput}}.}
  \item{SpParams}{A species parameter data frame, typically \code{\link{SpParamsMED}}.}
  \item{customParams}{A data frame with new parameter values (see details).}
  \item{subsetSpecies}{A flag to indicate that the only species mentioned in \code{customParams} should be returned.}
}
\details{
Data frame \code{customParams} should have as many rows as species or cohorts and as many columns as parameters to modify. An additional column is needed to identify species/cohorts. In the case of \code{modifySpParams} the column should be named 'SpIndex' (to match the corresponding column of\code{\link{SpParamsMED}}), while in \code{modifyCohortParams} the column should be named 'Cohort' and will be matched with the cohort names given by \code{\link{spwbInput}} or \code{\link{growthInput}}.
}
\value{
Function \code{modifySpParams} returns a modified species parameter data frame. Function \code{modifyCohortParams} returns a modified \code{\link{spwbInput}} object.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{
\code{\link{spwbInput}}, \code{\link{SpParamsMED}}
}
\examples{
## TO DO
}
