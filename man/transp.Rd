\encoding{UTF-8}
\name{transp}
\alias{transp.profitMaximization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transpiration submodel functions
}
\description{
Set of high-level functions used in the calculation of stomatal conductance and transpiration.
}
\usage{
transp.profitMaximization(supplyFunction, photosynthesisFunction, type = 1, kstemmax = NA)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{supplyFunction}{Water supply function (see \code{\link{hydraulics.supplyFunctionNetwork}}).}
  \item{photosynthesisFunction}{Function returned by \code{photo.photosynthesisFunction()}.}
  \item{type}{Optimization type (1,2 or 3).}
  \item{kstemmax}{Maximum stem conductance (flow per pressure drop).}
}
\details{
Details of the transpiration submodel are given in a vignette.
}
\value{
Values returned for each function are:
\itemize{
\item{\code{transp.profitMaximization}:Returns a list with the following elements:
\itemize{
  \item{\code{Cost}: Cost function [0-1].}
  \item{\code{Gain}: Gain function [0-1].}
  \item{\code{Profit}: Profit function [0-1].}
  \item{\code{iMaxProfit}: Index corresponding to maximum profit (starting from 0).}
}
}
}
}
\references{
Sperry, J. S., M. D. Venturas, W. R. L. Anderegg, M. Mencuccini, D. S. Mackay, Y. Wang, and D. M. Love. 2016. Predicting stomatal responses to the environment from the optimization of photosynthetic gain and hydraulic cost. Plant Cell and Environment.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
\code{\link{hydraulics.supplyFunctionNetwork}}, \code{\link{biophysics.leafTemperature}}, \code{\link{photo.photosynthesis}}, \code{\link{swb}}
}