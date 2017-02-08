\encoding{UTF-8}
\name{soil}
\alias{soil}
\alias{print.soil}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Soil definition}
\description{
Defines a soil object for its use in simulation functions.
}
\usage{
soil(SoilParams, W = as.numeric(c(1)))
\method{print}{soil}(x,  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SoilParams}{A list of soil parameters (see \code{\link{defaultSoilParams}}.}
  \item{W}{A numerical vector with the relative water content of each soil layer.}
  \item{x}{An object of class \code{soil}.}
  \item{...}{Additional parameters to \code{print}.}
}
\details{
Function \code{print} prompts a description of soil characteristics and status (water content). Volume at field capacity is calculated from assuming a soil water potential equal to -0.033 MPa.
}
\value{
An list of class \code{soil} with the following elements:
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
  \code{\link{soil.psi2theta}}, \code{\link{swb}}, \code{\link{defaultSoilParams}}
}
\examples{
s = soil(defaultSoilParams())
s
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.