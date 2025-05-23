% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/waterUseEfficiency.R
\name{waterUseEfficiency}
\alias{waterUseEfficiency}
\title{Water use efficiency}
\usage{
waterUseEfficiency(
  x,
  type = "Plant Ag/E",
  leaves = "average",
  freq = "days",
  draw = TRUE,
  ylim = NULL
)
}
\arguments{
\item{x}{An object of class \code{\link{spwb}}, \code{\link{pwb}}, \code{\link{growth}} or \code{\link{fordyn}}.}

\item{type}{A string to indicate the scale of WUE calculation. Either:
\itemize{
\item{\code{"Leaf iWUE"}: Leaf intrinsic WUE, i.e. instantaneous ratio between photosynthesis and stomatal conductance (only for simulations with \code{transpirationMode = "Sperry"} or \code{transpirationMode = "Sureau"} and \code{subdailyResults = TRUE}). }
\item{\code{"Leaf Ci"}: Leaf intercellular CO2 concentration (only for simulations with \code{transpirationMode = "Sperry"} or \code{transpirationMode = "Sureau"} and \code{subdailyResults = TRUE}).}
\item{\code{"Plant An/E"}: Plant (cohort) net photosynthesis over plant transpiration (only for simulations with \code{transpirationMode = "Sperry"} or \code{transpirationMode = "Sureau"})}
\item{\code{"Stand An/E"}: Stand net photosynthesis over stand transpiration (only for simulations with \code{transpirationMode = "Sperry"} or \code{transpirationMode = "Sureau"})}
\item{\code{"Plant Ag/E"}: Plant (cohort) gross photosynthesis over plant transpiration}
\item{\code{"Stand Ag/E"}: Stand gross photosynthesis over stand transpiration}
}}

\item{leaves}{Either \code{"sunlit"}, \code{"shade"} or \code{"average"}. Refers to the WUE of different leaf types or the average (with weights according to the LAI of sunlit and shade leaves). Only relevant for \code{type = "iWUE"}.}

\item{freq}{Frequency of summary statistics (see \code{\link{cut.Date}}).}

\item{draw}{A boolean flag to indicate that a plot should be returned.}

\item{ylim}{Range of values for y.}
}
\value{
If \code{draw=TRUE} a plot is returned. Otherwise, the function returns a matrix with WUE values,
where rows are dates (at the desired temporal scale), and columns are plant cohorts.
In the case of \code{type = "Plant Ag/E"}, \code{type = "Stand Ag/E"}, \code{type = "Plant An/E"} and \code{type = "Stand An/E"} values are in gC/L.
In the case of \code{type = "Leaf iWUE"} values are in micromol of carbon per mmol of water.
}
\description{
Calculates plant water use efficiency (WUE), at different temporal scales, from simulation results.
}
\details{
Temporal aggregation of WUE values is done differently depending on the value of \code{type}.
For \code{type = "Plant Ag/E"}, \code{type = "Stand Ag/E"}, \code{type = "Plant An/E"} and \code{type = "Stand An/E"} sums
or daily photosynthesis and transpiration are first calculated at the desired temporal scale and the ratio is calculated afterwards.
For \code{type = "Leaf iWUE"} intrinsic WUE values are first calculated at the daily scale (as averages of instantaneous An/gs ratios weighted by An)
and then they are aggregated to the desired scale by calculating weighted averages, where weights are given by daily photosynthesis.
}
\seealso{
\code{\link{droughtStress}}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
