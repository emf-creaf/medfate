\encoding{UTF-8}
\name{spwb_waterUseEfficiency}
\alias{spwb_waterUseEfficiency}
\title{
Water use efficiency
}
\description{
Calculates water use efficiency (WUE), at different temporal scales, from soil plant water balance calculations
}
\usage{
spwb_waterUseEfficiency(x, type = "An/E", leaves = "average", freq = "days")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{spwb} or \code{pwb}.}
  \item{type}{Either \code{"An/E"} (i.e. photosynthesis over transpiration) or \code{"iWUE"} (intrinsic WUE, i.e. instantaneous ratio between photosynthesis and stomatal conductance). \code{type = "iWUE"} is only valid for simulations with \code{transpirationMode = "Sperry"} and \code{subdailyResults = TRUE}.}
  \item{leaves}{Either \code{"sunlit"}, \code{"shade"} or \code{"average"}. Refers to the WUE of different leaf types or the average (with weights according to the LAI of sunlit and shade leaves). Only relevant for \code{type = "iWUE"}. }
  \item{freq}{Frequency of summary statistics (see \code{\link{cut.Date}}).}
}
\details{
Temporal aggregation of WUE values is done differently depending on the value of \code{type}. For \code{type = "An/E"} sums or daily photosynthesis and transpiration are first calculated at the desired temporal scale and the ratio is calculated afterwards. For \code{type = "iWUE"} intrinsic WUE values are first calculated at the daily scale (as averages of instantaneous An/gs ratios weighted by An) and then they are aggregated to the desired scale by calculating weighted averages, where weights are given by daily photosynthesis.
}
\value{
A matrix with WUE values, where rows are dates (at the desired temporal scale), and columns are plant cohorts. In the case of \code{type = "An/E"} values are in gC/L. In the case of \code{type = "iWUE"} values are in micromol of carbon per mmol of water.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{
\code{\link{spwb}}
}
