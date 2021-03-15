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
spwb_waterUseEfficiency(x, type = "Plant Ag/E", leaves = "average", freq = "days", 
                        draw = TRUE, ylim = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{spwb}, \code{pwb} or \code{growth}.}
  \item{type}{A string to indicate the scale of WUE calculation. Either:
    \itemize{
      \item{\code{"Leaf iWUE"}: Leaf intrinsic WUE, i.e. instantaneous ratio between photosynthesis and stomatal conductance (only for simulations with \code{transpirationMode = "Sperry"} and \code{subdailyResults = TRUE}). }
      \item{\code{"Leaf Ci"}: Leaf intercellular CO2 concentration (only for simulations with \code{transpirationMode = "Sperry"} and \code{subdailyResults = TRUE}).}
      \item{\code{"Plant An/E"}: Plant (cohort) net photosynthesis over plant transpiration (only for simulations with \code{transpirationMode = "Sperry"})}
      \item{\code{"Stand An/E"}: Stand net photosynthesis over stand transpiration (only for simulations with \code{transpirationMode = "Sperry"})}
      \item{\code{"Plant Ag/E"}: Plant (cohort) gross photosynthesis over plant transpiration}
      \item{\code{"Stand Ag/E"}: Stand gross photosynthesis over stand transpiration}
    }
  }
  \item{leaves}{Either \code{"sunlit"}, \code{"shade"} or \code{"average"}. Refers to the WUE of different leaf types or the average (with weights according to the LAI of sunlit and shade leaves). Only relevant for \code{type = "iWUE"}. }
  \item{freq}{Frequency of summary statistics (see \code{\link{cut.Date}}).}
  \item{draw}{A boolean flag to indicate that a plot should be returned.}
  \item{ylim}{Range of values for y.}
}
\details{
Temporal aggregation of WUE values is done differently depending on the value of \code{type}. For \code{type = "Plant Ag/E"}, \code{type = "Stand Ag/E"}, \code{type = "Plant An/E"} and \code{type = "Stand An/E"} sums or daily photosynthesis and transpiration are first calculated at the desired temporal scale and the ratio is calculated afterwards. For \code{type = "Leaf iWUE"} intrinsic WUE values are first calculated at the daily scale (as averages of instantaneous An/gs ratios weighted by An) and then they are aggregated to the desired scale by calculating weighted averages, where weights are given by daily photosynthesis.
}
\value{
If \code{draw=TRUE} a plot is returned. Otherwise, the function returns a matrix with WUE values, where rows are dates (at the desired temporal scale), and columns are plant cohorts. In the case of \code{type = "Plant Ag/E"}, \code{type = "Stand Ag/E"}, \code{type = "Plant An/E"} and \code{type = "Stand An/E"} values are in gC/L. In the case of \code{type = "Leaf iWUE"} values are in micromol of carbon per mmol of water.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{
\code{\link{spwb}}
}
