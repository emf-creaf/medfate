\encoding{UTF-8}
\name{defaultManagementFunction}
\alias{defaultManagementFunction}
\alias{defaultManagementArguments}
\title{
Default forest management actions
}
\description{
Function \code{defaultManagementFunction} implements management actions for 'regular' and 'irregular' monospecific stands, whereas function \code{defaultManagementArguments} returns a list with default values for the parameters regulating management. Both functions are meant to be used in simulations with \code{\link{fordyn}}.
}
\usage{
defaultManagementFunction(x, args, verbose = FALSE)
defaultManagementArguments()
}
\arguments{
  \item{x}{An object of class \code{\link{forest}}}
  \item{args}{A list of arguments regulating management actions, e.g. the list returned by \code{defaultManagementArguments}}
  \item{verbose}{A boolean flag enabling console printing.}
}
\details{
This function implements silvicultural actions following either 'regular' or 'irregular' management models. Irregular models are implemented by executing thinning operations only, whereas regular models include a set of final cuts. Thinning occurs anytime a stand-level metric (e.g. basal area) crosses a given threshold, and different kinds of thinning operations are allowed. Final cuts start whenever mean DBH exceeds a given threshold, and may include different cuts separated a number of years. Planting is only allowed under regular management models, and is applied after the last final cut. Understory clearings are assumed to occur anytime there is an intervention on trees, an only a residual shrub cover is left.
}
\value{
Function \code{defaultManagementFunction} returns a list with the following items:
  \itemize{
     \item{\code{"action"}: A string identifying the action performed (e.g. "thinning").}
     \item{\code{"N_tree_cut"}: A vector with the density of trees removed.}
     \item{\code{"Cover_shrub_cut"}: A vector with the cover of shrubs removed.} 
     \item{\code{"planted_forest"}: An object of class \code{\link{forest}} with the new plant cohorts resulting from tree/shrub planting.}
     \item{\code{"management_args"}: A list of management arguments to be used in the next call to the management function.}
   }
Function \code{defaultManagementArguments} returns a list with default arguments:
 \itemize{
     \item{\code{"type"}: Management model, either 'regular' or 'irregular'.}
     \item{\code{"thinning"}: Kind of thinning to be applied in irregular models or in regular models before the final cuts. Options are 'below', 'above', 'intermediate', 'below-intermediate', 'above-intermediate' or a string with the proportion of cuts to be applied to different diameter sizes.}
     \item{\code{"thinningMetric"}: The stand-level metric used to decide whether thinning is applied, either 'BA' (basal area), 'N' (density) or 'HB' (Hart-Becking index).} 
     \item{\code{"thinningThreshold"}: The threshold value of the stand-level metric causing the thinning decision.}
     \item{\code{"thinningPerc"}: Percentage of stand's basal area to be removed in thinning operations.}
     \item{\code{"finalMeanDBH"}: Mean DBH threshold to start final cuts.}
     \item{\code{"finalPerc"}: String with percentages of basal area to be removed in final cuts, separated by '-' (e.g. "40-60-100").}
     \item{\code{"finalPreviousStage"}: Integer state variable to store the stage of final cuts ('0' before starting final cuts).}
     \item{\code{"finalYearsBetweenCuts"}: Number of years separating final cuts.}
     \item{\code{"finalYearsToCut"}: State variable to count the years to be passed before new final cut is applied.}
     \item{\code{"plantingSpecies"}: Species code to be planted. If missing, planting does not occur and only natural regeneration is allowed.}
     \item{\code{"plantingDBH"}: Initial DBH (cm) of planted species.}
     \item{\code{"plantingHeight"}: Initial height (cm) of planted species.}
     \item{\code{"plantingDensity"}: Initial density (ind./ha) of the planted species.}
     \item{\code{"understoryMaxCover"}: Percentage of overall shrub cover to be left after any silvicultural intervention.}
   }
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF

Jose-Ramon Gonzalez Olabarria, CTFC
}
\seealso{
\code{\link{fordyn}}
}
