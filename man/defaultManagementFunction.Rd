\encoding{UTF-8}
\name{defaultManagementFunction}
\alias{defaultManagementFunction}
\alias{defaultManagementArguments}
\title{
Default forest management actions
}
\description{
Function \code{defaultManagementFunction} implements management actions for 'regular' and 'irregular' monospecific stands, whereas function \code{defaultManagementArguments} returns a list with default values for the parameters regulating management.
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
TO BE WRITTEN
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
Function \code{defaultManagementArguments} returns a list with default arguments.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF

Jose-Ramon Gonzalez Olabarria, CTFC
}
\seealso{
\code{\link{fordyn}}
}
