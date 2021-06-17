\encoding{UTF-8}
\name{forest}
\alias{forest}
\alias{emptyforest}
\alias{summary.forest}
\alias{print.summary.forest}

\title{Forest description}
\description{
Description of a forest stand.
}
\usage{
\method{summary}{forest}(object, SpParams, mode = "MED", detailed=FALSE, ...) 
\method{print}{summary.forest}(x, digits = getOption("digits"), ...)
emptyforest(ID="", patchsize=10000, ntree = 0, nshrub = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{forest} has the following structure:
  \itemize{
  \item{\code{ID}: An identifier of the forest stand (a string).}
  \item{\code{patchsize}: The area that forest stand description represents, in square meters.}
  \item{\code{treeData}: A data frame of tree cohorts (in rows) and the following columns:
   \itemize{
      \item{\code{Species}: Non-negative integer for tree species identity (i.e., 0,1,2,...).}
      \item{\code{Height}: Total height (in cm).}
      \item{\code{DBH}: Diameter at breast height (in cm).}
      \item{\code{N}: Density (number of individuals/cell).}
      \item{\code{Z50}: Depth (in mm) corresponding to 50\% of fine roots.}
      \item{\code{Z95}: Depth (in mm) corresponding to 95\% of fine roots.}
   }
  }
  \item{\code{shrubData}: A data frame of shrub cohorts (in rows) and the following columns:
   \itemize{
      \item{\code{Species}: Non-negative integer for shrub species identity (i.e., 0,1,2,...).}
      \item{\code{Height}: Total height (in cm).}
      \item{\code{Cover}: Percent cover.}
      \item{\code{Z50}: Depth (in mm) corresponding to 50\% of fine roots.}
      \item{\code{Z95}: Depth (in mm) corresponding to 95\% of fine roots.}
   }
  }
  \item{\code{seedBank}: A data frame containing the abundance of seeds for each species (in rows) and the following columns:
   \itemize{
      \item{\code{Species}: Non-negative integer for shrub species identity (i.e., 0,1,2,...).}
      \item{\code{Abundance}: Abundance class (0 - none; 1 - low; 2 - medium; 3 - high; 4 - very high).}
   }
  }
  \item{\code{herbCover}: Percent cover of the herb layer.}
  \item{\code{herbHeight}: Mean height (in cm) of the herb layer.}
  }
  }
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}}).}
  \item{mode}{Calculation mode, either "MED" or "US".}
  \item{detailed}{A boolean flag to indicate that a detailed summary is desired.}
  \item{x}{The object returned by \code{summary.forest}.}
  \item{digits}{Minimal number of significant digits.}
  \item{...}{Additional parameters for functions \code{\link{summary}} and \code{\link{print}}.}
  \item{ID}{An identifier of the forest stand (a string).}
  \item{patchsize}{The area of the forest stand, in square meters.}
  \item{ntree, nshrub}{Number of tree and shrub cohorts, respectively.}
}
\details{
  Function \code{summary.forest} can be used to summarize a \code{forest} object in the console. Function \code{emptyforest} creates an empty \code{forest} object.
}
\value{
  Function \code{summary.forest} returns a list with the basal area and LAI of the forest, either expressed as totals or divided among life stages and species. Function \code{emptyforest} returns an empty \code{forest} object.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{
 \code{\link{exampleforestMED}}, \code{\link{forest_mergeTrees}}
}
\examples{
data(exampleforestMED)
data(SpParamsMED)

summary(exampleforestMED, SpParamsMED)
}
