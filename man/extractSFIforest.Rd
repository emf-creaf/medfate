\encoding{UTF-8}
\name{extractSFIforest}
\alias{extractSFIforest}
\alias{translateSpeciesCodes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract forest from SFI data
}
\description{
Creates a \code{\link{forest}} object from Spanish Forest Inventory (SFI) data (DGCN 2005).
}
\usage{
extractSFIforest(SFItreeData, SFIshrubData, ID, SpParams,  
                 SFIherbData = NULL, SFIcodes=NULL, 
                 patchsize= 10000, setDefaults=TRUE)
translateSpeciesCodes(x, SFIcodes)                 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SFItreeData}{A data frame with measured tree data.}
  \item{SFIshrubData}{A data frame with measured shrub data.}
  \item{ID}{A string with the ID of the plot to be extracted.}
  \item{SpParams}{A data frame with species parameters (see details).}
  \item{SFIcodes}{A string vector (of length equal to the number of rows in \code{SpParams} of the SFI species codes that correspond to the model species codification. Each string may contain different coma-separated codes in order to merge SFI species into a single model species.}
  \item{SFIherbData}{A data frame with cover and mean height of the herb layer.}
  \item{patchsize}{The area of the forest stand, in square meters.}
  \item{setDefaults}{Initializes default values for missing fields in SFI data.}
  \item{x}{A data frame with a column called 'Especie'.}  
}
\details{
SFI data needs to be in a specific format. Function \code{extractSFIforest} calls \code{translateSpeciesCodes} internally.
}
\value{
An object of class \code{\link{forest}}.
}
\references{
DGCN (2005). Tercer Inventario Forestal Nacional (1997-2007): Catalunya. Dirección General de Conservación de la Naturaleza, Ministerio de Medio Ambiente, Madrid.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{\code{\link{forest}}}
