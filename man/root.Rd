\encoding{UTF-8}
\name{root}
\alias{root_ldrDistribution}
\alias{root_conicDistribution}
\alias{root_xylemConductanceProportions}
\alias{root_individualRootedGroundArea}
\alias{root_coarseRootLengths}
\alias{root_fineRootRadius}
\alias{root_fineRootArea}
\alias{root_fineRootBiomass}
\alias{root_fineRootSoilVolume}
\alias{root_specificRootSurfaceArea}
\alias{root_individualRootedGroundArea}
\alias{root_horizontalProportionsBasic}
\alias{root_horizontalProportionsAdvanced}
\title{Distribution of fine roots}
\description{
Functions to calculate the distribution of fine roots within the soil, given root system parameters and soil layer definition.
}
\usage{
root_conicDistribution(Zcone, d)
root_ldrDistribution(Z50, Z95, d)
root_xylemConductanceProportions(v, d, depthWidthRatio = 1)
root_coarseRootLengths(v, d, depthWidthRatio = 1.0)
root_fineRootRadius(specificRootLength, rootTissueDensity)
root_fineRootArea(vgrhizo_kmax, leafArea)
root_fineRootBiomass(vgrhizo_kmax, leafArea, 
                     specificRootLength, rootTissueDensity)
root_fineRootSoilVolume(fineRootBiomass, specificRootLength, rootLengthDensity = 10.0)
root_specificRootSurfaceArea(specificRootLength, rootTissueDensity)
root_individualRootedGroundArea(VolInd, V, d, bulkDensity)
root_horizontalProportionsBasic(poolProportions, V, LAIcell, poolOverlapFactor)
root_horizontalProportionsAdvanced(poolProportions, VolInd, N, V, 
                                   d, bulkDensity)
}
\arguments{
  \item{Z50}{A vector of depths (in mm) corresponding to 50\% of roots.}
  \item{Z95}{A vector of depths (in mm) corresponding to 95\% of roots.}
  \item{Zcone}{A vector of depths (in mm) corresponding to the root cone tip.}
  \item{d}{The width (in mm) corresponding to each soil layer.}
  \item{bulkDensity}{Bulk density (g/cm3) of each soil layer.}
  \item{v}{Vector with proportions of fine roots, as returned by functions \code{root_conicDistribution} or \code{root_ldrDistribution}.}
  \item{depthWidthRatio}{Ratio between radius of the soil layer with the largest radius and maximum rooting depth.}
  \item{specificRootLength}{Specific fine root length (length of fine roots over weight).}
  \item{rootTissueDensity}{Fine root tissue density (weight over volume at turgidity).}
  \item{vgrhizo_kmax}{Rhizosphere maximum conductance per leaf area (mmol·m-2·s-1·MPa-1).}
  \item{leafArea}{Leaf area (in m2).}
  \item{rootLengthDensity}{Fine root length density (length of fine roots over soil volume; cm/cm3)}
  \item{fineRootBiomass}{Biomass of fine roots (g).}
  \item{V}{Matrix of proportions of fine roots (cohorts x soil layers).}
  \item{LAIcell}{Leaf Area Index of the stand.}
  \item{poolOverlapFactor}{Factor to regulate the influence of LAI on root overlap.}
  \item{VolInd}{Volume of soil (in m3) occupied by fine roots per individual. }
  \item{N}{Density of individuals per hectare.}
  \item{poolProportions}{Division of the stand area among plant cohorts (proportions).}
}
\details{
  Function \code{root_conicDistribution} assumes a conic distribution of fine roots, whereas function \code{root_ldrDistribution} distributes fine roots according to the linear dose response model of Schenck & Jackson (2002). Function \code{root_xylemConductanceProportions} calculates the proportion of total root xylem conductance that can be attributed to each layer, according to layer widths and the proportion of fine roots (Sperry et al. 2016). Functions \code{root_horizontalProportionsBasic} and \code{root_horizontalProportionsAdvanced} calculate the (horizontal) proportion of roots of each cohort in the water pool corresponding to itself and that of other cohorts, for each soil layer.
}
\value{
  Functions \code{root_conicDistribution} and \code{root_ldrDistribution} return a matrix with as many rows as elements in \code{Z} (or \code{Z50}) and as many columns as soil layers.  Values in all cases correspond to the proportion of fine roots in each soil layer. Function \code{root_xylemConductanceProportions} returns a vector of proportions of the same length as the number of layers. Functions \code{root_horizontalProportionsBasic} and \code{root_horizontalProportionsAdvanced} return a list (with as many elements as cohorts) with each element being a matrix.
}
\references{
Schenk, H., Jackson, R., 2002. The global biogeography of roots. Ecol. Monogr. 72, 311–328.

Sperry, J. S., Y. Wang, B. T. Wolfe, D. S. Mackay, W. R. L. Anderegg, N. G. Mcdowell, and W. T. Pockman. 2016. Pragmatic hydraulic theory predicts stomatal responses to climatic water deficits. New Phytologist 212, 577–589.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{
 \code{\link{spwb}}, \code{\link{spwb_ldrOptimization}}, \code{\link{forest2spwbInput}}, \code{\link{soil}}
}
\examples{
#Load example plot plant data
data(exampleforestMED)

#Default species parameterization
data(SpParamsMED)

ntree = nrow(exampleforestMED$treeData)

#Initialize soil with default soil params
S = soil(defaultSoilParams())

#Calculate conic root system for trees
V1 = root_conicDistribution(Z=rep(2000,ntree), S$dVec)            
print(V1)

#Calculate LDR root system for trees (Schenck & Jackson 2002)
V2 = root_ldrDistribution(Z50 = rep(200,ntree), 
          Z95 = rep(1000,ntree), S$dVec)
print(V2)     

#Equal xylem conductance proportions for a cone distribution
#(assuming depth-width ratio 1)
root_xylemConductanceProportions(V1[1,], S$dVec)

#Xylem conductance proportions for LDR distribution
root_xylemConductanceProportions(V2[1,], S$dVec)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.