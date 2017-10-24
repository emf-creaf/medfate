\encoding{UTF-8}
\name{swb.ldrOptimization}
\alias{swb.ldrOptimization}

\title{Optimization of root distribution}

\description{The function \code{swb.ldrOptimization} optimizes the species root distribution within \code{swb}, given the arguments \code{x}, \code{meteo}, \code{soil} and \code{psi_crit}. The optimization is based on the eco-hydrological equilibrium hypothesis (Eagleson, 1982), which is formulated here as the root distribution for which plant transpiration is maximized while the plant water potential is close to the species-defined critical value \code{psi_crit} (Cabon et al., under review). Returns the optimized root distribution for each species, expressed as parameters of the function \code{root.ldrDistribution}.}

\usage{
swb.ldrOptimization(x, soil, meteo, psi_crit,  RZmin = 301, RZmax = 4000, 
                    V1min = 0.01, V1max = 0.94, resolution = 20, heat_stop = 0, 
                    transformation = "identity", explore_out = FALSE, 
                    verbose = FALSE)
}

\arguments{
  \item{x}{An object of class \code{\link{swbInput}}.}
  \item{soil}{A list containing the description of the soil (see \code{\link{soil}}).}
  \item{meteo}{A data frame with daily meteorological data series. When using the 'Simple' transpiration mode the following columns are required:
   \itemize{
     \item{\code{DOY}: Day of the year (Julian day).}
     \item{\code{Precipitation}: Precipitation (in mm).}
     \item{\code{MeanTemperature}: Mean temperature (in degrees Celsius).}
     \item{\code{PET}: Potential evapotranspiration (in mm).}
    }
  }
  \item{psi_crit}{A numerical vector of length iqual to the number of species in the plot containing the species values of water potential inducing hydraulic failure (in MPa). Use \code{NA} values to skip optimization for particular plant cohorts.}
  \item{RZmin}{The minimum value of RZ (the rooting depth) to be explored (in mm)}
  \item{RZmax}{The maximum value of RZ (the rooting depth) to be explored (in mm)}
  \item{V1min}{The minimum value of V1 (the root proportion in the first soil layer) to be explored}
  \item{V1max}{The maximum value of V1 (the root proportion in the first soil layer) to be explored}
  \item{resolution}{An integer defining the number of values to obtain by discretization of the root parameters RZ and V1. The number of parameter combinations and therefore the computation cost increases increase with the square of resolution}
  \item{transformation}{Function to modify the size of Z intervals to be explored (by default, bins are equal).}
  \item{heat_stop}{An integer defining the number of days during to discard from the calculation of the optimal root distribution. Usefull if the soil water content initialization is not certain}
  \item{explore_out}{A logical value. Are the values of average daily plant transpiration and minimum plant water potential to be returned?}
  \item{verbose}{A logical value. Print the internal messages of the function?}
}
\details{
For each combination of the parameters RZ and V1 the function runs \code{swb}, setting the total soil depth equal to RZ. The root proportion in each soil layer is derived from V1, the depth of the first soil layer (as defined in the argument \code{soil}) and RZ using the LDR root distribution model (Schenk and Jackson, 2002) and assuming that the depth containing 95 percent of the roots is equal to RZ. \code{psi_crit}, the species specific water potential inducing hydraulic failure, can be approached by the water potential inducing 50 percent of loss of conductance for the and gymnosperms and 88 percent for the angiosperms (Urli et al., 2013, Brodribb et al., 2010). Details of the hypothesis and limitations of the optimization method are given in Cabon et al. (under review).
}
\value{
If \code{explore_out = FALSE} returns a data frame with containing the species index used in medfate, \code{psi_crit} and the otpimized values of V1 and the LDR parameters Z50 and Z95 (see \code{\link{root.ldrDistribution}}) and as many rows as the number of species. If \code{explore_out = TRUE} returns a list containing \code{optim}, the aforementioned data frame, and \code{explore_out}, a list containing two array with the values of average daily plant transpiration and the minimum plant water potential for each species and parameter combination.
}
\references{
Brodribb, T.J., Bowman, D.J.M.S., Nichols, S., Delzon, S., Burlett, R., 2010. Xylem function and growth rate interact to determine recovery rates after exposure to extreme water deficit. New Phytol. 188, 533–542. doi:10.1111/j.1469-8137.2010.03393.x

Cabon, A., \enc{Martínez-Vilalta}{Martinez-Vilalta}, J., Poyatos, R., \enc{Martínez de Aragón}{Martinez de Aragon}, J., De \enc{Cáceres}{Caceres}, M., Applying the eco-hydrological equilibrium hypothesis to estimate root ditribution in water-limited forests. Under review in Agric. For. Meteorol.

Eagleson, P.S., 1982. Ecological optimality in water-limited natural soil-vegetation systems: 1. Theory and hypothesis. Water Resour. Res. 18, 325–340. doi:10.1029/WR018i002p00325

Schenk, H.J., Jackson, R.B., 2002. The Global Biogeography of Roots. Ecol. Monogr. 72, 311. doi:10.2307/3100092

Urli, M., Porte, A.J., Cochard, H., Guengant, Y., Burlett, R., Delzon, S., 2013. Xylem embolism threshold for catastrophic hydraulic failure in angiosperm trees. Tree Physiol. 33, 672–683. doi:10.1093/treephys/tpt030
}
\author{
Antoine Cabon, Centre Tecnologic Forestal de Catalunya

Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
 \code{\link{swb}}, \code{\link{soil}}, \code{\link{root.ldrDistribution}}
}

\examples{
\dontrun{
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Initialize soil with default soil params
examplesoil = soil(defaultSoilParams(2))

#Initialize control parameters
control = defaultControl()

#Initialize input
x = forest2swbInput(exampleforest,examplesoil, SpParamsMED, control)

#Run optimization
swb.ldrOptimization(x = x, soil = examplesoil, meteo = examplemeteo, psi_crit = c(-2,-3,-4,-4))
}

}

