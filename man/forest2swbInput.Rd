\encoding{UTF-8}
\name{swbInput}
\alias{forest2aboveground}
\alias{forest2growthInput}
\alias{forest2swbInput}
\alias{growthInput}
\alias{swbInput}

\title{
Input for simulation models
}

\description{
Functions \code{forest2swbInput} and \code{forest2growthInput} take an object of class \code{\link{forest}} and calculate input data for functions \code{\link{swb}} and \code{\link{growth}}, respectively. Functions \code{swbInput} and \code{growthInput} does the same from input data. Function \code{forest2aboveground} calculates aboveground variables that may be used in \code{swbInput} and \code{growthInput} functions.
}
\usage{
forest2aboveground(x, SpParams, gdd = NA)
forest2growthInput(x, soil, SpParams, control)
forest2swbInput(x, soil, SpParams, control)
growthInput(above,  Z, V,  soil, SpParams, control)
swbInput(above,  V,  soil, SpParams, control)
}
\arguments{
  \item{x}{An object of class \code{\link{forest}}.}
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}} and \code{\link{SpParamsMED}}).}
  \item{gdd}{Growth degree days to account for leaf phenology effects (in Celsius). This should be left \code{NA} in most applications.}
  \item{soil}{An object of class \code{\link{soil}}.}
  \item{control}{A list with default control parameters (see \code{\link{defaultControl}}).}
  \item{above}{A data frame with aboveground plant information (see the return value of \code{forest2aboveground} below). In the case of \code{swbInput} the variables should include \code{SP}, \code{LAI_live}, \code{LAI_dead}, \code{H} and \code{CR}.}
  \item{V}{A numeric matrix (with as many columns as soil layers and as many rows as the length as \code{SP}) containing the proportion of roots of each plant in each soil layer.}
  
}
\details{
Functions \code{forest2swbInput} and \code{forest2abovegroundInput} extracts height and species identity from plant cohorts of \code{x}, and calculate leaf area index and crown ratio.\code{forest2swbInput} also calculates the distribution of fine roots across soil. Both \code{forest2swbInput} and \code{swbInput} find parameter values for each plant cohort according to the parameters of its species as specified in \code{SpParams}. If \code{control$transpirationMode = "Complex"} the functions also estimate the maximum conductance of rhizosphere, root xylem and stem xylem elements.
}
\value{
Function \code{forest2aboveground} returns a data frame with the following columns:
\itemize{
  \item{\code{SP}: Species identity (an integer) (first species is 0).}
  \item{\code{N}: Cohort density (ind/ha) (see function \code{\link{plant.Density}}).}
  \item{\code{DBH}: Tree diameter at breast height (cm).}
  \item{\code{H}: Plant total height (cm).}
  \item{\code{CR}: Crown ratio (crown length to total height) (between 0 and 1).}
  \item{\code{LAI_live}: Live leaf area index (m2/m2) (one-side leaf area relative to plot area).}
  \item{\code{LAI_dead}: Dead leaf area index (m2/m2) (one-side leaf area relative to plot area).}
}

Functions \code{forest2swbInput} and \code{swbInput} return a list of class \code{swbInput} with the following elements:
\itemize{
  \item{\code{control}: List with control parameters (see \code{\link{defaultControl}}).}
  \item{\code{above}: A data frame with columns \code{SP}, \code{H}, \code{CR} and \code{LAI} (see function \code{forest2aboveground}).}
  \item{\code{below}: A list. If \code{control$transpirationMode = "Simple"} it contains a single element: 
     \itemize{
        \item{\code{V}: A matrix with the proportion of fine roots of each cohort (in rows) in each soil layer (in columns).}
     }
     If \code{control$transpirationMode = "Complex"} there are the following additional elements:
     \itemize{
        \item{\code{VGrhizo_kmax}: A matrix with maximum rhizosphere conductance values of each cohort (in rows) in each soil layer (in columns).}
        \item{\code{VGroot_kmax}: A matrix with maximum root xylem conductance values of each cohort (in rows) in each soil layer (in columns).}
     }
  }
  \item{\code{paramsBase}: A data frame with columns:
     \itemize{
       \item{\code{k}: PAR extinction coefficient.}
       \item{\code{g}: Canopy water retention capacity per LAI unit (mm/LAI).}
       \item{\code{Sgdd}: Growth degree days needed for leaf budburst (for winter decideous species).}
     }
  }
  \item{\code{paramsTransp}: A data frame with transpiration parameters. If \code{control$transpirationMode = "Simple"}, columns are:
    \itemize{
       \item{\code{Psi_Extract}: Water potential corresponding to 50\% relative conductance (in MPa).}
       \item{\code{WUE}: Water use efficiency for carbon assimilation (g C/mm water).} 
    }
    If \code{control$transpirationMode = "Complex"} columns are:
    \itemize{
       \item{\code{Gwmax}: Maximum stomatal conductance to water vapor (in mol H2O·m-2·s-1).}
       \item{\code{Vmax298}: Maximum Rubisco carboxilation rate at 25ºC (in micromol CO2·s-1·m-2).}
       \item{\code{Jmax298}: Maximum rate of electron transport at 25ºC (in micromol photons·s-1·m-2).}
       \item{\code{VCroot_c}, \code{VCroot_d}: Parameters of the root xylem vulnerability curve.}
       \item{\code{xylem_kmax}: Sapwood-specific hydraulic conductivity of stem xylem (in kg H2O·s-1·m-2).}
       \item{\code{VCstem_kmax}: Maximum stem xylem conductance values of each cohort.}
       \item{\code{VCstem_c}, \code{VCstem_d}: Parameters of the stem xylem vulnerability curve.}
    }
  }
  \item{\code{Transpiration}: Plant cohort transpiration of the current day (mm of water = L/m2; filled with zeroes before simulations).}
  \item{\code{Photosynthesis}: Plant cohort photosynthesis of the current day (g C/m2; filled with zeroes before simulations).}
}
Functions \code{forest2growthInput} and \code{growthInput} return a list of class \code{growthInput} with the same elements as \code{swbInput}, but with additional information. 

\itemize{
\item{Element \code{above} includes the following additional columns:
  \itemize{
    \item{\code{LA_live}: Live leaf area per individual (m2/ind).}
    \item{\code{LA_dead}: Dead leaf area per individual (m2/ind).}
    \item{\code{SA}: Live sapwood area per individual (cm2/ind).} 
  }
  }
  \item{\code{paramsGrowth}: A data frame with columns:
     \itemize{
       \item{\code{SLA}:Specific leaf area (mm2/mg = m2/kg).}
       \item{\code{Al2As}: Leaf area to sapwood area ratio (in m2·m-2).}
     }
  }
}

}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{\code{\link{swb}}, \code{\link{soil}},  \code{\link{forest}}, \code{\link{SpParamsMED}}, \code{\link{defaultSoilParams}}}

\examples{
#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Aboveground parameters
forest2aboveground(exampleforest, SpParamsMED)

#Initialize soil with default soil params
examplesoil = soil(defaultSoilParams())

#Initialize control parameters
control = defaultControl()

#Prepare swb input for simple hydraulic mode
forest2swbInput(exampleforest,examplesoil,SpParamsMED, control)

#Prepare growth input for simple hydraulic mode
forest2growthInput(exampleforest,examplesoil,SpParamsMED, control)


#Prepare input for Complex hydraulic mode
control$transpirationMode = "Complex"

forest2swbInput(exampleforest,examplesoil,SpParamsMED, control)
}
