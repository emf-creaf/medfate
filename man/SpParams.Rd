\encoding{UTF-8}
\name{SpParams}
\alias{SpParams}
\alias{SpParamsDefinition}
\alias{SpParamsMED}
\alias{SpParamsUS}
\docType{data}
\title{
Data tables with species parameter definition and values for Mediterranean and US
}
\description{
A data sets of species parameter definition and values, the latter resulting from bibliographic search, fit to empirical data or expert-based guesses.
}
\usage{
data("SpParamsDefinition")
data("SpParamsMED")
data("SpParamsUS")
}
\format{
  Data frame \code{SpParamsDefinition} contains parameters in rows and columns "ParameterName", "ParameterGroup", "Definition", "Type" and "Units". Data frames \code{SpParamsMED} and \code{SpParamsUS} have species or genus as rows and column names equal to "ParameterName" in \code{SpParamsDefinition}.
}
\details{
See details of parameterization in De Caceres et al. (2015) and De Caceres et al. (2021).
}
\source{
De \enc{Cáceres}{Caceres} M, \enc{Martínez}{Martinez}-Vilalta J, Coll L, Llorens P, Casals P, Poyatos R, Pausas JG, Brotons L. (2015) Coupling a water balance model with forest inventory data to predict drought stress: the role of forest structural changes vs. climate changes. Agricultural and Forest Meteorology (doi:10.1016/j.agrformet.2015.06.012).

De \enc{Cáceres}{Caceres} M, Mencuccini M, Martin-StPaul N, Limousin JM, Coll L, Poyatos R, Cabon A, Granda V, Forner A, Valladares F, \enc{Martínez}{Martinez}-Vilalta J (2021) Unravelling the effect of species mixing on water use and drought stress in holm oak forests: a modelling approach. Agricultural and Forest Meteorology 296 (doi:10.1016/j.agrformet.2020.108233).

De \enc{Cáceres}{Caceres} M, Casals P, Gabriel E, Castro X (2019) Scaling-up individual-level allometric equations to predict stand-level fuel loading in Mediterranean shrublands. Ann For Sci 76:87 (doi:10.1007/s13595-019-0873-4).


DGCN (2005). Tercer Inventario Forestal Nacional (1997-2007): Catalunya. Dirección General de Conservación de la Naturaleza, Ministerio de Medio Ambiente, Madrid.

}
\seealso{\code{\link{getSpParamsNames}},\code{\link{spwb}}}
\examples{
data(SpParamsDefinition)
data(SpParamsMED)
}
\keyword{datasets}
