\encoding{UTF-8}
\name{exampleforest}
\Rdversion{1.1}
\alias{exampleforest}
\alias{exampleforest2}
\docType{data}
\title{Example forest stands}
\description{
  Data set for illustration of model behaviour. Includes a description of the plant cohorts of a forest stand.
}
\usage{
data(exampleforest)
data(exampleforest2)
}
\format{
 An object of class \code{\link{forest}} containing the description of the woody (tree or shrub) cohorts and herb layer of a forest patch. \code{exampleforest} represents the standard forest inventory description, whereas \code{exampleforest2} is an alternative forest description where leaf area index and crown ratio are supplied instead of structural (density, DBH and cover) parameters. 
}
\source{
DGCN (2005). Tercer Inventario Forestal Nacional (1997-2007): Catalunya. \enc{Dirección}{Direccion} General de \enc{Conservación}{Conservacion} de la Naturaleza, Ministerio de Medio Ambiente, Madrid.
}
\seealso{ \code{\link{forest}}, \code{\link{spwb}}, \code{\link{spwbInput}}}
\examples{
data(exampleforest)
data(exampleforest2)
}
\keyword{datasets}
