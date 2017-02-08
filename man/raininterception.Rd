\encoding{UTF-8}
\name{swb.RainInterception}
\alias{swb.RainInterception}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rainfall interception}
\description{
Function \code{swb.RainInterception} calculates the amount of rainfall intercepted daily by the canopy, given a rainfall and canopy characteristics. Two canopy interception models are currently available: the sparse Gash (1995) model and the Liu (2001) model. In both cases the current implementation assumes no trunk interception.
}
\usage{
swb.RainInterception(Rainfall, Cm, p, ER=0.05, method="Gash1995") 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Rainfall}{A numeric vector of (daily) rainfall.}
  \item{Cm}{Canopy water storage capacity.}
  \item{p}{Proportion of throughfall (normally 1 - c, where c is the canopy cover).}
  \item{ER}{The ratio of evaporation rate to rainfall rate.}
  \item{method}{Rainfall interception method (either \code{"Gash1995"} or \code{Liu2001}).}
}
\details{
The function is prepared to accept both vectors or scalars for parameters \code{Cm}, \code{p} and \code{ER}. If they are supplied as vectors they should be of the same length as \code{Rainfall}.
}
\value{
Returns a vector of the same length as \code{Rainfall} containing intercepted rain values.
}
\references{
Liu (2001). Evaluation of the Liu model for predicting rainfall interception in forests world-wide. - Hydrol. Process. 15: 2341-2360.

Gash (1979). An analytical model of rainfall interception by forests. - Quarterly Journal of the Royal Meteorological Society                                       

Gash et al. (1995). Estimating sparse forest rainfall interception with an analytical model. - Journal of Hydrology  
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
 \code{\link{swb}}
}
\examples{
throughfallMatrixGash<-function(P = seq(1,50, by=1), Cm = seq(1,5, by=1), 
                                ER = 0.08,p=0.8) {
  m2<-P-swb.RainInterception(P,Cm[1],p,ER=ER)
  for(i in 2:length(Cm)) {
    m2<-rbind(m2,P-swb.RainInterception(P,Cm[i],p,ER=ER))
  }
  colnames(m2)<-P
  rownames(m2)<-Cm
  return(m2)
}

Cm = c(0.5,seq(1,4, by=1))
P = seq(1,50, by=1)

m2 = throughfallMatrixGash(P=P, p=0.2, Cm=Cm,ER = 0.05)
rt = sweep(m2,2,P,"/")*100
matplot(t(rt), type="l", axes=TRUE, ylab="Relative throughfall (\%)", 
        xlab="Gross rainfall (mm)", xlim=c(0,length(P)), 
        lty=1:length(Cm), col="black", ylim=c(0,100))
title(main="p = 0.2 E/R = 0.05")
legend("bottomright",lty=1:length(Cm), legend=paste("Cm =",Cm), bty="n")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.