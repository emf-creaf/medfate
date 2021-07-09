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
A data sets of species parameter definition and values, the latter resulting from existing databases, fit to empirical data or expert-based guesses.
}
\usage{
data("SpParamsDefinition")
data("SpParamsMED")
data("SpParamsUS")
}
\format{
Data frame \code{SpParamsDefinition} has parameters in rows and columns 'ParameterName', 'ParameterGroup', 'Definition', 'Type' and 'Units'. Data frames \code{SpParamsMED} and \code{SpParamsUS} have species or genus as rows and column names equal to parameter names in \code{SpParamsDefinition}.
}
\details{
Plant trait parameter sources are listed in the bibliography section. Details of the procedures used to obtain the species parameter tables can be found in the article 'Species parameterization for Spain' at https://emf-creaf.github.io/medfate/. 
}
\source{
Asse, D., Randin, C. F., Bonhomme, M., Delestrade, A., & Chuine, I. (2020). Process-based models outcompete correlative models in projecting spring phenology of trees in a future warmer climate. Agricultural and Forest Meteorology, 285–286(January), 107931. https://doi.org/10.1016/j.agrformet.2020.107931

Bartlett, Megan K, Christine Scoffoni, and Lawren Sack. 2012. The determinants of leaf turgor loss point and prediction of drought tolerance of species and biomes: a global meta-analysis. Ecology Letters 15 (5): 393–405. https://doi.org/10.1111/j.1461-0248.2012.01751.x. 

Choat, Brendan, Steven Jansen, Tim J Brodribb, \enc{Hervé}{Herve} Cochard, Sylvain Delzon, Radika Bhaskar, Sandra J Bucci, et al. 2012. Global convergence in the vulnerability of forests to drought. Nature 491 (7426): 752–55. https://doi.org/10.1038/nature11688.

De \enc{Cáceres}{Caceres} M, Casals P, Gabriel E, Castro X (2019) Scaling-up individual-level allometric equations to predict stand-level fuel loading in Mediterranean shrublands. Ann For Sci 76:87 (doi:10.1007/s13595-019-0873-4).

Delpierre, N., \enc{Dufrêne}{Dufrene}, E., Soudani, K., Ulrich, E., Cecchini, S., \enc{Boé}{Boe}, J., & \enc{François}{Francois}, C. (2009). Modelling interannual and spatial variability of leaf senescence for three deciduous tree species in France. Agricultural and Forest Meteorology, 149(6–7), 938–948. https://doi.org/10.1016/j.agrformet.2008.11.014

DGCN (2005). Tercer Inventario Forestal Nacional (1997-2007): Catalunya. Dirección General de Conservación de la Naturaleza, Ministerio de Medio Ambiente, Madrid.

Duursma, Remko A., Christopher J. Blackman, Rosana \enc{López}{Lopez}, Nicolas K. Martin-StPaul, \enc{Hervé}{Herve} Cochard, and Belinda E. Medlyn. 2018. On the minimum leaf conductance: its role in models of plant water use, and ecological and environmental controls. New Phytologist, no. September. https://doi.org/10.1111/nph.15395.

Hoshika, Yasutomo, Yutaka Osada, Alessandra de Marco, Josep \enc{Peñuelas}{Penuelas}, and Elena Paoletti. 2018. Global diurnal and nocturnal parameters of stomatal conductance in woody plants and major crops. Global Ecology and Biogeography 27 (2): 257–75. https://doi.org/10.1111/geb.12681.

Kattge, J., S. \enc{Díaz}{Diaz}, S. Lavorel, I. C. Prentice, P. Leadley, G. Bonisch, E. Garnier, et al. 2011. TRY - a global database of plant traits. Global Change Biology 17 (9): 2905–35. https://doi.org/10.1111/j.1365-2486.2011.02451.x.

Manzoni, Stefano, Giulia Vico, Amilcare Porporato, and Gabriel Katul. 2013. Biological constraints on water transport in the soil-plant-atmosphere system. Advances in Water Resources 51: 292–304. https://doi.org/10.1016/j.advwatres.2012.03.016.

Martin-StPaul, Nicolas, Sylvain Delzon, and \enc{Hervé}{Herve} Cochard. 2017. Plant resistance to drought depends on timely stomatal closure. Ecology Letters 20 (11): 1437–47. https://doi.org/10.1111/ele.12851.

Mencuccini, Maurizio, Teresa Rosas, Lucy Rowland, Brendan Choat, Hans Cornelissen, Steven Jansen, Koen Kramer, et al. 2019. Leaf economics and plant hydraulics drive leaf : wood area ratios. New Phytologist. https://doi.org/10.1111/nph.15998.

Pausas, Juli G., R. Brandon Pratt, Jon E. Keeley, Anna L. Jacobsen, Aaron R. Ramirez, Alberto Vilagrosa, Susana Paula, Iolana N. Kaneakua-Pia, and Stephen D. Davis. 2016. Towards understanding resprouting at the global scale. New Phytologist 209 (3): 945–54. https://doi.org/10.1111/nph.13644.

Sanchez‐Martinez, Pablo, Jordi Martinez‐Vilalta, Kyle G. Dexter, Ricardo A. Segovia, and Maurizio Mencuccini. 2020. Adaptation and coordinated evolution of plant hydraulic traits. Ecology Letters, August, ele.13584. https://doi.org/10.1111/ele.13584.

\enc{Tavşanoǧlu}{Tavsanoglu}, \enc{Çaǧatay}{Cagatay}, and Juli G. Pausas. 2018. A functional trait database for Mediterranean Basin plants. Scientific Data 5: 1–18. https://doi.org/10.1038/sdata.2018.135.

Vitasse, Y., \enc{François}{Francois}, C., Delpierre, N., \enc{Dufrêne}{Dufrene}, E., Kremer, A., Chuine, I., & Delzon, S. (2011). Assessing the effects of climate change on the phenology of European temperate trees. Agricultural and Forest Meteorology, 151(7), 969–980. https://doi.org/10.1016/j.agrformet.2011.03.003

Yebra, Marta, Gianluca Scortechini, Abdulbaset Badi, María Eugenia Beget, Matthias M. Boer, Ross Bradstock, Emilio Chuvieco, et al. 2019. Globe-LFMC, a global plant water status database for vegetation ecophysiology and wildfire applications. Scientific Data 6 (1): 155. https://doi.org/10.1038/s41597-019-0164-9.

Zanne, A. E., G. Lopez-Gonzalez, D. A. Coomes, J. Ilic, S. Jansen, S. L. Lewis, R. B. Miller, N. G. Swenson, M. C. Wiemann, and J. Chave. 2009. Global wood density database. http://datadryad.org/handle/10255/dryad.235. 
}
\seealso{\code{\link{spwb}}}
\examples{
data(SpParamsMED)
}
\keyword{datasets}
