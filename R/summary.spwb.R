.mergeVectorOfMatrices<-function(vec){
  out<- NULL
  for(i in 1:length(vec)) {
    out_i <- vec[[i]]
    if(is.null(out)) out = out_i
    else {
      cno <- colnames(out)
      cno_i <- colnames(out_i)
      if(length(cno)==length(cno_i) && sum(cno==cno_i)==length(cno)) {
        out<-rbind(out, out_i)
      } else {
        cn_all = unique(c(cno, cno_i))
        cn_old = cno[which(!(cno %in% cno_i))]
        cn_new = cno_i[which(!(cno_i %in% cno))]
        if(length(cn_new)>0) {
          out <- as.data.frame(out)
          for(n in cn_new) out[[n]] = NA
        }
        if(length(cn_old)>0) {
          out_i <- as.data.frame(out_i)
          for(n in cn_old) out_i[[n]] = NA
        }
        out <- as.matrix(out[,cn_all, drop= FALSE])
        out_i <- as.matrix(out_i[,cn_all, drop= FALSE])
        out<-rbind(out, out_i)
      }
    }
  }
  return(out)
}
.summarysim<-function(object, freq="years", output="WaterBalance", FUN=sum, bySpecies = FALSE,
                      months = NULL, ...){  
  dates <- as.Date(rownames(object$WaterBalance))
  month_dates <- as.character(as.numeric(format(dates, "%m")))
  ndaysTotal <- length(dates)
  date.factor <- cut(dates, breaks=freq)
  
  if("spwbInput" %in% names(object)) input = object$spwbInput
  else input = object$growthInput
  object_names = names(object)
  output_vec = strsplit(output, "\\$")[[1]]
  if(!(output_vec[1] %in% object_names)) {
    # Try to complete
    found = FALSE
    object_names_search = object_names
    object_names_search = object_names_search[!(object_names_search %in% c("Soil", "Stand", "Temperature",
                                                                           "CarbonBalance", "WaterBalance", "EnergyBalance", "BiomassBalance",
                                                                           "FireHazard"))]
    for(nm in object_names_search) {
      if(!found) {
        if(output_vec[1] %in% names(object[[nm]])) {
          output_vec = c(nm, output_vec[1])
          found = TRUE
        }
      }
    }
    if(!found) {
      stop(paste0("Output table '", output_vec[1], "' was not generated in simulation results"))
    }
  } else if(length(output_vec)==2) {
    if(!(output_vec[2] %in% names(object[[output_vec[1]]]))) stop(paste0("Unrecognized output string: '", output ,"'\n"))
  } else if(length(output_vec)==1 && output_vec[1]=="LabileCarbonBalance") {
    output_vec = rep(output_vec, 2)
  }
  
  if(output_vec[1] %in% c("Soil", "Stand", "Temperature", 
                          "CarbonBalance", "WaterBalance", "EnergyBalance", "BiomassBalance",
                          "FireHazard")) {
    OM = object[[output_vec[1]]]
  } else if(output_vec[1]=="Plants" && output_vec[2]=="LAI") {
    OM = object$Plants$LAI
    if(bySpecies) {
      OM = t(apply(OM,1, tapply, as.factor(input$cohorts$Name), sum, na.rm=T))
      if(nrow(OM)==1) rownames(OM) = levels(as.factor(input$cohorts$Name))
    } 
  } else {
    OM = object[[output_vec[1]]][[output_vec[2]]]
    if(bySpecies && ncol(OM)>0) {
      lai1 = t(apply(object$Plants$LAI,1, tapply, as.factor(input$cohorts$Name), sum, na.rm=T))
      if(nrow(lai1)==1) rownames(lai1) = levels(as.factor(input$cohorts$Name))
      m1 = t(apply(object$Plants$LAI * OM,1, tapply, as.factor(input$cohorts$Name), sum, na.rm=T))
      if(nrow(m1)==1) rownames(m1) = levels(as.factor(input$cohorts$Name))
      OM = m1/lai1
      OM[lai1==0] = NA
    } else if(bySpecies) {
      colnames(OM) = levels(as.factor(input$cohorts$Name))
    }
  }
  if(ncol(OM)==length(date.factor) && nrow(OM)==1) OM = t(OM)
  
  # Subset dates if 'months' is specified
  if(!is.null(months)) {
    months <- match.arg(as.character(months), as.character(1:12), several.ok = TRUE)
    sel_months <- month_dates %in% months
    date.factor <- date.factor[sel_months]
    OM <- OM[sel_months,,drop = FALSE]
  }
  #Perform summary at the desired temporal scale
  M <- apply(OM,2,tapply, INDEX=date.factor, FUN, na.rm=T)
  if(length(M)==0) {
    l <- levels(date.factor)
    M <- matrix(nrow = length(l), ncol=0)
    rownames(M) <- l
  } else if(is.vector(M)) {
    M = t(as.matrix(M))
    rownames(M) <- levels(date.factor)
  }
  # if(sum(is.na(M[nrow(M), drop=FALSE]))==ncol(M)) M = M[-nrow(M), drop=FALSE] #Remove empty row
  ncases = table(date.factor)
  if(length(M)>0) {
    M = M[ncases>0, ,drop = FALSE]
  }
  return(M)
}

#' Summarize simulation results
#' 
#' Function \code{summary} summarizes the model's output in different temporal steps (i.e. weekly, annual, ...).
#' 
#' @param object An object of class \code{spwb}, \code{pwb}, \code{growth} or \code{fordyn}.
#' @param freq Frequency of summary statistics (see \code{\link{cut.Date}}).
#' @param output The data table to be summarized. Accepted values are the path to data tables in \code{object}, such as 'WaterBalance', 'Soil', 'Stand' or 'Plants$LAI'. It is also possible to use strings like 'Transpiration' and the function will interpret it as 'Plants$Transpiration'.
#' @param FUN The function to summarize results (e.g., \code{sum}, \code{mean}, ...)
#' @param bySpecies Allows aggregating output by species before calculating summaries (only has an effect with some values of \code{output}). Aggregation can involve a sum (as for plant lai or transpiration) or a LAI-weighted mean (as for plant stress or plant water potential).
#' @param months A vector of month numbers (1 to 12) to subset the season where summaries apply.
#' @param ... Additional parameters for function \code{summary}.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @return A matrix with dates as row names and the desired summaries in columns
#' 
#' @note When applied to \code{\link{fordyn}} objects, the summary function can be used to gather the results of different yearly steps into a single table while keeping a daily resolution (i.e. using \code{freq = "days"}.
#' 
#' @seealso \code{\link{spwb}}, \code{\link{pwb}}, \code{\link{growth}}, \code{\link{fordyn}}, \code{\link{plot.spwb}}, \code{\link{extractSubdaily}}
#' 
#' @examples
#' #Load example daily meteorological data
#' data(examplemeteo)
#' 
#' #Load example plot plant data
#' data(exampleforestMED)
#' 
#' #Default species parameterization
#' data(SpParamsMED)
#' 
#' #Initialize soil with default soil params (2 layers)
#' examplesoil = soil(defaultSoilParams(2))
#' 
#' #Initialize control parameters
#' control = defaultControl("Granier")
#' 
#' #Initialize input
#' x = forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)
#' 
#' #Call simulation function
#' S1<-spwb(x, examplemeteo, latitude = 41.82592, elevation = 100)
#' 
#' #Monthly summary (averages) of soil status
#' summary(S1, freq="months",FUN=mean, output="Soil")
#' 
#' #Queries the tables in 'Plants'
#' names(S1$Plants)
#' 
#' #Monthly summary (averages) of plant stress
#' summary(S1, freq="months",FUN=mean, output="Plants$PlantStress", 
#'         bySpecies = TRUE)
#' 
#' 
#' @name summary.spwb
summary.spwb<-function(object, freq="years", output="WaterBalance", FUN=sum, bySpecies = FALSE, months = NULL, ...){  
  .summarysim(object = object, freq = freq, output = output, FUN = FUN, bySpecies = bySpecies, months = months, ...)
}

#' @rdname summary.spwb
summary.pwb<-function(object, freq="years", output="WaterBalance", FUN=sum, bySpecies = FALSE,  months = NULL,...){  
  .summarysim(object = object, freq = freq, output = output, FUN = FUN, bySpecies = bySpecies,  months = months,...)
}

#' @rdname summary.spwb
summary.growth<-function(object, freq="years", output="WaterBalance", FUN=sum, bySpecies = FALSE,  months = NULL,...){  
  .summarysim(object = object, freq = freq, output = output, FUN = FUN, bySpecies = bySpecies, months = months, ...)
}

#' @rdname summary.spwb
summary.fordyn<-function(object, freq="years", output="WaterBalance", FUN=sum, bySpecies = FALSE,  months = NULL, ...){
  vec<-vector("list", length(object$GrowthResults))
  for(i in 1:length(object$GrowthResults)) {
    vec[[i]] <- .summarysim(object = object$GrowthResults[[i]], 
                           freq = freq, output = output, FUN = FUN, bySpecies = bySpecies, months = months, ...)
  }
  return(.mergeVectorOfMatrices(vec))
}