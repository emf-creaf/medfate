.summarysim<-function(object, freq="years", output="WaterBalance", FUN=sum, bySpecies = FALSE, ...){  
  dates = as.Date(rownames(object$WaterBalance))
  ndaysTotal = length(dates)
  date.factor = cut(dates, breaks=freq)
  input = object$spwbInput
  object_names = names(object)
  output_vec = strsplit(output, "\\$")[[1]]
  if(!(output_vec[1] %in% object_names)) {
    # Try to complete
    found = FALSE
    object_names_search = object_names
    object_names_search = object_names_search[!(object_names_search %in% c("Soil", "Stand", "Temperature",
                                                                           "WaterBalance", "EnergyBalance"))]
    for(nm in object_names_search) {
      if(!found) {
        if(output_vec[1] %in% names(object[[nm]])) {
          output_vec = c(nm, output_vec[1])
          found = TRUE
        }
      }
    }
    if(!found) {
      stop(paste0("Unrecognized output string: '", output ,"'\n"))
    }
  }
  else if(length(output_vec)==2) {
    if(!(output_vec[2] %in% names(object[[output_vec[1]]]))) stop(paste0("Unrecognized output string: '", output ,"'\n"))
  }
  
  if(output_vec[1] %in% c("Soil", "Stand", "Temperature", 
                          "WaterBalance", "EnergyBalance")) {
    OM = object[[output_vec[1]]]
  }
  else if(output_vec[1]=="Plants" && output_vec[2]=="LAI") {
    OM = object$Plants$LAI
    if(bySpecies) {
      OM = t(apply(OM,1, tapply, input$cohorts$Name, sum, na.rm=T))
    } 
  } 
  else {
    OM = object[[output_vec[1]]][[output_vec[2]]]
    if(bySpecies) {
      lai1 = t(apply(object$Plants$LAI,1, tapply, input$cohorts$Name, sum, na.rm=T))
      m1 = t(apply(object$Plants$LAI * OM,1, tapply, input$cohorts$Name, sum, na.rm=T))
      OM = m1/lai1
      OM[lai1==0] = NA
    } 
  }
  if(ncol(OM)==length(date.factor) && nrow(OM)==1) OM = t(OM)
  
  #Perform summary at the desired temporal scale
  M <- apply(OM,2,tapply, INDEX=date.factor, FUN)
  
  if(is.vector(M)) {
    M = t(as.matrix(M))
    rownames(M) <- levels(date.factor)
  }
  # if(sum(is.na(M[nrow(M), drop=FALSE]))==ncol(M)) M = M[-nrow(M), drop=FALSE] #Remove empty row
  ncases = table(date.factor)
  M = M[ncases>0, ,drop = FALSE]
  
  return(M)
}

summary.spwb<-function(object, freq="years", output="WaterBalance", FUN=sum, bySpecies = FALSE, ...){  
  .summarysim(object = object, freq = freq, output = output, FUN = FUN, bySpecies = bySpecies, ...)
}

summary.pwb<-function(object, freq="years", output="WaterBalance", FUN=sum, bySpecies = FALSE, ...){  
  .summarysim(object = object, freq = freq, output = output, FUN = FUN, bySpecies = bySpecies, ...)
}

summary.growth<-function(object, freq="years", output="WaterBalance", FUN=sum, bySpecies = FALSE, ...){  
  .summarysim(object = object, freq = freq, output = output, FUN = FUN, bySpecies = bySpecies, ...)
}

summary.fordyn<-function(object, freq="years", output="WaterBalance", FUN=sum, bySpecies = FALSE, ...){
  out<-NULL
  for(i in 1:length(object$GrowthResults)) {
    out_i <- .summarysim(object = object$GrowthResults[[i]], 
                     freq = freq, output = output, FUN = FUN, bySpecies = bySpecies, ...)
    if(is.null(out)) out = out_i
    else {
      out<-rbind(out, out_i)
    }
  }
  return(out)
}