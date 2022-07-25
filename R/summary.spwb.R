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
.summarysim<-function(object, freq="years", output="WaterBalance", FUN=sum, bySpecies = FALSE, ...){  
  dates = as.Date(rownames(object$WaterBalance))
  ndaysTotal = length(dates)
  date.factor = cut(dates, breaks=freq)
  if("spwbInput" %in% names(object)) input = object$spwbInput
  else input = object$growthInput
  object_names = names(object)
  output_vec = strsplit(output, "\\$")[[1]]
  if(!(output_vec[1] %in% object_names)) {
    # Try to complete
    found = FALSE
    object_names_search = object_names
    object_names_search = object_names_search[!(object_names_search %in% c("Soil", "Stand", "Temperature",
                                                                            "CarbonBalance", "WaterBalance", "EnergyBalance", "BiomassBalance"))]
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
  } else if(length(output_vec)==2) {
    if(!(output_vec[2] %in% names(object[[output_vec[1]]]))) stop(paste0("Unrecognized output string: '", output ,"'\n"))
  } else if(length(output_vec)==1 && output_vec[1]=="LabileCarbonBalance") {
    output_vec = rep(output_vec, 2)
  }
  
  if(output_vec[1] %in% c("Soil", "Stand", "Temperature", 
                           "CarbonBalance", "WaterBalance", "EnergyBalance", "BiomassBalance")) {
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
  vec<-vector("list", length(object$GrowthResults))
  for(i in 1:length(object$GrowthResults)) {
    vec[[i]] <- .summarysim(object = object$GrowthResults[[i]], 
                           freq = freq, output = output, FUN = FUN, bySpecies = bySpecies, ...)
  }
  return(.mergeVectorOfMatrices(vec))
}