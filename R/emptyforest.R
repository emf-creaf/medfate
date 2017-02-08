emptyforest<-function(ID="", patchsize=10000) {
  l = list()
  l$ID = ID
  l$patchsize=10000
  l$treeData = data.frame(Species=numeric(0),DBH=numeric(0), 
                         Height=numeric(0), N=numeric(0),
                         Z50 = numeric(0), Z95=numeric(0))
  l$shrubData = data.frame(Species=numeric(0), Height=numeric(0), 
                          Cover = numeric(0), Z = numeric(0))
  l$seedBank = data.frame(Species=numeric(0),Abundance=numeric(0))
  l$herbCover = 0;
  l$herbHeight = 0;
  class(l)<-c("forest","list")
  return(l)
}