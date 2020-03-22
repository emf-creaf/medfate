emptyforest<-function(ID="", patchsize=10000, ntree = 0, nshrub = 0) {
  l = list()
  l$ID = ID
  l$patchsize=10000
  l$treeData = data.frame(Species=numeric(ntree),DBH=numeric(ntree), 
                         Height=numeric(ntree), N=numeric(ntree),
                         Z50 = numeric(ntree), Z95=numeric(ntree))
  l$shrubData = data.frame(Species=numeric(nshrub), Height=numeric(nshrub), 
                          Cover = numeric(nshrub), 
                          Z50 = numeric(nshrub), Z95=numeric(nshrub))
  l$herbCover = 0;
  l$herbHeight = 0;
  class(l)<-c("forest","list")
  return(l)
}