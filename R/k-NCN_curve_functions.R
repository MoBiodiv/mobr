centroid_accumulate<-function(x,focal_sample= 1, n=NULL, coords=NULL){
  require(rgeos)
  require(maptools)
  require(pracma)
  require(sp)
  if(class(x)== "mob_in"){
      sites=x$comm
      coords=x$spat
  }else{sites=x}
  sites<-as.data.frame(sites)
  rownames(sites)=NULL
  rownames(coords)=NULL
  if(is.null(n)) n=dim(sites)[1]
  included=focal_sample
  S_accumulated=as.numeric()
  for(i in 1:n){
    accumulated_plots<-sites[included,]
    coordinates(accumulated_plots)<-coords[included,]
    S_accumulated[i]<-vegan::specnumber(colSums(sites[included,,drop=F]))
    centroid<-gCentroid(accumulated_plots)@coords
    candidates<-coords[-included,]
    if(dim(candidates)[1]==0) break
    closest<-candidates[order(distmat(centroid, as.matrix(candidates)),  runif(dim(candidates)[1]) )[1],]
    included=c(included, as.numeric(rownames(closest)))
  }
  return(S_accumulated)
}

kNCN_average<-function(x, n=NULL, coords=NULL){
  if(class(x)== "mob_in"){
    sites=x$comm
    coords=x$spat
  }else{sites=x}
  rownames(sites)=NULL
  rownames(coords)=NULL
  if(is.null(n)) n=dim(sites)[1]
  samples=1:n
  require(pbapply)
  out<-pbsapply(samples,function(sample)  mobr:::centroid_accumulate(x = x,focal_sample = sample, coords = coords) )
  return(rowMeans(out))
}


compare_curves<-function(x){
  
  centroid_curve<-kNCN_average(x)
  spatcurve<-rarefaction(x,method = "spat")
  samplecurve<-rarefaction(x,method = "samp")
  plot(samplecurve, type= "l", xlab="Samples", ylab="Expected species richness")
  lines(centroid_curve,col=2)
  lines(spatcurve, col=3)
  legend("bottomright", legend = c("sample based", "k-NCN","k-NN"), col=1:3, lty=1)
}

