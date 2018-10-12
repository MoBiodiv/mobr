centroid_accumulate<-function(x,focal_sample= 1, n=dim(x$comm)[1]){
  require(rgeos)
  require(maptools)
  sites<-x$comm
  included=focal_sample
  S_accumulated=as.numeric()
  for(i in 1:n){
    accumulated_plots<-sites[included,]
    coordinates(accumulated_plots)<-x$spat[included,]
    S_accumulated[i]<-vegan::specnumber(colSums(sites[included,,drop=F]))
    centroid<-gCentroid(accumulated_plots)@coords
    rownames(centroid)<-"centroid"
    candidates<-x$spat[-included,]
    
    closest<-candidates[order(as.matrix(dist(rbind(centroid,candidates)))[-1,1],  runif(dim(candidates)[1]) )[1],]
    if(rownames(closest)=="NA") break
    included=c(included, as.numeric(rownames(closest)))
  }
  return(S_accumulated)
}

kNCN_average<-function(x, n=dim(x$comm)[1]){
  samples=1:n
  require(pbapply)
  out<-pbsapply(samples,function(sample)  centroid_accumulate(x = x,focal_sample = sample,n = n) )
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
