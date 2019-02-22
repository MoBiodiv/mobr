# Internal function for kNCN_average (). k-NCN algorithm starting with a specified focal sample

centroid_accumulate<-function(x,focal_sample= 1, n=NULL, coords=NULL){
  require(rgeos)
  require(maptools)
  require(pracma)
  require(sp)
  if(class(x)== "mob_in"){
      sites=x$comm
      coords=x$spat
  }else{sites=x}
  sites<-as.matrix(sites)
  coords<-as.matrix(coords)
  rownames(sites)=NULL
  rownames(coords)=NULL
  if(is.null(n)) n=dim(sites)[1]
  included=focal_sample
  pool= colSums(sites[included,,drop=F])
  S_accumulated=as.numeric()
  sp_object<-SpatialPoints(coords)
  candidates<-cbind(coords, 1:dim(sites)[1])
  for(i in 1:n){
    S_accumulated[i]<-vegan::specnumber(pool)
    if(i==n) break() 
    accumulated_plots<-sp_object[included,]
    centroid<-gCentroid(accumulated_plots)@coords
    candidates2<-candidates[-included,, drop=F]
    closest<-candidates2[order(distmat(centroid, candidates2[,-3, drop=F]),  runif(dim(candidates2)[1]) )[1],3,drop=T]
    included=c(included, closest)
    pool<- pool + sites[closest,]
    
  }
  return(S_accumulated)
}

#' Construct spatially constrained sample-based rarefaction (sSBR) curve using the k-Nearest-Centroid-neighbour (k-NCN) algorithm 
#'
#' This function accumulates samples according their proximity to all previously included samples (their centroid) as opposed
#' to the proximity to the inital focal sample. This ensures that included samples mutually close to each other and not all over the place.
#' 
#' Internally the function constructs one curve per sample whereby each sample serves as the initial sample once. Finally, the average curve is returned.
#' @param x a mob_in object or a sites x species matrix
#' @param n number of sites to include. 
#' @param coords spatial coordinates of the samples. If x is a mob_in object, the fuction uses its "spat" table as coords.
#'
#' @return a numeric vector of estimated species richness
#' @export
#'
#' @examples
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
  out<-pbsapply(samples,function(sample)  centroid_accumulate(x = x,focal_sample = sample, coords = coords) )
  return(rowMeans(out))
}


#' Compare all sample-based curves (random, spatiallyconstrained-k-NN, spatially constrained-k-NCN)
#'
#'This is just plotting all curves.
#'
#' @param x a mob_in object
#'
#' @return a plot 
#' @export
#'
#' @examples
compare_curves<-function(x){
  
  centroid_curve<-kNCN_average(x)
  spatcurve<-rarefaction(x,method = "spat")
  samplecurve<-rarefaction(x,method = "samp")
  plot(samplecurve, type= "l", xlab="Samples", ylab="Expected species richness")
  lines(centroid_curve,col=2)
  lines(spatcurve, col=3)
  legend("bottomright", legend = c("sample based", "k-NCN","k-NN"), col=1:3, lty=1)
}

