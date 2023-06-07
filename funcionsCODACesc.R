#



Encerts=function(A,tt){
if(diff(dim(A))!=0){stop("A no és quadrada")}
n=dim(A)[1]
if (n==1){optim=0} 
else {
perms=gtools::permutations(n,n)
valors=c(sum(diag(A)),sum((A-diag(tt))^2))
for (j in 2:dim(perms)[1]){
x=perms[j,]  
valors=rbind(valors,c(sum(diag(A[,x])),sum((A[,x]-diag(tt))^2)))
optim=min(valors[valors[1,]==max(valors[1,]),2])
}
}
  return(optim)
}

#
#' Adaptada de EasyCODA
#' 
codaSeq.outlier=function(y, plot.me=TRUE, col=rgb(1,0,0,0.3)){
  pcx=prcomp(y)
  mv=sum(pcx$sdev^2)
  sample.var= apply(pcx$x,1,function(y){sum(y^2/mv)})
  cut=median(apply(pcx$x,1,function(x){sum(x^2/mv)})) + 2 * IQR(apply(pcx$x,1,function(x){sum(x^2/mv)}))
  bad=which(apply(pcx$x,1,function(x){sum(x^2/mv)}) > cut) ##
  good=which(apply(pcx$x,1,function(x){sum(x^2/mv)}) <= cut) ##
  if(plot.me == TRUE){
    hist(sample.var, breaks=100)
    boxplot(sample.var, horizontal=TRUE, col=col, add=TRUE)
    abline(v=cut, lty=2)
  }
  return(list(sample.var=sample.var, bad=bad, good=good) )
}
#
#' Adaptada de CODA4microbiome
#' Fa servir la AUC de predicció multinomial de 
#' https://link.springer.com/article/10.1023/A:1010920819831   
#'           
ELRmulti.HT=function (x, y){
  k <- ncol(x)
  lrmatrix <- coda4microbiome::logratios_matrix(x)
  lrX <- lrmatrix[[1]]
  idlrX <- lrmatrix[[2]]
  logratio_cor <- matrix(rep(0, k * k), nrow = k)
  s = 0
  for (i in (1:(k - 1))) {
    for (j in ((i + 1):k)) {
      s <- s + 1
      lr_ij <- lrX[, s]
      model1=nnet::multinom(y ~ lr_ij, trace = F)
      result1= predict(model1, lr_ij, type='probs')
      logratio_cor[i, j]=
        HandTill2001::auc(multcap(response = y,
                                  predicted = as.matrix(result1)))
      logratio_cor[j, i] <- logratio_cor[i, j]
    }
  }
  o <- order(colSums(abs(logratio_cor)), decreasing = T)
  M <- logratio_cor[o, o]
  colnames(M) <- o
  rownames(M) <- colnames(M)
  #
  results <- list(
    `max log-ratio` = colnames(M)[which(M == max(abs(M)), arr.ind = TRUE)[(2:1)]],
    `names max log-ratio` = colnames(x)[as.numeric(colnames(M)[which(M ==max(abs(M)), arr.ind = TRUE)[(2:1)]])], 
    `order of importance` = o, 
    `most important variables` = colnames(x)[o], 
    `association matrix` = M
  )
  return(results)
}
#
#' Separat el gràfic del càlcul de l'anterior
#' X el resultat de ELRmulti.HT()
Gràfic.CP=function(X,shownames = FALSE, maxrow = 15, maxcol = 15){
  title <- "AUC multinomial regression y~log(xi/xj)"
  GnBu8 <- c("#f7fcf0", "#e0f3db", "#ccebc5", "#a8ddb5", "#7bccc4", 
                      "#4eb3d3", "#2b8cbe", "#08589e")
Blues4 <- c("#eff3ff", "#bdd7e7", "#6baed6", "#2171b5")
col2 <- grDevices::colorRampPalette(GnBu8, space = "Lab")
k <- ncol(X$`association matrix`)
if (maxrow > k) maxrow <- k
if (maxcol > k) maxcol <- k
o=X$`order of importance`
M=X$`association matrix`
if (shownames == TRUE) {rownames(M) <- colnames(x)[o]}
corrplot::corrplot(M[(1:maxrow), (1:maxcol)],
                    tl.pos = "lt", tl.col="black", title = title, mar =c(0, 0, 1, 0), 
                    col.lim = c(min(M), max(M)), col = col2(200), 
                    is.corr = FALSE)
}
#
#'



dendo.més.barplot=function(X,cols){
X.CLR=easyCODA::CLR(X)
hc=easyCODA::WARD(X.CLR,weight=TRUE)
dend=as.dendrogram(hc)
labels_colors(dend)=cols[hc$order]
XOr=X[hc$order,]
#' Reorden les mostres per dibuixar els barplot en el mateix ordre
XOr.CLR=acomp(XOr)
d.names=colnames(X)[order(apply(X, 2, sum), decreasing=T) ]
nb.cols=dim(X)[2]
colors.OTU=colorRampPalette(brewer.pal(length(d.names),"Spectral"))(nb.cols)
#' Resultat
results=list(`hc`=hc,`den`=dend)
return(results)
#' Dibuix
layout(matrix(c(1,3,2,3),2,2, byrow=T), widths=c(6,2), height=c(4,4))
par(mar=c(2,1,1,1)+0.1,cex=0.75)
plot(dend, main = "")
barplot(XOr.CLR, legend.text=F, col=colors.OTU, axisnames=F, border=NA, xpd=T,)
par(mar=c(0,1,1,1)+0.1,cex=1)
plot(1,2, pch = 1, lty = 1, ylim=c(-20,20), type = "n", axes = FALSE, ann = FALSE)
legend(x="center", legend=d.names, col=colors.OTU, lwd=5, cex=.6, border=NULL)
}

#
#' Adaptacions d'ALDEx
#' 

#'rdirichlet usual, però no pot donar 0
#'

rdirichlet.nn <- function (n, alpha)
{
  n <- as.integer(n)
  if(is.vector(alpha)) alpha <- t(alpha)
  l <- dim(alpha)[2]
  x <- matrix(rgamma(l * n, t(alpha)), ncol = l, byrow=TRUE)  
  x[which(x==0)]=10^(-100)
  return(x / rowSums(x))
}

#' aldexCesc.clr.function d'ALDEx, però (a) empra la rdirichlet anterior;
#' (b) s'aplica a matriu amb zeros ja imputats
#'

aldexCesc.clr <- function( reads, conds, mc.samples=128, denom="all", verbose=FALSE, useMC=FALSE, summarizedExperiment=NULL ) {
  
  # INPUT
  # The 'reads' data.frame MUST have row
  # and column names that are unique, and
  # looks like the following:
  #
  #              T1a T1b  T2  T3  N1  N2
  #   Gene_00001   0   0   2   0   0   1
  #   Gene_00002  20   8  12   5  19  26
  #   Gene_00003   3   0   2   0   0   0
  #       ... many more rows ...
  #
  # ---------------------------------------------------------------------
  
  # OUTPUT
  # The output returned is a list (x) that contains Monte-Carlo instances of
  # the centre log-ratio transformed values for each sample
  # Access to values
  # sample IDs: names(x)
  # number of features (genes, OTUs): length(x[[1]][,1])
  # number of Monte-Carlo Dirichlet instances: length(x[[1]][1,])
  # feature names: rownames(x[[1]])
  
  # coerce SummarizedExperiment reads into data.frame
  if (summarizedExperiment) {
    reads <- data.frame(as.list(assays(reads,withDimnames=TRUE)))
    if (verbose) {
      print("converted SummarizedExperiment read count object into data frame")
    }
  }
  
  # Fully validate and coerce the data into required formats
  # make sure that the multicore package is in scope and return if available
  has.BiocParallel <- FALSE
  if ("BiocParallel" %in% rownames(installed.packages()) & useMC){
    print("multicore environment is is OK -- using the BiocParallel package")
    #require(BiocParallel)
    has.BiocParallel <- TRUE
  }
  else {
    print("operating in serial mode")
  }
  
  # make sure that mc.samples is an integer, despite it being a numeric type value
  as.numeric(as.integer(mc.samples))
  
  #  remove all rows with reads less than the minimum set by minsum
  minsum <- 0
  
  # remove any row in which the sum of the row is 0
  z <- as.numeric(apply(reads, 1, sum))
  reads <- as.data.frame( reads[(which(z > minsum)),]  )
  
  if (verbose) print("removed rows with sums equal to zero")
  
  
  #  SANITY CHECKS ON THE DATA INPUT
  # if ( any( round(reads) != reads ) ) stop("not all reads are integers")
  if ( any( reads < 0 ) )             stop("one or more reads are negative")
  
  for ( col in names(reads) ) {
    if ( any( ! is.finite( reads[[col]] ) ) )  stop("one or more reads are not finite")
  }
  
  if ( length(rownames(reads)) == 0 ) stop("rownames(reads) cannot be empty")
  if ( length(colnames(reads)) == 0 ) stop("colnames(reads) cannot be empty")
  
  if ( length(rownames(reads)) != length(unique(rownames(reads))) ) stop ("row names are not unique")
  if ( length(colnames(reads)) != length(unique(colnames(reads))) ) stop ("col names are not unique")
  if ( mc.samples < 128 ) warning("values are unreliable when estimated with so few MC smps")
  
  # add a prior expection to all remaining reads that are 0
  # this should be by a Count Zero Multiplicative approach, but in practice
  # this is not necessary because of the large number of features
  prior <- 0
  
  # This extracts the set of features to be used in the geometric mean computation
  feature.subset <- aldex.set.mode(reads, conds, denom)
  
  
  reads <- reads + prior
  
  if (verbose == TRUE) print("data format is OK")
  
  # ---------------------------------------------------------------------
  # Generate a Monte Carlo instance of the frequencies of each sample via the Dirichlet distribution,
  # returns frequencies for each feature in each sample that are consistent with the
  # feature count observed as a proportion of the total counts per sample given
  # technical variation (i.e. proportions consistent with error observed when resequencing the same library)
  
  nr <- nrow( reads )
  rn <- rownames( reads )
  
  #this returns a list of proportions that are consistent with the number of reads per feature and the
  #total number of reads per sample
  
  # environment test, runs in multicore if possible
  if (has.BiocParallel){
    p <- bplapply( reads ,
                   function(col) {
                     q <- t( rdirichlet.nn( mc.samples, col ) ) ;
                     rownames(q) <- rn ;
                     q })
    names(p) <- names(reads)
  }
  else{
    p <- lapply( reads ,
                 function(col) {
                   q <- t( rdirichlet.nn( mc.samples, col ) ) ;
                   q[q==0]=10^(-100);
                   rownames(q) <- rn ; q } )
  }
  
  # sanity check on the data, should never fail
  for ( i in 1:length(p) ) {
    if ( any( ! is.finite( p[[i]] ) ) ) stop("non-finite frequencies estimated")
  }
  
  if (verbose == TRUE) print("dirichlet samples complete")
  
  # ---------------------------------------------------------------------
  # Take the log2 of the frequency and subtract the geometric mean log2 frequency per sample
  # i.e., do a centered logratio transformation as per Aitchison
  
  # apply the function over elements in a list, that contains an array
  
  # DEFAULT
  if(length(feature.subset) == nr)
  {
    # Default ALDEx2
    if (has.BiocParallel){
      l2p <- bplapply( p, function(m) {
        apply( log2(m), 2, function(col) { col - mean(col) } )
      })
      names(l2p) <- names(p)
    }
    else{
      l2p <- lapply( p, function(m) {
        apply( log2(m), 2, function(col) { col - mean(col) } )
      })
    }
  } else {
    ## IQLR or ZERO
    feat.result <- vector("list", length(unique(conds))) # Feature Gmeans
    condition.list <- vector("list", length(unique(conds)))    # list to store conditions
    
    for (i in 1:length(unique(conds)))
    {
      condition.list[[i]] <- which(conds == unique(conds)[i]) # Condition list
      feat.result[[i]] <- lapply( p[condition.list[[i]]], function(m) {
        apply(log2(m), 2, function(x){mean(x[feature.subset[[i]]])})
      })
    }
    set.rev <- unlist(feat.result, recursive=FALSE) # Unlist once to aggregate samples
    p.copy <- p
    for (i in 1:length(set.rev))
    {
      p.copy[[i]] <- as.data.frame(p.copy[[i]])
      p[[i]] <- apply(log2(p.copy[[i]]),1, function(x){ x - (set.rev[[i]])})
      p[[i]] <- t(p[[i]])
    }
    l2p <- p    # Save the set in order to generate the aldexCesc.clr variable
  }
  
  
  # sanity check on data
  for ( i in 1:length(l2p) ) {
    if ( any( ! is.finite( l2p[[i]] ) ) ) stop("non-finite log-frequencies were unexpectedly computed")
  }
  if (verbose == TRUE) print("clr transformation complete")
  
  return(new("aldex.clr",reads=reads,conds=conds,mc.samples=mc.samples,verbose=verbose,useMC=useMC,analysisData=l2p))
}
