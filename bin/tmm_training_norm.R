#Necessary functions for performing TMM normalization on a DGEList using training samples only as the reference samples

.calcFactorRLE <- function(data)
  #	Scale factors as in Anders et al (2010)
  #	Mark Robinson
  #	Created 16 Aug 2010
{
  gm <- exp(rowMeans(log(data)))
  apply(data, 2, function(u) median((u/gm)[gm > 0]))
}

.calcFactorQuantile <- function (data, lib.size, p=0.75)
  #	Generalized version of upper-quartile normalization
  #	Mark Robinson
  #	Created 16 Aug 2010. Last modified 2 Jun 2020.
{
  y <- t(t(data)/lib.size)
  f <- apply(y,2,function(x) quantile(x,probs=p))
}

.calcFactorTMM <- function(obs, ref, libsize.obs=NULL, libsize.ref=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10)
  #	TMM between two libraries
  #	Mark Robinson
{
  obs <- as.numeric(obs)
  ref <- as.numeric(ref)
  
  if( is.null(libsize.obs) ) nO <- sum(obs) else nO <- libsize.obs
  if( is.null(libsize.ref) ) nR <- sum(ref) else nR <- libsize.ref
  
  logR <- log2((obs/nO)/(ref/nR))          # log ratio of expression, accounting for library size
  absE <- (log2(obs/nO) + log2(ref/nR))/2  # absolute expression
  v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref   # estimated asymptotic variance
  
  #	remove infinite values, cutoff based on A
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
  
  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]
  
  if(max(abs(logR)) < 1e-6) return(1)
  
  #	taken from the original mean() function
  n <- length(logR)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS
  
  #	keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
  #	a fix from leonardo ivan almonacid cardenas, since rank() can return
  #	non-integer values when there are a lot of ties
  keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)
  
  if(doWeighting)
    f <- sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE)
  else
    f <- mean(logR[keep], na.rm=TRUE)
  
  #	Results will be missing if the two libraries share no features with positive counts
  #	In this case, return unity
  if(is.na(f)) f <- 0
  2^f
}

.calcFactorTMMwsp <- function(obs, ref, libsize.obs=NULL, libsize.ref=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10)
  #	TMM with pairing of singleton positive counts between the obs and ref libraries
  #	Gordon Smyth
  #	Created 19 Sep 2018. Last modified 23 April 2019.
{
  obs <- as.numeric(obs)
  ref <- as.numeric(ref)
  
  #	epsilon serves as floating-point zero
  eps <- 1e-14
  
  #	Identify zero counts
  pos.obs <- (obs > eps)
  pos.ref <- (ref > eps)
  npos <- 2L * pos.obs + pos.ref
  
  #	Remove double zeros and NAs
  i <- which(npos==0L | is.na(npos))
  if(length(i)) {
    obs <- obs[-i]
    ref <- ref[-i]
    npos <- npos[-i]
  }
  
  #	Check library sizes	
  if(is.null(libsize.obs)) libsize.obs <- sum(obs)
  if(is.null(libsize.ref)) libsize.ref <- sum(ref)
  
  #	Pair up as many singleton positives as possible
  #	The unpaired singleton positives are discarded so that no zeros remain
  zero.obs <- (npos == 1L)
  zero.ref <- (npos == 2L)
  k <- (zero.obs | zero.ref)
  n.eligible.singles <- min( sum(zero.obs), sum(zero.ref))
  if(n.eligible.singles > 0L) {
    refk <- sort(ref[k],decreasing=TRUE)[1:n.eligible.singles]
    obsk <- sort(obs[k],decreasing=TRUE)[1:n.eligible.singles]
    obs <- c(obs[!k],obsk)
    ref <- c(ref[!k],refk)
  } else {
    obs <- obs[!k]
    ref <- ref[!k]
  }
  
  #	Any left?
  n <- length(obs)
  if(n==0L) return(1)
  
  #	Compute M and A values
  obs.p <- obs / libsize.obs
  ref.p <- ref / libsize.ref
  M <- log2( obs.p / ref.p )
  A <- 0.5 * log2( obs.p * ref.p )
  
  #	If M all zero, return 1
  if(max(abs(M)) < 1e-6) return(1)
  
  #	M order, breaking ties by shrunk M
  obs.p.shrunk <- (obs+0.5) / (libsize.obs+0.5)
  ref.p.shrunk <- (ref+0.5) / (libsize.ref+0.5)
  M.shrunk <- log2( obs.p.shrunk / ref.p.shrunk )
  o.M <- order(M, M.shrunk)
  
  #	A order
  o.A <- order(A)
  
  #	Trim
  loM <- as.integer(n * logratioTrim) + 1L
  hiM <- n + 1L - loM
  keep.M <- rep.int(FALSE,n)
  keep.M[o.M[loM:hiM]] <- TRUE
  loA <- as.integer(n * sumTrim) + 1L
  hiA <- n + 1L - loA
  keep.A <- rep.int(FALSE,n)
  keep.A[o.A[loA:hiA]] <- TRUE
  keep <- keep.M & keep.A
  M <- M[keep]
  
  #	Average the M values
  if(doWeighting) {
    obs.p <- obs.p[keep]
    ref.p <- ref.p[keep]
    v <- (1-obs.p)/obs.p/libsize.obs + (1-ref.p)/ref.p/libsize.ref
    w <- (1+1e-6) / (v+1e-6)
    TMM <- sum(w*M) / sum(w)
  } else {
    TMM <- mean(M)
  }
  
  2^TMM
}

#Code to retrieve the reference sample for TMM

getTMMref <- function(object, samples.for.training = NULL){
  
  #	Check object
  x <- as.matrix(object)
  if(any(is.na(x))) stop("NA counts not permitted")
  nsamples <- ncol(x)
  
  #	Check lib.size
  lib.size <- colSums(x)
  
  #	Remove all zero rows
  allzero <- .rowSums(x>0, nrow(x), nsamples) == 0L
  if(any(allzero)) x <- x[!allzero,,drop=FALSE]
  
  #	Degenerate cases
  if(nrow(x)==0 || nsamples==1){
    warning("Only one sample or gene")
    refColumn <- 1L
  } 
  
  # Ensure that refgroup samples come from training if training samples provided
  if (!is.null(samples.for.training)) {
    refGroupSamples <- rownames(object$samples)[rownames(object$samples) %in% samples.for.training]
  } else {
    refGroupSamples <- rownames(object$samples)
  }
  
  #return(refGroupSamples)
  
  # MB: select lib.size values of the selected samples
  lib.size.refGroup <- object$samples$lib.size[which(rownames(object$samples) %in% refGroupSamples)]
  #return(lib.size.refGroup)
  
  #First part of TMM method, with sample filter if training only
  f75 <- .calcFactorQuantile(data=x[,refGroupSamples],
                             lib.size=lib.size.refGroup,
                             p=0.75)
  
  #Set refColumn, filter for training samples in line above
  refColumn <- which.min(abs(f75-mean(f75)))
  if(length(refColumn)==0L | refColumn < 1 | refColumn > nsamples) refColumn <- 1L
  
  data.frame(ref.sample = names(refColumn),
             col.index = refColumn[[1]])
  
}

# https://rdrr.io/bioc/edgeR/src/R/calcNormFactors.R
# custom calcNormFactorsThromboseq that includes TMM-reference sample selection from only
# Similar/identical to calcNormFactorsThromboseq from ThromboSeq package,
# but avoids some situational production of NA values
# Updated for most recent version of edgeR code (July 2020)
calcNormFactorsTraining <- function(object, lib.size=NULL,
                                    #method=c("TMM","TMMwsp","RLE","upperquartile","none"), #Method is always TMM for our purposes
                                    method = "TMM",
                                    refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE,
                                    Acutoff=-1e10, p=0.75,
                                    normalize.on.training.series = TRUE, samples.for.training = NULL, #New params
                                    ...){
  #	Check object
  x <- as.matrix(object)
  if(any(is.na(x))) stop("NA counts not permitted")
  nsamples <- ncol(x)
  
  #	Check lib.size
  if(is.null(lib.size)) {
    lib.size <- colSums(x)
  } else {
    if(anyNA(lib.size)) stop("NA lib.sizes not permitted")
    if(length(lib.size) != nsamples) {
      if(length(lib.size) > 1L) warning("calcNormFactors: length(lib.size) doesn't match number of samples",call.=FALSE)
      lib.size <- rep_len(lib.size,nsamples)
    }
  }
  
  
  #	Remove all zero rows
  allzero <- .rowSums(x>0, nrow(x), nsamples) == 0L
  if(any(allzero)) x <- x[!allzero,,drop=FALSE]
  
  #	Degenerate cases
  if(nrow(x)==0 || nsamples==1) method="none"
  
  # Check parameter congruity
  if (normalize.on.training.series == TRUE & is.null(samples.for.training)){
    stop("normalize.on.training.series is true, but samples.for.training was not provided")
  }
  
  # Ensure that refgroup samples come from training if training samples provided
  if (!is.null(samples.for.training)) {
    refGroupSamples <- rownames(object$samples)[rownames(object$samples) %in% samples.for.training]
  } else {
    refGroupSamples <- rownames(object$samples)
  }
  
  #return(refGroupSamples)
  stopifnot(length(refGroupSamples) == length(samples.for.training))
  
  # lib.size from filtered training samples
  lib.size.refGroup <- object$samples$lib.size[which(rownames(object$samples) %in% refGroupSamples)]
  #return(lib.size.refGroup)
  stopifnot(length(lib.size.refGroup) == length(samples.for.training))
  
  #	Calculate factors
  # filter for refGroupSamples
  f75 <- .calcFactorQuantile(data=x[,refGroupSamples],
                             lib.size=lib.size.refGroup,
                             p=0.75)
  #return(f75)
  #Set refColumn, filter for training samples in line above
  if( is.null(refColumn) ) refColumn <- which.min(abs(f75-mean(f75)))
  if(length(refColumn)==0L | refColumn < 1 | refColumn > nsamples) refColumn <- 1L
  
  f <- rep_len(NA_real_,nsamples)
  for(i in 1:nsamples)
    f[i] <- .calcFactorTMM(obs=x[,i],
                           ref=x[,refColumn],
                           libsize.obs=lib.size[i],
                           libsize.ref=lib.size[refColumn],
                           logratioTrim=logratioTrim,
                           sumTrim=sumTrim,
                           doWeighting=doWeighting,
                           Acutoff=Acutoff)
  
  #return(f) 
  
  #f <- f/exp(mean(log(f)))
  if (normalize.on.training.series == TRUE) {
    f <- f/exp(mean(log(f[which(colnames(object) %in% 
                                  samples.for.training)])))
  }
  else {
    f <- f/exp(mean(log(f)))
  }
  
  #	Output
  names(f) <- colnames(x)
  #return(f) #Add it directly to dge instead
  object$samples$norm.factors <- f
  object
  
  
}

#Wrapper function for calculating norm factors and then performing cpm and log2 transformation
normalize.tmm <- function(
  dge,
  training.samples,
  tmm.ref.from.training = T,
  refCol = dge$TMMref$col.index
){
  
  dge <- calcNormFactorsTraining(object = dge,
                                 normalize.on.training.series = tmm.ref.from.training,
                                 samples.for.training = training.samples,
                                 refColumn=refCol)
  
  # calculate counts-per-million matrix, log-transformed and normalized via the TMM-normalization factor
  # (cpm with normalized.lib.size=T after calculating norm factors yields a TMM-normalized matrix)
  cpm(dge, log = T, normalized.lib.sizes = T) 
  
}
