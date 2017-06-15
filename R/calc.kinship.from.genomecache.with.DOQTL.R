#' @export
calc.kinship.from.genomecache.with.DOQTL <- function(genomecache, 
                                                     model=c("additive", "full", "dominant")){
  model <- model[1]
  h <- DiploprobReader$new(genomecache)
  mapping <- straineff.mapping.matrix()
  
  n <- length(h$getSubjects())
  h <- length(h$getFounders())
  p <- length(h$getLoci())
  if(model == "additive"){
    probs <- array(NA, dim=c(n, h, p))
  }
  else if(model == "full"){
    probs <- array(NA, dim=c(n, h + choose(h, 2), p))
  }
  else if(model == "dominant"){
    probs <- array(NA, dim=c(n, choose(h, 2), p))
  }
  
  loci <- h$getLoci()
  
  for(i in 1:length(loci)){
    diplotypes <- h$getLocusMatrix(loci[i], model="full")
    if(any(diplotypes < 0)){
      diplotypes[diplotypes < 0] <- 0
      diplotypes <- t(apply(diplotypes, 1, function(x) x/sum(x)))
    }
    if(model == "additive"){
      use.probs <- (diplotypes %*% mapping)/2
    }
    else if(model == "full"){
      use.probs <- diplotypes
    }
    else if(model == "dominant"){
      use.probs <- diplotypes[,-(1:8)]
    }
    probs[,,i] <- use.probs
  }
  
  K <- DOQTL::kinship.probs(probs)
  colnames(K) <- rownames(K) <- h$getSubjects()
  return(K)
}