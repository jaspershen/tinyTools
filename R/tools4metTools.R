#-----------------------------------------------------------------------------
#@title GetDiffMZppm
#@description GetDiffMZppm
#@author Hao Li, Yandong Yin, Kai Weng
#\email{shenxt@@sioc.ac.cn}
#@param mz mz
#@param mz.ppm.thr mz.ppm.thr
#@export
#' @title GetDiffMZppm
#' @export

GetDiffMZppm <- function(mz, mz.ppm.thr = NULL) {
  mz.diff <- diff(mz) / mz[-1] * 1e6
  if (!is.null(mz.ppm.thr)) {
    idx <- which(mz[-1] <= mz.ppm.thr)
    mz.diff[idx] <- mz.diff[idx] * mz[-1] / mz.ppm.thr
  }
  mz.diff
}

# getSpectraMatchScore(pk.spec, x, 
#                      ppm.tol = ppm.ms2match,
#                      mz.ppm.thr = mz.ppm.thr)
####20190530
# exp.spectrum <- exp.spectra
# lib.spectrum <- exp.spectra
###getSpectraMatchScore is used to get two MS2 spectra match score, see MS-DIAL
#' @title getSpectraMatchScore
#' @export
getSpectraMatchScore <- function(exp.spectrum, 
                                 lib.spectrum,
                                 ppm.tol = 30,
                                 mz.ppm.thr = 400){
  exp.spectrum <- as.data.frame(exp.spectrum)
  lib.spectrum <- as.data.frame(lib.spectrum)
  
  exp.spectrum$intensity <- exp.spectrum$intensity/max(exp.spectrum$intensity)  
  lib.spectrum$intensity <- lib.spectrum$intensity/max(lib.spectrum$intensity)  
  
  match.matrix <- ms2Match(exp.spectrum = exp.spectrum, 
                           lib.spectrum = lib.spectrum, 
                           ppm.tol = ppm.tol, 
                           mz.ppm.thr = mz.ppm.thr)
  
  fraction <- sum(!is.na(match.matrix$Lib.index) & !is.na(match.matrix$Exp.index))/nrow(match.matrix)
  
  dp.forward <- getDP(exp.int = match.matrix$Exp.intensity, 
                      lib.int = match.matrix$Lib.intensity)
  dp.reverse <- getDP(exp.int = match.matrix$Exp.intensity[which(match.matrix$Lib.intensity > 0)], 
                      lib.int = match.matrix$Lib.intensity[which(match.matrix$Lib.intensity > 0)])
  dp.forward[is.na(dp.forward)] <- 0
  dp.reverse[is.na(dp.reverse)] <- 0
  score <- dp.forward/3 + dp.reverse/3 + fraction/3
  return(score)
}




#####getDP is used to calculate dot product
#' @title getDP
#' @export
getDP <- function(exp.int, lib.int){
  exp.weight <- sapply(exp.int, function(x){
    1/(1+x/(sum(exp.int) - 0.5))
  })
  lib.weight <- sapply(lib.int, function(x){
    1/(1+x/(sum(lib.int) - 0.5))
  })
  x <- exp.weight * exp.int
  y <- lib.weight * lib.int
  return(sum(x * y) ^ 2 / (sum(x ^ 2) * sum(y ^ 2)))
}


###ms2Match is used to match two MS2 spectra
#' @title ms2Match
#' @export

ms2Match <- function(exp.spectrum, 
                     lib.spectrum,
                     ppm.tol = 30,
                     mz.ppm.thr = 400){
  ##remove noisy fragments
  exp.spectrum <- removeNoise(spec = exp.spectrum, 
                              ppm.ms2match = ppm.tol, 
                              mz.ppm.thr = mz.ppm.thr)
  lib.spectrum <- removeNoise(spec = lib.spectrum, 
                              ppm.ms2match = ppm.tol, 
                              mz.ppm.thr = mz.ppm.thr)
  
  ##for each fragment in lib.spectrum, its matched fragments index in exp.spectrum  
  match.idx <- lapply(lib.spectrum$mz, function(x){
    diff.mz <- abs(x - exp.spectrum$mz)
    x[x < mz.ppm.thr] <- mz.ppm.thr
    mz.error <- diff.mz*10^6/x
    temp.idx <- which(mz.error < ppm.tol)
    if(length(temp.idx) == 0) return(NA)
    if(length(temp.idx) > 1) {
      return(temp.idx[which.max(exp.spectrum$intensity[temp.idx])])
    }
    return(temp.idx)
  })
  
  match.idx <- do.call(rbind, match.idx)
  match.idx <- cbind(1:nrow(match.idx), match.idx)
  colnames(match.idx) <- c("Lib", "Exp")
  
  non.idx2 <- setdiff(c(1:nrow(exp.spectrum)), match.idx[,2][!is.na(match.idx[,2])])
  
  if(length(non.idx2) != 0){
    match.idx2 <- data.frame(NA, non.idx2, stringsAsFactors = FALSE)
    colnames(match.idx2) <- c('Lib', 'Exp')
  }else{
    match.idx2 <- NULL
  }
  
  match.matrix <- as.data.frame(rbind(match.idx, match.idx2), stringsAsFactors = FALSE)
  
  
  match.matrix <- data.frame(match.matrix, 
                             lib.spectrum[match.matrix$Lib,c(1,2)],
                             exp.spectrum[match.matrix$Exp,c(1,2)])
  colnames(match.matrix) <- c("Lib.index", "Exp.index", "Lib.mz", "Lib.intensity",
                              "Exp.mz", "Exp.intensity")
  
  match.matrix$Lib.intensity[is.na(match.matrix$Lib.intensity)] <- 0
  match.matrix$Exp.intensity[is.na(match.matrix$Exp.intensity)] <- 0
  rownames(match.matrix) <- NULL
  match.matrix
}


###for a spectrum, if two fragments with the similar m/z, the fragment with the low fragment should be removed from
##the spectrum
#' @title removeNoise
#' @export
removeNoise <- function(spec, ppm.ms2match = 30, 
                        mz.ppm.thr = 400){
  if(nrow(spec) == 1) return(spec)
  spec <- spec[order(spec[,1]),]
  mz <- spec[,1]
  mz <- mz[-1]
  diff.mz <- diff(spec[,1])
  mz[which(mz < mz.ppm.thr)] <- mz.ppm.thr
  mz.error <- diff.mz*10^6/mz
  temp.idx <- which(mz.error < ppm.ms2match)
  if(length(temp.idx) > 0){
    remove.idx <- lapply(temp.idx, function(idx){
      c(idx, idx + 1)[which.min(spec[c(idx, idx + 1), 2])]
    })
    
    remove.idx <- unique(unlist(remove.idx))
    spec <- spec[-remove.idx, , drop = FALSE]
  }else{
    return(spec)
  }
}
