#-----------------------------------------------------------------------------
setGeneric(name = "removeNoise",
           def = function(spec, mz.tol = 30){
             spec <- matrix(spec, ncol = 2)
             colnames(spec) <- c("mz", "intensity")
             if(nrow(spec) == 1) return(spec)
             spec <- spec[order(spec[,1]),]
             mz <- as.numeric(spec[,1])
             
             new.spec <- NULL
             i = 1
             while(i < length(mz)){
               # cat(i); cat(" ")
               temp.mz <- mz[i]
               mz.error <- abs(temp.mz - mz)*10^6/ifelse(temp.mz >= 400, temp.mz, 400)
               temp.idx <- which(mz.error <= mz.tol)
               temp.spec <- spec[temp.idx,]
               if(length(temp.idx) == 1) {
                 new.spec <- rbind(new.spec, temp.spec)
                 i <- max(temp.idx) + 1
                 next()
               }
               temp.mz <- median(temp.spec[,1])
               temp.int <- max(temp.spec[,2])
               new.spec <- rbind(new.spec, c(temp.mz, temp.int))
               # spec <- spec[-temp.idx,]
               i <- max(temp.idx) + 1
             }
             
             row.names(new.spec) <- NULL
             colnames(new.spec) <- c("mz", "intensity")
             return(new.spec)
           })



#-----------------------------------------------------------------------------
###
###GetMatchResult
#@title GetMatchResult
#@description GetMatchResult
#@author Hao Li, Yandong Yin, Kai Weng
#\email{shenxt@@sioc.ac.cn}
#@param spec.exp spec.exp
#@param spec.lib spec.lib
#@param weight.int weight.int
#@param weight.mz weight.mz
#@param ppm.ms2match ppm.ms2match
#@param mz.ppm.thr mz.ppm.thr
#@param ppm.sanity.check ppm.sanity.check
#@param is.sanity.check is.sanity.check
#@param direction direction
#@param ... other parameters
#@export
#'

GetMatchResult <- function(spec.exp, spec.lib,
                           weight.int = 1,
                           weight.mz = 0,
                           ppm.ms2match = 30,
                           mz.ppm.thr = 400,
                           ppm.sanity.check = 100,
                           is.sanity.check = FALSE,
                           direction = direction,
                           ...) {
  # ()
  if (is.sanity.check) {
    switch(direction,
           'forward' = {
             if (any(c(GetDiffMZppm(spec.exp[, 'mz']), GetDiffMZppm(spec.lib[, 'mz'])) <= ppm.sanity.check)) {
               stop('Difference between m/z is too small!!')
             }
           },
           'reverse' = {
             if (any(GetDiffMZppm(spec.lib) <= ppm.sanity.check)) {
               stop('Difference between m/z is too small!!')
             }
           },
           stop('Error setup for parameter: direction!!!'))
  }
  #
  spec2match <- GetSpec2Match(spec.exp, spec.lib,
                              ppm.ms2match = ppm.ms2match,
                              mz.ppm.thr = mz.ppm.thr,
                              direction = direction)
  int.weighted.pk  <- GetWeightedInt(spec2match$exp, weight.mz, weight.int)
  int.weighted.lib <- GetWeightedInt(spec2match$lib, weight.mz, weight.int)
  match.score <- GetDotProduct(int.weighted.pk, int.weighted.lib)
  attr(match.score, 'spec') <- spec2match
  attr(match.score, 'spec.compared') <- cbind('exp' = int.weighted.pk,
                                              'lib' = int.weighted.lib)
  return(match.score)
}



#-----------------------------------------------------------------------------
###GetSpec2Match
#@title GetSpec2Match
#@description GetSpec2Match
#@author Hao Li, Yandong Yin, Kai Weng
#\email{shenxt@@sioc.ac.cn}
#@param spec.exp spec.exp
#@param spec.lib spec.lib
#@param ppm.ms2match ppm.ms2match
#@param mz.ppm.thr mz.ppm.thr
#@param direction direction
#@export
#'

GetSpec2Match <- function(spec.exp, spec.lib,
                          ppm.ms2match = 30,
                          mz.ppm.thr = 400,
                          direction = c('reverse', 'forward')) {
  #
  direction = match.arg(direction)
  mz.pool   <- sort(c(spec.exp[, 'mz'], spec.lib[, 'mz']))
  spec.temp <- cbind('mz' = mz.pool, 'intensity' = 0)
  
  spec.exp.temp  <- MatchFromTemp(spec.exp, spec.temp)
  spec.lib.temp <- MatchFromTemp(spec.lib, spec.temp)
  
  # combine nearby peaks
  pk.spec  <- MatchSpec(spec.exp.temp,  ppm.ms2match = ppm.ms2match, mz.ppm.thr = mz.ppm.thr)
  lib.spec <- MatchSpec(spec.lib.temp, ppm.ms2match = ppm.ms2match, mz.ppm.thr = mz.ppm.thr)
  
  if (direction == 'reverse') {
    idx.rm <- which(lib.spec[, 'intensity'] == 0)
    if (length(idx.rm > 0)) {
      pk.spec  <- pk.spec[-idx.rm, , drop = FALSE]
      lib.spec <- lib.spec[-idx.rm, , drop = FALSE]
    }
  }
  
  return(list('exp' = pk.spec, 'lib' = lib.spec))
}



#-----------------------------------------------------------------------------
###MatchFromTemp
#@title MatchFromTemp
#@description MatchFromTemp
#@author Hao Li, Yandong Yin, Kai Weng
#\email{shenxt@@sioc.ac.cn}
#@param spec spec
#@param temp temp
#@export
#'

MatchFromTemp <- function(spec, temp) {
  temp[match(spec[, 'mz'], temp[, 'mz']), 'intensity'] <- spec[, 'intensity']
  temp
}




#-----------------------------------------------------------------------------
###MatchSpec
#@title MatchSpec
#@description MatchSpec
#@author Hao Li, Yandong Yin, Kai Weng
#\email{shenxt@@sioc.ac.cn}
#@param spec spec
#@param ppm.ms2match ppm.ms2match
#@param mz.ppm.thr mz.ppm.thr
#@export
#'


MatchSpec <- function(spec, ppm.ms2match = 30, mz.ppm.thr = 400) {
  while (TRUE) {
    mz.diff.ppm <- GetDiffMZppm(spec[, 'mz'], mz.ppm.thr = mz.ppm.thr)
    idx <- which(mz.diff.ppm < ppm.ms2match)
    if (length(idx) > 0) {
      i <- tail(idx, 1)
      j <- which.max(spec[c(i, i + 1), 'intensity'])
      spec[i, 'intensity'] <- spec[i + j - 1, 'intensity']
      i2 <- i + 1
      spec[i, 'mz'] <- spec[i2, 'mz']
      spec <- spec[-i - 1, , drop = FALSE]
    } else {
      break
    }
  }
  return(spec)
}



#-----------------------------------------------------------------------------
#@title GetDiffMZppm
#@description GetDiffMZppm
#@author Hao Li, Yandong Yin, Kai Weng
#\email{shenxt@@sioc.ac.cn}
#@param mz mz
#@param mz.ppm.thr mz.ppm.thr
#@export

GetDiffMZppm <- function(mz, mz.ppm.thr = NULL) {
  mz.diff <- diff(mz) / mz[-1] * 1e6
  if (!is.null(mz.ppm.thr)) {
    idx <- which(mz[-1] <= mz.ppm.thr)
    mz.diff[idx] <- mz.diff[idx] * mz[-1] / mz.ppm.thr
  }
  mz.diff
}



#-----------------------------------------------------------------------------
#@title GetWeightedInt
#@description GetWeightedInt
#@author Hao Li, Yandong Yin, Kai Weng
#\email{shenxt@@sioc.ac.cn}
#@param spec spec
#@param weight.mz weight.mz
#@param weight.int weight.int
#@export
#'


GetWeightedInt <- function(spec, weight.mz = 0, weight.int = 1) {
  return(spec[, 'mz'] ^ weight.mz * spec[, 'intensity'] ^ weight.int)
}


#-----------------------------------------------------------------------------
#'@title GetDotProduct
#'@description GetDotProduct
#'@author Hao Li, Yandong Yin, Kai Weng
#'\email{shenxt@@sioc.ac.cn}
#'@param x x
#'@param y y
#'@export

GetDotProduct <- function(x, y) {
  # return(sum(x * y) ^ 2 / (sum(x ^ 2) * sum(y ^ 2)))
  return(sum(x * y) / sqrt((sum(x ^ 2) * sum(y ^ 2))))
}


#-----------------------------------------------------------------------------

TuneMS2 <- function(spec,
                    mz.precursor,
                    is.include.precursor = TRUE,
                    snthr = 3,
                    noise.ms2 = 3,
                    mz.range.ms2 = NULL,
                    int.ms2.min.abs,
                    int.ms2.min.relative = 0.01,
                    is.apply.ms2.min.relative = TRUE,
                    is.check.sanity = TRUE,
                    int.check.sanity = 50,
                    ppm.precursor.filter = 20, ...) {
  # re-orgnize spec with increasing mz
  spec <- spec[order(spec[, 'mz']), , drop = FALSE]
  # define the lowest non-noise signal
  if (missing(int.ms2.min.abs)) {
    int.ms2.min.abs <- noise.ms2 * snthr
  }
  #
  # provent noise over estimation
  if (is.check.sanity & int.ms2.min.abs > int.check.sanity) {
    stop('Noise estimation is too high!')
  }
  
  # considering precursor ion
  if (missing(mz.precursor)) {
    mz.precursor <- max(spec[, 'mz'])
  }
  mz.precursor.range <- GetRangePPM(mz.precursor, ppm.precursor.filter)
  idx.mz.precursor.range <- ifelse(is.include.precursor, 2, 1)
  mz.cutoff <- mz.precursor.range[idx.mz.precursor.range]
  spec <- spec[spec[,'mz'] < mz.cutoff, , drop = FALSE]
  
  if (nrow(spec) == 0) {
    return()
  }
  
  if (!is.null(mz.range.ms2)) {
    nr.keep <- which(spec[, 'mz'] >= mz.range.ms2[1] &
                       spec[, 'mz'] <= mz.range.ms2[2])
    if (length(nr.keep) > 0) {
      spec <- spec[nr.keep, , drop = FALSE]
    }
    else {
      return()
    }
  }
  
  if ((max(spec[, 'intensity']) * int.ms2.min.relative) < int.ms2.min.abs) {
    is.warning.lowspec <- TRUE
  }
  
  # discarding low intensity spec (1% highest int and int.ms2.min.abs)
  int.cutoff <- max(max(spec[, 'intensity']) * int.ms2.min.relative,
                    int.ms2.min.abs)
  spec <- spec[spec[, 'intensity'] >= int.cutoff, , drop = FALSE]
  if (nrow(spec) == 0) {
    return()
  }
  
  # discarding ring effects
  spec <- RemoveRingEffect(spec)
}


#-----------------------------------------------------------------------------

RemoveRingEffect <- function(spec, mz.diff.thr = 0.3, int.rel.thr = 0.2) {
  nr.ring <- nrow(spec) + 1
  mz <- spec[, 'mz']
  
  mz.diff <- diff(mz)
  idx.mzdiff <- which(mz.diff <= mz.diff.thr)
  if (length(idx.mzdiff) == 0) {
    return(spec)
  }
  
  nr.ring.possible <- unique(c(idx.mzdiff, idx.mzdiff + 1))
  while (TRUE) {
    
    idx.int.max <- which.max(spec[nr.ring.possible, 2])
    nr.int.max <- nr.ring.possible[idx.int.max]
    int.thr <- spec[nr.int.max, 2] * int.rel.thr
    
    mz.diff <- abs(mz[nr.ring.possible[-idx.int.max]] - mz[nr.int.max])
    int <- spec[nr.ring.possible[-idx.int.max], 2]
    nr.ring <- append(nr.ring, nr.ring.possible[-idx.int.max][which(mz.diff <= mz.diff.thr & int <= int.thr)])
    nr.ring.possible <- nr.ring.possible[!nr.ring.possible %in% c(nr.ring, nr.int.max)]
    if (length(nr.ring.possible) == 0) {
      break
    }
  }
  
  return(spec[-nr.ring, , drop = FALSE])
}


#-----------------------------------------------------------------------------

NormalizeSpec <- function(spec, ref = max(abs(spec[, 2])), pos = 'top') {
  spec[, 2] <- spec[, 2] / ref * ifelse(pos == 'down', -1, 1)
  return(spec)
}



#-----------------------------------------------------------------------------

GetRangePPM <- function(data, ppm) {
  t(sapply(data, function(dt) dt * (1 + c(-1, 1) * ppm * 1e-6)))
}
#-----------------------------------------------------------------------------
