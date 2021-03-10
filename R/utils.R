###getSpectraMatchScore is used to get two MS2 spectra match score, see MS-DIAL
#' @title getSpectraMatchScore
#' @param exp.spectrum exp.spectrum
#' @param lib.spectrum lib.spectrum
#' @param ppm.tol ppm.tol
#' @param mz.ppm.thr mz.ppm.thr
#' @param fraction.weight fraction.weight
#' @param dp.forward.weight dp.forward.weight
#' @param dp.reverse.weight dp.reverse.weight
#' @export
getSpectraMatchScore <- function(exp.spectrum,
                                 lib.spectrum,
                                 ppm.tol = 30,
                                 mz.ppm.thr = 400,
                                 fraction.weight = 0.2,
                                 dp.forward.weight = 0.7,
                                 dp.reverse.weight = 0.1) {
  exp.spectrum <- as.data.frame(exp.spectrum)
  lib.spectrum <- as.data.frame(lib.spectrum)
  
  exp.spectrum$intensity <-
    exp.spectrum$intensity / max(exp.spectrum$intensity)
  lib.spectrum$intensity <-
    lib.spectrum$intensity / max(lib.spectrum$intensity)
  
  match.matrix <- ms2Match(
    exp.spectrum = exp.spectrum,
    lib.spectrum = lib.spectrum,
    ppm.tol = ppm.tol,
    mz.ppm.thr = mz.ppm.thr
  )
  
  fraction <-
    sum(!is.na(match.matrix$Lib.index) &
          !is.na(match.matrix$Exp.index)) / nrow(match.matrix)
  
  dp.forward <- getDP(exp.int = match.matrix$Exp.intensity,
                      lib.int = match.matrix$Lib.intensity)
  dp.reverse <-
    getDP(exp.int = match.matrix$Exp.intensity[which(match.matrix$Lib.intensity > 0)],
          lib.int = match.matrix$Lib.intensity[which(match.matrix$Lib.intensity > 0)])
  dp.forward[is.na(dp.forward)] <- 0
  dp.reverse[is.na(dp.reverse)] <- 0
  score <-
    dp.forward * dp.forward.weight + dp.reverse * dp.reverse.weight + fraction * fraction.weight
  return(score)
}




#####getDP is used to calculate dot product
#' @title getDP
#' @param exp.int exp.int
#' @param lib.int lib.int
#' @export
getDP <- function(exp.int, lib.int) {
  exp.weight <- sapply(exp.int, function(x) {
    1 / (1 + x / (sum(exp.int) - 0.5))
  })
  lib.weight <- sapply(lib.int, function(x) {
    1 / (1 + x / (sum(lib.int) - 0.5))
  })
  x <- exp.weight * exp.int
  y <- lib.weight * lib.int
  return(sum(x * y) ^ 2 / (sum(x ^ 2) * sum(y ^ 2)))
}


###ms2Match is used to match two MS2 spectra
#' @title ms2Match
#' @param exp.spectrum exp.spectrum
#' @param lib.spectrum lib.spectrum
#' @param ppm.tol ppm.tol
#' @param mz.ppm.thr mz.ppm.thr
#' @export

ms2Match <- function(exp.spectrum,
                     lib.spectrum,
                     ppm.tol = 30,
                     mz.ppm.thr = 400) {
  ##remove noisy fragments
  exp.spectrum <- removeNoise(spec = exp.spectrum,
                              ppm.ms2match = ppm.tol,
                              mz.ppm.thr = mz.ppm.thr)
  lib.spectrum <- removeNoise(spec = lib.spectrum,
                              ppm.ms2match = ppm.tol,
                              mz.ppm.thr = mz.ppm.thr)
  
  ##for each fragment in lib.spectrum, its matched fragments index in exp.spectrum
  match.idx <- lapply(lib.spectrum$mz, function(x) {
    diff.mz <- abs(x - exp.spectrum$mz)
    x[x < mz.ppm.thr] <- mz.ppm.thr
    mz.error <- diff.mz * 10 ^ 6 / x
    temp.idx <- which(mz.error < ppm.tol)
    if (length(temp.idx) == 0)
      return(NA)
    if (length(temp.idx) > 1) {
      return(temp.idx[which.max(exp.spectrum$intensity[temp.idx])])
    }
    return(temp.idx)
  })
  
  match.idx <- do.call(rbind, match.idx)
  match.idx <- cbind(1:nrow(match.idx), match.idx)
  colnames(match.idx) <- c("Lib", "Exp")
  
  non.idx2 <-
    setdiff(c(1:nrow(exp.spectrum)), match.idx[, 2][!is.na(match.idx[, 2])])
  
  if (length(non.idx2) != 0) {
    match.idx2 <- data.frame(NA, non.idx2, stringsAsFactors = FALSE)
    colnames(match.idx2) <- c('Lib', 'Exp')
  } else{
    match.idx2 <- NULL
  }
  
  match.matrix <-
    as.data.frame(rbind(match.idx, match.idx2), stringsAsFactors = FALSE)
  
  
  match.matrix <- data.frame(match.matrix,
                             lib.spectrum[match.matrix$Lib, c(1, 2)],
                             exp.spectrum[match.matrix$Exp, c(1, 2)])
  colnames(match.matrix) <-
    c("Lib.index",
      "Exp.index",
      "Lib.mz",
      "Lib.intensity",
      "Exp.mz",
      "Exp.intensity")
  
  match.matrix$Lib.intensity[is.na(match.matrix$Lib.intensity)] <- 0
  match.matrix$Exp.intensity[is.na(match.matrix$Exp.intensity)] <- 0
  rownames(match.matrix) <- NULL
  match.matrix
}


###for a spectrum, if two fragments with the similar m/z, the fragment with the low fragment should be removed from
##the spectrum
#' @title removeNoise
#' @param spec spec
#' @param ppm.ms2match ppm.ms2match
#' @param mz.ppm.thr mz.ppm.thr
#' @export
removeNoise <- function(spec,
                        ppm.ms2match = 30,
                        mz.ppm.thr = 400) {
  if (nrow(spec) == 1)
    return(spec)
  spec <- spec[order(spec[, 1]),]
  mz <- spec[, 1]
  mz <- mz[-1]
  diff.mz <- diff(spec[, 1])
  mz[which(mz < mz.ppm.thr)] <- mz.ppm.thr
  mz.error <- diff.mz * 10 ^ 6 / mz
  temp.idx <- which(mz.error < ppm.ms2match)
  if (length(temp.idx) > 0) {
    remove.idx <- lapply(temp.idx, function(idx) {
      c(idx, idx + 1)[which.min(spec[c(idx, idx + 1), 2])]
    })
    
    remove.idx <- unique(unlist(remove.idx))
    spec <- spec[-remove.idx, , drop = FALSE]
  } else{
    return(spec)
  }
}


# # NOTE: keep in sync with install-patRoon version
# findPWizPath <- function() {
#   # try to find ProteoWizard
#   # order: options --> win registry --> PATH
#   # the PATH is searched last because OpenMS might have added its own old version.
#   
#   path <- getOption("patRoon.path.pwiz")
#   if (!is.null(path) && nzchar(path))
#     return(path)
#   
#   if (Sys.info()[["sysname"]] == "Windows")
#   {
#     # Inspired by scan_registry_for_rtools() from pkgload
#     key <- "Directory\\shell\\Open with SeeMS\\command"
#     reg <-
#       tryCatch(
#         utils::readRegistry(key, "HCR"),
#         error = function(e)
#           NULL
#       )
#     
#     # not sure if this might occur
#     if (is.null(reg))
#       reg <-
#       tryCatch(
#         utils::readRegistry(key, "HLM"),
#         error = function(e)
#           NULL
#       )
#     
#     if (!is.null(reg))
#     {
#       path <-
#         tryCatch(
#           dirname(sub("\"([^\"]*)\".*", "\\1", reg[[1]])),
#           error = function(e)
#             NULL
#         )
#       if (!is.null(path) &&
#           file.exists(file.path(path, "msconvert.exe")))
#         # extra check: see if msconvert is there
#         return(path)
#     }
#   }
#   
#   # check PATH
#   msc <-
#     if (Sys.info()[["sysname"]] == "Windows")
#       "msconvert.exe"
#   else
#     "msconvert"
#   path <- dirname(Sys.which(msc))
#   if (nzchar(path))
#     return(path)
#   
#   return(NULL)
# }
