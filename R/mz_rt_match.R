#' @title mz_rt_match
#' @description Match peaks according to m/z and RT.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param data1 First data for matching, first column must be mz
#' and seconod column must be rt.
#' @param data2 Second data for matching, first column must be mz
#' and seconod column must be rt.
#' @param mz.tol mz tol for ms1 and ms2 data matching.
#' @param rt.tol RT tol for ms1 and ms2 data matching.
#' @param rt.error.type RT error is calculated with relative or absolute.
#' @return Return a result which give the matching result of data1 and database.
#' @export

mz_rt_match =
  function(data1,
           data2,
           mz.tol,
           #rt.tol is relative
           rt.tol = 30,
           rt.error.type = c("relative", "abs")) {
    rt.error.type <- match.arg(rt.error.type)
    #
    if (nrow(data1) == 0 | nrow(data2) == 0) {
      result <- NULL
      return(result)
    }
    # mz1 <- as.numeric(data1[, 1])
    # rt1 <- as.numeric(data1[, 2])
    info1 <- data1[, c(1, 2)]
    info1 <- apply(info1, 1, list)
    
    mz2 <- as.numeric(data2[, 1])
    rt2 <- as.numeric(data2[, 2])
    
    result <- pbapply::pblapply(info1, function(x) {
      temp.mz1 <- x[[1]][[1]]
      temp.rt1 <- x[[1]][[2]]
      mz.error <- abs(temp.mz1 - mz2) * 10 ^ 6 / temp.mz1
      if (rt.error.type == "relative") {
        rt.error <- abs(temp.rt1 - rt2) * 100 / temp.rt1
      } else{
        rt.error <- abs(temp.rt1 - rt2)
      }
      
      j <- which(mz.error <= mz.tol & rt.error <= rt.tol)
      if (length(j) == 0) {
        matrix(NA, ncol = 7)
      } else{
        cbind(j, temp.mz1, mz2[j], mz.error[j], temp.rt1, rt2[j], rt.error[j])
      }
    })
    
    if (length(result) == 1) {
      result <- cbind(1, result[[1]])
    } else{
      result <- mapply(function(x, y) {
        list(cbind(x, y))
      },
      x <- 1:length(info1),
      y = result)
      result <- do.call(rbind, result)
    }
    
    result <-
      matrix(result[which(!apply(result, 1, function(x)
        any(is.na(x)))), ], ncol = 8)
    if (nrow(result) == 0)
      return(NULL)
    colnames(result) <-
      c("Index1",
        "Index2",
        "mz1",
        "mz2",
        "mz error",
        "rt1",
        "rt2",
        "rt error")
    result <- as.data.frame(result)
    result
  }


#' @title keep_one
#' @description Remove multiple vs. one in result from mz_rt_match function.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param result result from mz_rt_match function.
#' @param according.to According to mz error or rt error?
#' @return Return a result without multiple vs. one.
#' @export

keep_one = function(result,
                    according.to = c("mz.error", "rt.error")) {
  according.to <- match.arg(according.to)
  if (is.null(result))
    return(result)
  if (class(result) != "matrix" &
      class(result) != "data.frame")
    stop("result must be matrix or data.frame.")
  if (ncol(result) != 8)
    stop("result must from mz_rt_match.")
  if (paste(colnames(result), collapse = ";") != "index1;index2;mz1;mz2;mz.error;rt1;rt2;rt.error") {
    stop("result must from mz_rt_match.")
  }
  
  duplicated.idx <- unique(result$index1[duplicated(result$index1)])
  if (length(duplicated.idx) == 0)
    return(result)
  
  for (i in 1:length(duplicated.idx)) {
    temp.idx <- which(result$index1 == duplicated.idx[i])
    temp.result <- result[temp.idx, ]
    if (according.to == "mz.error") {
      temp.idx1 <- temp.idx[which.min(temp.result$mz.error)]
      temp.idx2 <- setdiff(temp.idx, temp.idx1)
      result <- result[-temp.idx2, ]
    }
    
    if (according.to == "rt.error") {
      temp.idx1 <- temp.idx[which.min(temp.result$rt.error)]
      temp.idx2 <- setdiff(temp.idx, temp.idx1)
      result <- result[-temp.idx2, ]
    }
    
  }
  result
}
