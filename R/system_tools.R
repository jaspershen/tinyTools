##------------------------------------------------------------------------------
#' @title get_os
#' @description Convert from wide data to long data.
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @export
#' @examples 
#' get_os()

get_os <- function() {
  sysinf <- Sys.info()
  if (!is.null(sysinf)) {
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else {
    ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}