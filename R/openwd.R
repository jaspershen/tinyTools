#' @title openwd
#' @description Open current work directory.
#' @author Xiaotao Shen
#' \email{shenxt@@163.com}
#' @return inder.
#' @export

openwd <- function(){
  system(sprintf('open %s', shQuote(getwd())))
}
