#' @title setwd_project
#' @description Set work directory where Rproj object in.
#' @author Xiaotao Shen
#' \email{shenxt@@163.com}
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom stringr str_split
#' @importFrom utils head
#' @importFrom utils tail
#' @export

setwd_project <- function() {
  currect_wd <-
    getwd()
  
  candidate_wd <-
    currect_wd %>%
    stringr::str_split("/") %>%
    unlist()
  
  if (length(candidate_wd) == 1) {
    candidate_wd <- currect_wd
  } else{
    candidate_wd <-
      lapply(2:length(candidate_wd), function(i) {
        paste(candidate_wd[1:i], collapse = "/")
      })
  }
  
  candidate_wd <-
    rev(candidate_wd)
  
  for (i in 1:length(candidate_wd)) {
    wd <- candidate_wd[[i]]
    file_name <-
      list.files(wd,
                 recursive = ifelse(wd == currect_wd, TRUE, FALSE),
                 full.names = TRUE)
    project_index <-
      grep(".Rproj", file_name)
    
    if (length(project_index) != 0) {
      project_wd <-
        file_name[project_index[1]] %>%
        stringr::str_split("/") %>%
        unlist() %>%
        head(-1) %>%
        paste(collapse = "/")
      cat(
        "The project name is:",
        file_name[project_index[1]] %>%
          stringr::str_split("/") %>%
          unlist() %>%
          tail(1),
        "\n"
      )
      cat("The project wd is:",
          project_wd, "\n")
      
      setwd(project_wd)
      break()
    }
  }
  if (length(project_index) == 0) {
    cat("There are no .Rproj in your file. No change for wd.\n")
  }
}