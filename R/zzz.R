.onAttach <- function(libname, pkgname) {
  needed <- core[!is_attached(core)]
  if (length(needed) == 0)
    return()
  
  crayon::num_colors(TRUE)
  tinyTools_attach()
  
  
  packageStartupMessage(
    crayon::green(
      "metPath,
More information can be found at https://jaspershen.github.io/tinyTools/
Authors: Xiaotao Shen (shenxt@stanford.edu)
Maintainer: Xiaotao Shen"
    )
  )
}

is_attached <- function(x) {
  paste0("package:", x) %in% search()
}



globalVariables(names = c(
  "Exp.intensity",
  "Exp.mz",
  "Lib.intensity",
  "Lib.mz",
  "intensity",
  "mz"
))
