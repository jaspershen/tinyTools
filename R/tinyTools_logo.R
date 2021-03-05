#' @title tinyTools_logo
#' @description Get the detailed of metPath package.
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @importFrom stringr str_replace str_split str_replace_all str_trim
#' @importFrom dplyr filter mutate select everything
#' @importFrom rstudioapi isAvailable hasFun getThemeInfo
#' @importFrom utils packageDescription
#' @importFrom cli rule symbol
#' @importFrom crayon green blue col_align red black white style make_style num_colors
#' @importFrom plotly ggplotly
#' @importFrom pbapply pblapply
#' @import ggplot2
#' @importFrom methods .hasSlot new
#' @importFrom stats p.adjust
#' @importFrom utils data str
#' @importFrom magrittr %>%
#' @export
#' @examples
#' tinyTools_logo()

tinyTools_logo <- function(){
  cat(crayon::green("Thank you for using metPath!\n"))
  cat(crayon::green("Version 0.9.0 (20210303)\n"))
  cat(crayon::green("More information can be found at https://jaspershen.github.io/tinyTools/\n"))
  cat(crayon::green(
    c("  _   _          _______          _     ", " | | (_)        |__   __|        | |    ", 
      " | |_ _ _ __  _   _| | ___   ___ | |___ ", " | __| | '_ \\| | | | |/ _ \\ / _ \\| / __|", 
      " | |_| | | | | |_| | | (_) | (_) | \\__ \\", "  \\__|_|_| |_|\\__, |_|\\___/ \\___/|_|___/", 
      "               __/ |                    ", "              |___/                     "
    )
    
  ), sep = "\n")
}



# library(cowsay)
# # https://onlineasciitools.com/convert-text-to-ascii-art
# # writeLines(capture.output(say("Hello"), type = "message"), con = "ascii_art.txt")
# art <- readLines("logo.txt")
# dput(art)
# metid_logo <-
#   c("                _    _____  ___ ", " _ __ ___   ___| |_  \\_   \\/   \\",
#     "| '_ ` _ \\ / _ \\ __|  / /\\/ /\\ /", "| | | | | |  __/ |_/\\/ /_/ /_// ",
#     "|_| |_| |_|\\___|\\__\\____/___,'  ", "                                "
#   )
# cat(metid_logo, sep = "\n")