##------------------------------------------------------------------------------
#' @title sum_formula
#' @description Get the total formula after add adduct.
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param formula Chemical formula.
#' @param adduct Adduct.
#' @export
#' @examples
#' sum_formula(formula = "C9H11NO2", adduct = 'M+H')
#' sum_formula(formula = "C9H11NO2", adduct = 'M+')
#' sum_formula(formula = "C9H11NO2", adduct = 'M+CH3COOH')
#' sum_formula(formula = "C9H11", adduct = 'M-H20')

sum_formula <-
  function(formula = "C9H11NO2",
           adduct = "M-H2O+H") {
    if (is.na(formula)) {
      return(NA)
    }
    
    if (is.na(adduct)) {
      return(formula)
    }
    
    if (adduct == "M+" | adduct == "M-") {
      return(formula)
    }
    
    formula1 <- split_formula(formula)
    adduct1 <-
      strsplit(x = adduct, split = "\\-|\\+")[[1]][-1]
    
    polymer <-
      as.numeric(gsub(
        pattern = "M",
        replacement = "",
        strsplit(x = adduct, split = "\\-|\\+")[[1]][1]
      ))
    
    if (is.na(polymer)) {
      polymer <- 1
    }
    
    plusorminus <- strsplit(x = adduct, split = "")[[1]]
    plusorminus <-
      grep("\\+|\\-", plusorminus, value = TRUE)
    
    formula1$number <- formula1$number * polymer
    
    adduct1 <- mapply(function(x, y) {
      temp <- split_formula(x)
      temp$number <- temp$number * ifelse(y == "+", 1, -1)
      list(temp)
    },
    x = adduct1,
    y = plusorminus)
    
    adduct1 <- do.call(rbind, adduct1)
    
    formula <- rbind(formula1, adduct1)
    rownames(formula) <- NULL
    
    unique.element <- unique(formula$element.name)
    if (length(unique.element) == nrow(formula)) {
      if (any(formula$number < 0)) {
        return(NA)
      } else{
        formula$number[formula$number == 1] <- "W"
        formula <-
          paste(paste(formula$element.name, formula$number, sep = ""),
                collapse = "")
        formula <- strsplit(formula, split = "")[[1]]
        formula[formula == "W"] <- ""
        formula <- paste(formula, collapse = "")
        return(formula)
      }
    } else{
      formula <- lapply(unique.element, function(x) {
        formula[formula$element.name == x, , drop = FALSE]
      })
      
      formula <- lapply(formula, function(x) {
        data.frame(unique(x$element.name),
                   sum(x$number),
                   stringsAsFactors = FALSE)
      })
      
      formula <- do.call(rbind, formula)
      formula <- formula[formula[, 2] != 0, ]
      colnames(formula) <- c("element.name", "number")
      if (any(formula$number < 0)) {
        return(NA)
      } else{
        formula$number[formula$number == 1] <- "W"
        formula <-
          paste(paste(formula$element.name, formula$number, sep = ""),
                collapse = "")
        formula <- strsplit(formula, split = "")[[1]]
        formula[formula == "W"] <- ""
        formula <- paste(formula, collapse = "")
        return(formula)
      }
    }
  }


##------------------------------------------------------------------------------
#' @title split_formula
#' @description Split a formula into element and number.
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param formula Chemical formula.
#' @return A splited formula.
#' @export
#' @examples
#' split_formula(formula = "C9H11NO2")

split_formula <-
  function(formula = "C9H11NO2") {
    temp.formula <- strsplit(formula, split = "")[[1]]
    
    number <- NULL
    for (i in 1:length(temp.formula)) {
      if (length(grep("[0-9]{1}", temp.formula[i])) == 0) {
        break
      }
      number[i] <- temp.formula[i]
    }
    
    if (!is.null(number)) {
      number <- as.numeric(paste(number, collapse = ""))
    } else{
      number <- 1
    }
    ##first select the Na, Cl and so on element
    idx1 <- gregexpr("[A-Z][a-z][0-9]*", formula)[[1]]
    len1 <- attributes(idx1)$match.length
    ##no double element
    if (idx1[1] == -1) {
      double.formula <- matrix(NA, ncol = 2)
      formula1 <- formula
    } else{
      double.letter.element <- NULL
      double.number <- NULL
      remove.idx <- NULL
      for (i in 1:length(idx1)) {
        double.letter.element[i] <-
          substr(formula, idx1[i], idx1[i] + len1[i] - 1)
        if (nchar(double.letter.element[i]) == 2) {
          double.number[i] <- 1
        } else{
          double.number[i] <-
            as.numeric(substr(
              double.letter.element[i],
              3,
              nchar(double.letter.element[i])
            ))
        }
        double.letter.element[i] <-
          substr(double.letter.element[i], 1, 2)
        remove.idx <-
          c(remove.idx, idx1[i]:(idx1[i] + len1[i] - 1))
      }
      
      double.formula <- data.frame(double.letter.element,
                                   double.number, stringsAsFactors = FALSE)
      formula1 <- strsplit(formula, split = "")[[1]]
      formula1 <- formula1[-remove.idx]
      formula1 <- paste(formula1, collapse = "")
    }
    
    ## no one element
    if (formula1 == "") {
      one.formula <- matrix(NA, ncol = 2)
    } else{
      idx2 <- gregexpr("[A-Z][0-9]*", formula1)[[1]]
      len2 <- attributes(idx2)$match.length
      one.letter.element <- NULL
      one.number <- NULL
      for (i in 1:length(idx2)) {
        one.letter.element[i] <-
          substr(formula1, idx2[i], idx2[i] + len2[i] - 1)
        if (nchar(one.letter.element[i]) == 1) {
          one.number[i] <- 1
        } else{
          one.number[i] <-
            as.numeric(substr(one.letter.element[i], 2, nchar(one.letter.element[i])))
        }
        one.letter.element[i] <- substr(one.letter.element[i], 1, 1)
      }
      one.formula <- data.frame(one.letter.element, one.number,
                                stringsAsFactors = FALSE)
    }
    
    colnames(double.formula) <-
      colnames(one.formula) <- c("element.name", "number")
    formula <- rbind(double.formula, one.formula)
    formula <-
      formula[!apply(formula, 1, function(x)
        any(is.na(x))), ]
    
    formula <- formula[order(formula$element.name), ]
    formula$number <- formula$number * number
    unique.element <- unique(formula$element.name)
    if (length(unique.element) == nrow(formula)) {
      return(formula)
    } else{
      formula <- lapply(unique.element, function(x) {
        formula[formula$element.name == x, , drop = FALSE]
      })
      
      formula <- lapply(formula, function(x) {
        data.frame(unique(x$element.name),
                   sum(x$number),
                   stringsAsFactors = FALSE)
      })
      
      formula <- do.call(rbind, formula)
      colnames(formula) <- c("element.name", "number")
      return(formula)
    }
  }
