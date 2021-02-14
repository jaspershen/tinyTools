##------------------------------------------------------------------------------
#' @title convert2long
#' @description Convert from wide data to long data.
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param expression_data expression_data
#' @param sample_info sample_info
#' @param variable_info variable_info
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join
#' @importFrom dplyr select
#' @importFrom dplyr everything
#' @export

convert2long <-
  function(expression_data,
           sample_info,
           variable_info) {
    final_data <-
      expression_data %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "variable_id") %>%
      tidyr::pivot_longer(
        cols = -variable_id,
        names_to = "sample_id",
        values_to = "value"
      ) %>%
      dplyr::left_join(sample_info, by = "sample_id") %>%
      dplyr::left_join(variable_info, by = "variable_id") %>%
      dplyr::select(sample_id, variable_id, value, dplyr::everything())
    return(final_data)
  }


##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
#' @title convert2wide
#' @description Convert from wide data to long data.
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param data data
#' @param sample_info_name sample_info_name
#' @param variable_info_name variable_info_name
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr distinct
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join
#' @importFrom dplyr select
#' @importFrom dplyr everything
#' @export

# sample_info_name = c(
#   "sample_id",
#   "subject_id",
#   "accurate_time",
#   "day",
#   "time",
#   "time_window",
#   "seconds",
#   "hour"
# )
#
# variable_info_name = c("variable_id", "mol_name", "class")

convert2wide <-
  function(data,
           sample_info_name = c("sample_id"),
           variable_info_name = c("variable_id")) {
    variable_info <-
      data %>%
      dplyr::select(variable_info_name) %>%
      dplyr::distinct(.keep_all = TRUE) %>%
      as.data.frame() %>%
      dplyr::select(variable_id, dplyr::everything())
    
    expression_data <-
      data %>%
      dplyr::select(sample_id, variable_id, value) %>%
      tidyr::pivot_wider(names_from = sample_id,
                         values_from = value) %>%
      as.data.frame() %>%
      tibble::column_to_rownames(var = "variable_id")
    
    sample_info <-
      data %>%
      dplyr::select(sample_info_name) %>%
      dplyr::distinct(.keep_all = TRUE) %>%
      as.data.frame() %>%
      dplyr::select(sample_id, dplyr::everything())
    
    expression_data =
      expression_data[, sample_info$sample_id]
    expression_data =
      expression_data[variable_info$variable_id, ]
    
    
    final_data <-
      list(
        sample_info = sample_info,
        variable_info = variable_info,
        expression_data = expression_data
      )
    return(final_data)
  }