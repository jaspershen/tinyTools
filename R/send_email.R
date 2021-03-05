#' ##------------------------------------------------------------------------------
#' #' @title convert2long
#' #' @description Convert from wide data to long data.
#' #' @author Xiaotao Shen
#' #' \email{shenxt@@stanford.edu}
#' #' @param expression_data expression_data
#' #' @param sample_info sample_info
#' #' @param variable_info variable_info
#' #' @importFrom magrittr %>%
#' #' @importFrom tibble rownames_to_column
#' #' @importFrom tidyr pivot_longer
#' #' @importFrom dplyr left_join
#' #' @importFrom dplyr select
#' #' @importFrom dplyr everything
#' #' @export
#' 
#' 
#' 
#' setGeneric(
#'   name = "sendMail",
#'   def = function(email.address = "zutaoshen@gmail.com",
#'                  from = "kagawabale@163.com",
#'                  user.name = "kagawabale",
#'                  pw = "19900521@2016") {
#'     to <- email.address
#'     
#'     subject <- "metID noticement"
#'     
#'     body <- "<html> <b>Dear user:</b>
#'          <br/>
#'          <br/>Thank you for using MetDNA.
#'          <br/>The analysis of your data is done and the results can be found in attachment.
#'          <br/>The attachment contains a analysis report (HTML) and a zip file for analysis result.
#'          <br/>Please don't be hesitate to send email to us (shenxt@sioc.ac.cn) if you have any questions.
#'          <br/>
#'          <br/> <b>Your sincerely,</b>
#'          <br/> <b>XiaotaoShen and Zhengjiaing Zhu<b>
#'          <br/><b>Website:</b> <a href = \"http://www.zhulab.cn/\"> www.zhulab.cn </a>
#'          <br/>
#'          <br/>
#'          <img src = \"http://a3.qpic.cn/psb?/V12nMOGs2VfuNZ/L0OOyQhz2rQKFWLs0NXqS2UyWhXP8g21XcJXoZDAE80!/b/dD0BAAAAAAAA&bo=sQd2AgAAAAADAOc!&rf=viewer_4\" heigh = '600', width = '600'>
#'          </html>"
#'     
#'     smtp <- list(
#'       host.name = "smtp.163.com",
#'       port = 465,
#'       user.name = user.name,
#'       passwd = pw,
#'       ssl = TRUE
#'     )
#' 
#'       mailR::send.mail(
#'         from = from,
#'         to = to,
#'         subject = subject,
#'         body = body,
#'         smtp = smtp,
#'         authenticate = TRUE,
#'         send = TRUE,
#'         html = TRUE
#'       )
#'     
#'     
#'     
#'   }
#' )