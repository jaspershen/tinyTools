#' @title ms2_plot2
#' @description Plot ms2 spectra or ms2 match plot.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param spectrum1 Spectrum 1.
#' @param spectrum2 Spectrum 2.
#' @param spectrum1_name Name of spectrum1
#' @param spectrum2_name Name of spectrum2
#' @param range.mz range.mz
#' @param ppm.tol ppm.tol
#' @param mz.ppm.thr mz.ppm.thr
#' @param xlab xlab.
#' @param ylab ylab.
#' @param col1 Color 1.
#' @param col2 Color 2.
#' @param title.size title.size
#' @param lab.size lab.size
#' @param axis.text.size axis.text.size.
#' @param legend.title.size legend.title.size
#' @param legend.text.size legend.text.size
#' @param interactive_plot Interactive plot or not.
#' @return Return a MS2 spectrum.
#' @export
#' @examples 
#' spectrum1 = data.frame(
#' mz = c(
#'   87.50874,
#'   94.85532,
#'   97.17808,
#'   97.25629,
#'   103.36186,
#'   106.96647,
#'   107.21461,
#'   111.00887,
#'   113.79269,
#'   118.70564
#' ),
#' intensity =
#'   c(
#'     8356.306,
#'     7654.128,
#'     9456.207,
#'     8837.188,
#'     8560.228,
#'     8746.359,
#'     8379.361,
#'     169741.797,
#'     7953.080,
#'     8378.066
#'   )
#' )
#' spectrum2 = spectrum1
#' ms2_plot(spectrum1, spectrum2)
#' ms2_plot(spectrum1, spectrum2, interactive_plot = TRUE)
#' ms2_plot(spectrum1)
#' ms2_plot(spectrum1, interactive_plot = TRUE)
#' 
#' 

ms2_plot = function(spectrum1,
                    spectrum2,
                    spectrum1_name = "spectrum1",
                    spectrum2_name = "spectrum2",
                    range.mz,
                    ppm.tol = 30,
                    mz.ppm.thr = 400,
                    xlab = "Mass to charge ratio (m/z)",
                    ylab = "Relative intensity",
                    col1 = "red",
                    col2 = "black",
                    title.size = 15,
                    lab.size = 15,
                    axis.text.size = 15,
                    legend.title.size = 15,
                    legend.text.size = 15,
                    interactive_plot = FALSE) {
  if (missing(spectrum1) & missing(spectrum2)) {
    stop("No spectrum1 and spectrum2")
  }
  if (!missing(spectrum1)) {
    spectrum1 =
      apply(spectrum1, 2, function(x) {
        as.numeric(x)
      }) %>%
      as.data.frame()
    
    spectrum1[, 2] <- spectrum1[, 2] / max(spectrum1[, 2])
  }
  
  if (!missing(spectrum2)) {
    spectrum2 =
      apply(spectrum2, 2, function(x) {
        as.numeric(x)
      }) %>%
      as.data.frame()
    spectrum2[, 2] <- spectrum2[, 2] / max(spectrum2[, 2])
  }
  
  ##two spectrum
  if (!missing(spectrum1) & !missing(spectrum2)) {
    if (missing(range.mz)) {
      range.mz <- c(min(spectrum1[, 1], spectrum2[, 1]),
                    max(spectrum1[, 1], spectrum2[, 1]))
      
    }
    
    matched.spec <- ms2Match(spectrum1,
                             spectrum2,
                             ppm.tol = ppm.tol,
                             mz.ppm.thr = mz.ppm.thr)
    
    matched.idx <- which(matched.spec[, "Lib.intensity"] > 0 &
                           matched.spec[, "Exp.intensity"] > 0)
    
    plot <- ggplot(matched.spec) +
      geom_segment(
        mapping = aes(
          x = Exp.mz,
          y = Exp.intensity - Exp.intensity,
          xend = Exp.mz,
          yend = Exp.intensity
        ),
        colour = col2
      ) +
      geom_point(
        data = matched.spec[matched.idx, , drop = FALSE],
        mapping = aes(x = Exp.mz, y = Exp.intensity),
        colour = col2
      ) +
      xlim(range.mz[1], range.mz[2]) +
      labs(x = xlab, y = ylab) +
      scale_y_continuous(
        limits = c(-1, 1),
        breaks = c(-1, -0.5, 0, 0.5, 1),
        labels = c("1", "0.5", "0", "0.5", "1")
      ) +
      theme_bw() +
      theme(
        # axis.line = element_line(arrow = arrow()),
        plot.title = element_text(
          color = "black",
          size = title.size,
          face = "plain",
          hjust = 0.5
        ),
        axis.title = element_text(
          color = "black",
          size = lab.size,
          face = "plain"
        ),
        axis.text = element_text(
          color = "black",
          size = axis.text.size,
          face = "plain"
        ),
        legend.title = element_text(
          color = "black",
          size = legend.title.size,
          face = "plain"
        ),
        legend.text = element_text(
          color = "black",
          size = legend.text.size,
          face = "plain"
        )
      )
    
    plot <- plot +
      annotate(
        geom = "text",
        x = Inf,
        y = Inf,
        label = spectrum1_name,
        color = col2,
        hjust = 1,
        vjust = 1
      ) +
      annotate(
        geom = "text",
        x = Inf,
        y = -Inf,
        label = spectrum2_name,
        color = col1,
        hjust = 1,
        vjust = -1
      )
    
    plot <- plot +
      geom_segment(
        data = matched.spec,
        mapping = aes(
          x = Lib.mz,
          y = Lib.intensity - Lib.intensity,
          xend = Lib.mz,
          yend = -Lib.intensity
        ),
        colour = col1
      ) +
      geom_point(
        data = matched.spec[matched.idx, , drop = FALSE],
        mapping = aes(x = Lib.mz, y = -Lib.intensity),
        colour = col1
      )
    if (interactive_plot) {
      plot = plotly::ggplotly(plot)
    }
    return(plot)
  }
  
  if (!missing(spectrum1) & missing(spectrum2)) {
    spectrum = spectrum1
  }
  
  if (!missing(spectrum2) & missing(spectrum1)) {
    spectrum = spectrum2
  }
  
  if (missing(range.mz)) {
    range.mz <- c(min(spectrum[, 1]), max(spectrum[, 1]))
  }
  
  plot <- ggplot(spectrum) +
    geom_segment(mapping = aes(
      x = mz,
      y = 0,
      xend = mz,
      yend = intensity
    ),
    colour = col1) +
    xlim(range.mz[1], range.mz[2]) +
    labs(x = xlab, y = ylab) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_bw() +
    theme(
      # axis.line = element_line(arrow = arrow()),
      plot.title = element_text(
        color = "black",
        size = title.size,
        face = "plain",
        hjust = 0.5
      ),
      axis.title = element_text(
        color = "black",
        size = lab.size,
        face = "plain"
      ),
      axis.text = element_text(
        color = "black",
        size = axis.text.size,
        face = "plain"
      ),
      legend.title = element_text(
        color = "black",
        size = legend.title.size,
        face = "plain"
      ),
      legend.text = element_text(
        color = "black",
        size = legend.text.size,
        face = "plain"
      )
    )
  
  if (interactive_plot) {
    plot = plotly::ggplotly(plot)
  }
  plot
}
