###########################
# specificty distribution #
###########################
#' Density plot of transformed specificity
#'
#' Plot density curves of 9 methods. Transformation methods in this function,
#' refers to paper, please!
#'
#' @param lst List of 3 data.frame. one of them contains psi values with their
#' specificty values of 9 methods: "Tau", "Gini", "Tsi", "Counts", "Ee", "Pem",
#' "Hg", "Z", "Spm"(named "raw"), the second contains rank values and binary
#' index(named "rank"), the third contains binary pattern values and
#' index binary "Ib"(named "bin"), which were generated from \code{ts_psi} or
#' \code{ts_expr}.
#' @param ymax Numeric. maximum of y limitation. Default 12.
#' @importFrom gplots rich.colors
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom stats density
#' @importFrom stats quantile
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices palette
#' @export
#' @examples
#' \dontrun{
#' result <- ts_psi(demo_psi,
#'                  tissues = c("sample_A", "sample_B", "sample_C",
#'                              "sample_D", "sample_E", "sample_F",
#'                              "sample_G", "sample_H", "sample_I",
#'                              "sample_J", "sample_K", "sample_L",
#'                              "sample_M", "sample_N", "sample_O",
#'                              "sample_P", "sample_Q"),
#'                              identifier = "AS_events")
#' plot_density(result)
#' }
plot_density <- function(lst,
                         ymax = 12) {
  #1:Backgroundcolor for all graphs, 2: Foregroundcolor for all graphs (E6E6E6),
  #3: Fill for histograms, 4: Red, for boxplots, 5: Blue, for boxplots, 6: Green,
  #for boxplots, 7: Light gray
  my.col <- colorRampPalette(c("#FFFFFF",
                               "black",
                               "blue",
                               "#FA8072",
                               "#00A2FF",
                               "#00CC00",
                               "#E0E0E0"))(7)
  if (!is.numeric(ymax) | (ymax <= 0)) {
    stop("ymax must be a numeric and larger than 0!")
  }
  opar<-par(no.readonly = TRUE)
  # pdf(file= filename, height=9, width=12)
  par(cex.main = 0.95,
      bg = my.col[1],
      fg = my.col[2],
      col.axis = my.col[2],
      col.lab = my.col[2],
      col.main = my.col[2])
  palette(rev(rich.colors(10)))
  #palette(rev(blues9))

  df <- lst$raw
  plot(density(df[, "Tau"], n = 1000),
       main = " ",
       xlab = "Tissue specificity",
       col = (1), lwd = 4, lty = 1,
       ylim = c(0, ymax),
       xlim = c(-0.1, 1.1)
  )
  lines(density(df[, "Gini"], n = 1000), col = (2), lwd = 4, lty = 2)
  lines(density(df[, "Tsi"], n = 1000), col = (3), lwd = 4, lty = 1)
  lines(density(df[, "Counts"], n = 1000), col = (4), lwd = 4, lty = 2)
  lines(density(df[, "Ee"], n = 1000), col = (5), lwd = 4, lty = 1)
  lines(density(df[, "Hg"], n = 1000), col = (6), lwd = 4, lty = 2)
  lines(density(df[, "Zscore"], n = 1000), col = (7), lwd = 4, lty = 1)
  lines(density(df[, "Spm"], n = 1000), col = (8), lwd = 4, lty = 2)
  lines(density(df[, "Pem"], n = 1000), col = (9), lwd = 4, lty = 1)

  legend("topright",
         c("Tau", "Gini", "TSI", "Counts", "EE", "Hg", "Zscore", "SPM", "PEM"),
         col = (1:11),
         lwd = 4, lty = c(1,2),
         bty = "n",
         seg.len =4)

  par(opar)
  # dev.off()
}



###########################
#         heatpam         #
###########################
#' Heatmap for subset specific-expression
#'
#' Wrapper function derived from \code{pheatmap}. Use the results of
#' \code{ts_expr} or \code{ts_psi} as Input. Using pheatmap youself may be
#' better for yourwork sometimes.
#'
#' @param lst List of 3 data.frame. one of them contains psi values with their
#' specificty values of 9 methods: "Tau", "Gini", "Tsi", "Counts", "Ee", "Pem",
#' "Hg", "Z", "Spm"(named "raw"), the second contains rank values and binary
#' index(named "rank"), the third contains binary pattern values and
#' index binary "Ib"(named "bin"), which were generated from \code{ts_psi} or
#' \code{ts_expr}.
#' @param dat_type "raw", "rank" or c('raw', 'rank'). Plot "raw"(psi/FPKM/RPKM/TPM) or
#' "rank"(0-n+1). Default "raw". "raw" work with \code{specificity}, and
#' "rank" work with \code{Ib}. If specify c('raw', 'rank'), then plot the heatmap of
#' raw value(psi or tpm) of Ib.
#' @param specificity Vector of range, numeric. Region of specificity to plot,
#' work with \code{dat_type}="raw". Default \code{c(0.8, 1)}.
#' @param ts_method Specificity methods. "Tau", "Gini", "TSI", "Counts", "EE",
#' "Hg", "Zscore", "Spm" or "Pem".
#' @param Ib Integer. Cutoff of binary index, work with \code{dat_type}="rank".
#' @param angle_col refer to the same parameter in \code{pheatmap}.
#' @param cluster_cols refer to the same parameter in \code{pheatmap}.
#' @param cluster_rows refer to the same parameter in \code{pheatmap}.
#' @param show_rownames refer to the same parameter in \code{pheatmap}.
#' @param ... parameters of \code{pheatmap} pkg.
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices palette
#' @export
#' @examples
#' result <- ts_psi(demo_psi,
#'                  tissues = c("sample_A", "sample_B", "sample_C",
#'                              "sample_D", "sample_E", "sample_F",
#'                              "sample_G", "sample_H", "sample_I",
#'                              "sample_J", "sample_K", "sample_L",
#'                              "sample_M", "sample_N", "sample_O",
#'                              "sample_P", "sample_Q"),
#'                              identifier = "AS_events")
#' plot_heatmap(result)
#' plot_heatmap(result, dat_type = "rank", Ib = 2)
plot_heatmap <- function(lst,
                         dat_type = "raw",
                         specificity = c(0.8, 1),
                         ts_method = "Tau",
                         Ib = 1,
                         angle_col = 45,
                         cluster_cols = FALSE,
                         cluster_rows = FALSE,
                         show_rownames = FALSE,
                         ...) {
  color <- colorRampPalette(c("#000099",
                              "#3399FF",
                              "#3399CC",
                              "#FFFF00",
                              "#FF9900",
                              "#FF0000"))(100)
  if (!is.list(lst) | !(length(names(lst)) == 3)) {
    stop("lst maybe fault input data!")
  }

  if (length(dat_type) == 1) {
    if (dat_type == "raw") {
      reg <- range(specificity)
      if (is.numeric(specificity) & (length(reg) == 2) & (reg[1] >= 0) & (reg[2] <= 1)) {
        df <- lst$raw
        tissue_num <- ncol(df) - 11
        df <- df[order(-as.vector(df[, ts_method])), ]
        df2plot <- df[which((df[, ts_method] >= reg[1]) & (df[, ts_method] <= reg[2])), 1:tissue_num]
        main = paste("specificity of ", ts_method, "(", reg[1], "-", reg[2], ") and nrow =", nrow(df2plot))
        pheatmap(mat = df2plot,
                 ...,
                 angle_col = angle_col,
                 cluster_cols = cluster_cols,
                 cluster_rows = cluster_rows,
                 show_rownames = show_rownames,
                 color = color,
                 main = main)
      } else {
        stop("specificity must be numeric vector and range from 0 to 1!")
      }
    } else if (dat_type == "rank") {
      df <- lst$rank
      tissue_num <- ncol(df) - 2
      if (is.numeric(Ib) & (Ib > 0) & (Ib < tissue_num)) {
        df2plot <- df[which(df$Ib == Ib), 1:tissue_num]
        if (Ib < tissue_num/2) {
          df2plot <- sort_dat_de(df2plot)
        } else {
          df2plot <- sort_dat_in(df2plot)
        }
        main = paste("Ib =", Ib, "and nrow =", nrow(df2plot))
        pheatmap(mat = df2plot,
                 ...,
                 angle_col = angle_col,
                 cluster_cols = cluster_cols,
                 cluster_rows = cluster_rows,
                 show_rownames = show_rownames,
                 color = color,
                 main = main)
      } else {
          stop("Ib must be integer and range from 1 to tissue numbers!")
      }
    }
  } else if (length(dat_type) == 2) {
    if (all(dat_type %in% c("rank", "raw"))) {
      df1 <- lst$raw
      df2 <- lst$rank
      tissue_num <- ncol(df2) - 2
      if (is.numeric(Ib) & (Ib > 0) & (Ib < tissue_num)) {
        df2plot <- df1[which(df2$Ib == Ib), 1:tissue_num]
        if (Ib < tissue_num/2) {
          df2plot <- sort_dat_de(df2plot)
        } else {
          df2plot <- sort_dat_in(df2plot)
        }
        main = paste("Ib =", Ib, "and nrow =", nrow(df2plot))
        pheatmap(mat = df2plot,
                 ...,
                 angle_col = angle_col,
                 cluster_cols = cluster_cols,
                 cluster_rows = cluster_rows,
                 show_rownames = show_rownames,
                 color = color,
                 main = main)
      } else {
        stop("Ib must be integer and range from 1 to tissue numbers!")
      }
    }
  } else {
    stop("dat_type must be one of 'raw', 'rank' or c('raw', 'rank')!")
  }
}

