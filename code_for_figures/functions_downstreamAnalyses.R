library(RColorBrewer)

# Prepare data for the volcano plot ======================================
prep_volcano_plot <- function(deseq_res = deseq2_WT_vs_KO_all,
                              padj_threshold = 0.05){
  out <- deseq_res
  out <- out[with(out, order(log2FoldChange)),]
  out <- out[which(!is.na(out$padj)),]
  out$threshold  <- FALSE
  out[out$padj < padj_threshold,]$threshold <- TRUE
  out$gene <- row.names(out)
  return(out)
}

# Scale rows =========================================================
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

# Heatmap =============================================================
#' ABC's modified pheatmap
#' @param val_mat 
#' @param n_quant_breaks Integer specifying the number of color breaks. Set to 
#' \code{NULL} to turn off. Default: 150
#' @param colors vector of color names.
#' @param zero_color If zero is supposed to have a certain color assigned to it, 
#' use this parameter. Default is "white". Turn of by setting to \code{NULL}.
#' @param scale whether values should be shown as is or as z-scales.
#' @param ... additional parameters for pheatmap
abch <- function(val_mat, n_quant_breaks = 150, 
                 colors = rev(brewer.pal(n = 7, name = "RdYlBu")),
                 zero_color = "white",
                 scale = "none", ...){
  mat_plt <- as.matrix(val_mat)
  
  if (!is.null(n_quant_breaks)) {
    if (scale != "none") {
      mat4breaks <- pheatmap:::scale_mat(mat_plt, scale)
    }
    else {
      mat4breaks <- mat_plt
    }
    mat_breaks <- quantile(mat4breaks, probs = seq(0, 1, length.out = n_quant_breaks))
    mat_breaks <- mat_breaks[!duplicated(mat_breaks)]
    
    if(is.null(zero_color)){
      cols <- colorRampPalette(colors)(length(mat_breaks) - 1)
    }else{
      n_cols <- length(colors)
      grad.low <- gplots::colorpanel( sum( mat_breaks[-1] <=0 ),
                                      colors[1],
                                      colors[floor(n_cols/2)],
                                      zero_color)
      grad.high = gplots::colorpanel( sum( mat_breaks[-1] >= 0 ), 
                                      zero_color,
                                      colors[(ceiling(n_cols/2)+1)],
                                      colors[n_cols] )
      cols <- c(grad.low, grad.high)
    }
  }
  else {
    cols = colorRampPalette(colors)(100)
    mat_breaks <- NA
  }
  pheatmap::pheatmap(mat_plt, breaks = mat_breaks, color = cols, 
                     scale = scale, ...)
}
