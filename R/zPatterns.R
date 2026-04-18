zPatterns <- function(X, label = NULL, plot = TRUE,
                       axis.labels = c("Component", "Pattern ID"),
                       bar.ordered = as.character(c(FALSE, FALSE)),
                       bar.colors = c("red3", "red3"), bar.labels = FALSE, 
                       show.means = FALSE, round.means = 2, cex.means = 1,
                       type.means = c("cgm", "am"),
                       cell.colors = c("dodgerblue", "white"), 
                       cell.labels = c(label, paste("No", label)),
                       cex.axis = 1.1, grid.color = "black",
                       grid.lty = "dotted", legend = TRUE, suppress.print = FALSE, ...) {
  
  type.means <- match.arg(type.means)
  
  if (any(X < 0, na.rm = TRUE)) stop("X contains negative values")  
  if (is.vector(X)) stop("X must be a matrix or data.frame class object")
  if (is.null(label)) stop("A value for label must be given")
  
  if (!is.na(label)) {
    if (!any(X == label, na.rm = TRUE)) stop(paste("Label", label, "was not found in the data set"))
    if (label != 0 && any(X == 0, na.rm = TRUE)) warning("Unidentified zero values were found and will be ignored")
    if (any(is.na(X))) warning("Unidentified NA values were found in the data set and will be ignored")
  } else {
    if (any(X == 0, na.rm = TRUE)) warning("Unidentified zero values were found in the data set and will be ignored")
    if (!any(is.na(X))) stop(paste("Label", label, "was not found in the data set"))
  }
  
  X_mat <- as.matrix(X)
  n <- nrow(X_mat)
  p <- ncol(X_mat)
  
  if (is.na(label)) {
    miss <- is.na(X_mat)
  } else {
    miss <- !is.na(X_mat) & X_mat == label
  }
  
  pat_chars <- do.call(paste, c(as.data.frame(miss * 1L), sep = ""))
  pat_table <- table(pat_chars)
  unique_pats <- names(pat_table)
  
  pat_ID_vec <- match(pat_chars, unique_pats)
  pat_ID_factor <- factor(pat_ID_vec, levels = 1:length(unique_pats))
  
  pat_freq <- round(as.vector(pat_table) / n * 100, 2)
  prop_col <- round(colSums(miss * 1L) / n * 100, 2)
  prop_total <- round(sum(miss * 1L) / (n * p) * 100, 2)
  
  unique_idx <- match(unique_pats, pat_chars)
  tab_num <- miss[unique_idx, , drop = FALSE] * 1L
  rownames(tab_num) <- as.character(1:length(unique_pats))
  colnames(tab_num) <- colnames(X)
  
  pat_ID_levels <- rownames(tab_num)
  
  bar_ordered_bool <- as.logical(bar.ordered)
  if (is.na(bar_ordered_bool[1])) bar_ordered_bool[1] <- FALSE
  if (is.na(bar_ordered_bool[2])) bar_ordered_bool[2] <- FALSE
  
  if (bar_ordered_bool[1]) {
    ord1 <- order(pat_freq, decreasing = TRUE)
    tab_num <- tab_num[ord1, , drop = FALSE]
    pat_ID_levels <- pat_ID_levels[ord1]
    pat_freq <- pat_freq[ord1]
  }
  
  if (bar_ordered_bool[2]) {
    ord2 <- order(prop_col, decreasing = TRUE)
    tab_num <- tab_num[, ord2, drop = FALSE]
    X_mat <- X_mat[, ord2, drop = FALSE]
    prop_col <- prop_col[ord2]
  }
  
  X_NA_ignored <- X_mat
  if (!is.na(label)) {
    X_NA_ignored[!is.na(X_mat) & X_mat == label] <- NA
  }
  X_NA_ignored[!is.na(X_NA_ignored) & X_NA_ignored == 0] <- NA
  
  if (plot) {
    zones <- matrix(c(2, 4, 1, 3), ncol = 2, byrow = TRUE)
    layout(zones, widths = c(4/5, 1.5/5), heights = c(2/5, 3.5/5))
    
    par(mar = c(3, 3, 0.5, 0.5))
    a <- tab_num[nrow(tab_num):1, , drop = FALSE]
    image(1:ncol(a), 1:nrow(a), t(a), col = rev(cell.colors), axes = FALSE)
    mtext(side = 1, text = axis.labels[1], line = 1.75)
    mtext(side = 2, text = axis.labels[2], line = 1.75)
    par(mgp = c(3, 0.3, 0))
    axis(side = 1, at = 1:ncol(a), labels = colnames(a), tck = 0, cex.axis = cex.axis)
    axis(side = 2, at = 1:nrow(a), labels = rownames(a), las = 2, tck = 0, cex.axis = cex.axis)
    box()
    grid(ncol(a), nrow(a), col = grid.color, lty = grid.lty)
    
    if (show.means) {
      agg_list <- by(X_NA_ignored, pat_ID_factor, function(chunk) {
        if (type.means == "cgm") {
          ms <- apply(chunk, 2, function(x) if (all(is.na(x))) NA else exp(mean(log(x), na.rm = TRUE)))
          ms[is.na(ms)] <- 0
          s <- sum(ms)
          if (s > 0) round(ms / s * 100, round.means) else rep(0, length(ms))
        } else {
          apply(chunk, 2, function(x) if (all(is.na(x))) NA else mean(x, na.rm = TRUE))
        }
      })
      
      agg_mat <- do.call(rbind, agg_list)
      if (type.means == "am") agg_mat <- round(agg_mat, round.means)
      
      agg_mat <- agg_mat[match(pat_ID_levels, rownames(agg_mat)), , drop = FALSE]
      agg_mat_rev <- agg_mat[nrow(agg_mat):1, , drop = FALSE]
      
      if (type.means == "cgm") {
        agg_mat_rev[agg_mat_rev == 0] <- NA
      }
      
      x_coords <- col(agg_mat_rev)
      y_coords <- row(agg_mat_rev)
      labels_vec <- as.vector(agg_mat_rev)
      valid <- !is.na(labels_vec)
      
      text(x_coords[valid], y_coords[valid], labels = labels_vec[valid], cex = cex.means)
    }
    
    par(mar = c(0, 3.25, 1, 0.75))
    barplot_top <- barplot(as.vector(prop_col), axes = FALSE, col = bar.colors[1], xaxs = "i",
                           ylim = c(0, max(prop_col) * 1.2))
    if (bar.labels) {
      text(barplot_top, as.vector(prop_col), labels = as.vector(prop_col), cex = 0.85, pos = 3)
    }
    
    par(mar = c(3.25, 0, 0.75, 0.75))
    barplot_side <- barplot(rev(as.vector(pat_freq)), horiz = TRUE, axes = FALSE, col = bar.colors[2], yaxs = "i",
                            xlim = c(0, max(pat_freq) * 1.3))
    if (bar.labels) {
      text(rev(as.vector(pat_freq)), barplot_side, labels = rev(as.vector(pat_freq)), cex = 0.85, pos = 4)
    }
    
    par(mar = c(0, 0, 3, 0))
    plot.new()
    if (legend) {
      if (any(is.na(cell.labels))) cell.labels[is.na(cell.labels)] <- "NA"
      legend("topleft", legend = cell.labels, pch = c(22, 22), bty = "n",
             pt.bg = cell.colors, pt.cex = 2, cex = 1.1)
    }
  }
  
  tab_char <- matrix(ifelse(tab_num == 1, "+", "-"), nrow = nrow(tab_num), dimnames = dimnames(tab_num))
  tab_print <- data.frame(Patt.ID = pat_ID_levels, tab_char, No.Unobs = rowSums(tab_num == 1), 
                          Patt.Perc = pat_freq, stringsAsFactors = FALSE)
  
  if (!suppress.print) {
    cat("Patterns ('+' means ", cell.labels[1], ", '-' means ", cell.labels[2], ") \n\n", sep = "")
    print(tab_print, row.names = FALSE)
    cat("\nPercentage cells by component \n")
    print(prop_col)
    cat(sprintf("\nOverall percentage cells: %s%% \n", prop_total))
  }
  
  invisible(factor(pat_ID_vec, levels = pat_ID_levels))
}
