
# data manipulation libraries 
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(tidyr) 
library(magrittr)
library(jsonlite)
library(purrr)
library(rlang) # for as_name

# plotting 
library(ggplot2)
library(scales) # help with axises 
library(GGally) # ggpairs for pairwise comparisons 

# interactive html outputs
library(plotly) # ggplotly to turn any ggplot object into an interactive plot 
library(DT) # data tables for displaying tables 

# stats package
library(DescTools) # CCC, Calculates Lin's concordance correlation coefficient


### Some helpful R utilities 


#| code-fold: true

# commonly used functions  
`%!in%` <- Negate(`%in%`)
is.notna <-function(x){
  return(!is.na(x))
}

# ggplot themes 
transparentBackground = theme(
  panel.background = element_rect(fill='transparent'), #transparent panel bg
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  # panel.grid.major = element_blank(), #remove major gridlines
  # panel.grid.minor = element_blank(), #remove minor gridlines
  legend.background = element_rect(color = NA, fill='transparent'), #transparent legend bg
  legend.box.background = element_rect(color = NA, fill='transparent') #transparent legend panel
)
custom_ggplot_theme_noTransparentBackground = theme_bw() +
  # theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank() )+
  theme(axis.line.x = element_line(color="black", linewidth = 0.3),axis.line.y =
          element_line(color="black", linewidth = 0.3))+
  theme(text=element_text(size=12, family="Helvetica"))+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.position = "bottom") + 
  theme(plot.title = element_text(hjust = 0.5))

custom_ggplot_theme = custom_ggplot_theme_noTransparentBackground +
  transparentBackground
custom_ggplot_theme_backgroundTransparent = custom_ggplot_theme

custom_ggplot_theme_xRotate_noBackgroundTransparent = theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank() )+
  theme(axis.line.x = element_line(color="black", linewidth = 0.3),axis.line.y =
          element_line(color="black", linewidth = 0.3))+
  theme(text=element_text(size=12, family="Helvetica"))+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.position = "bottom") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(size=12, angle = -90, vjust = 0.5, hjust = 0)) 

custom_ggplot_theme_xRotate = custom_ggplot_theme_xRotate_noBackgroundTransparent + 
  transparentBackground

custom_ggplot_theme_xRotate_backgroundTransparent = custom_ggplot_theme_xRotate

# color blind safe distinct color palettes 
colorPalette_08 = c("#2271B2","#F748A5","#359B73","#F0E442","#D55E00","#3DB7E9","#E69F00","#000000")
colorPalette_12 = c("#E20134","#FF6E3A","#008DF9","#8400CD","#FFC33B","#9F0162","#009F81","#FF5AAF","#00FCCF","#00C2F9","#FFB2FD","#A40122")
colorPalette_15 = c("#F60239","#003C86","#EF0096","#9400E6","#009FFA","#008169","#68023F","#00DCB5","#FFCFE2","#FF71FD","#7CFFFA","#6A0213","#008607","#00E307","#FFDC3D")


# creating a datatable html widget
create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All"))), 
                filter = "top")
}

# for creating tabsets html outputs
create_tabsetOfHtmlWidgets <- function(htmlObjectsList) {
  zz <- textConnection("foo", "w")
  sink(zz)
  cat("::: {.panel-tabset}\n")
  iwalk(htmlObjectsList, ~ {
    cat('## ', .y, '\n\n')
    tempList = list()
    tempList[["item"]] = .x
    print(htmltools::tagList(tempList))
    cat('\n\n')
  })
  cat(":::\n")
  sink()
  close(zz)
  paste0(foo, collapse = "\n")
}

create_tabsetOfGgplotObjects <- function(ggplotObjectsList) {
  cat("\n\n::: {.panel-tabset}\n\n")
  iwalk(ggplotObjectsList, ~ {
    cat('## ', .y, '\n\n')
    print(.x)
    cat('\n\n')
  })
  cat(":::\n")
}




#' Add various contingency table metrics, input must have TP, FP, TN, FN columns defined already
#'
#' @param input A data table with at least TP, FP, TN, FN
#' @param summarize_first Whether or not to sum up the columns first, will lose any non-grouping columns though
#'
#' @returns the input table plus various performance metrics 
#' @export
add_contigency_metrics <- function(input, summarize_first = F){
  output = input
  if(summarize_first){
    output = output  %>%
      summarise(TP = sum(TP), 
                FN = sum(FN), 
                FP = sum(FP), 
                TN = sum(TN))
  }
  output = output %>%
    mutate(
      sensitivity           = if_else(TP + FN > 0, TP / (TP + FN), NA_real_),
      recall                = sensitivity,
      specificity           = if_else(TN + FP > 0, TN / (TN + FP), NA_real_),
      ppv                   = if_else(TP + FP > 0, TP / (TP + FP), NA_real_),
      precision             = ppv,
      npv                   = if_else(TN + FN > 0, TN / (TN + FN), NA_real_),
      accuracy              = (TP + TN) / (TP + TN + FP + FN),
      prevalence            = (TP + FN) / (TP + TN + FP + FN),
      overlap_coefficient   = ifelse(TP + FP + FN > 0, TP / (TP + FP + FN), NA_real_),
      tpr                   = sensitivity,
      fpr                   = if_else(FP + TN > 0, FP / (FP + TN), NA_real_),
      fdr                   = if_else(TP + FP > 0, FP / (TP + FP), NA_real_),
      f1_score              = if_else((2 * TP + FP + FN) > 0, (2 * TP) / (2 * TP + FP + FN), NA_real_),
      balanced_accuracy     = (sensitivity + specificity) / 2
    )
  return(output)
}


#' Add various contingency table metrics but only if they don't depend on TN, input must have TP, FP, FN columns defined already
#'
#' @param input A data table with at least TP, FP, FN
#' @param summarize_first Whether or not to sum up the columns first, will lose any non-grouping columns though
#' @returns the input table plus various performance metrics 
#' @export
add_contigency_metrics_no_TN <- function(input, summarize_first = F){
  output = input
  if(summarize_first){
    output = output  %>%
      summarise(TP = sum(TP), 
                FN = sum(FN), 
                FP = sum(FP), 
                n = n())
  }
  output = output %>%
    mutate(
      
      sensitivity           = if_else(TP + FN > 0, TP / (TP + FN), NA_real_),
      recall                = sensitivity,
      ppv                   = if_else(TP + FP > 0, TP / (TP + FP), NA_real_),
      precision             = ppv,
      overlap_coefficient   = ifelse(TP + FP + FN > 0, TP / (TP + FP + FN), NA_real_),
      tpr                   = sensitivity,
      fdr                   = if_else(TP + FP > 0, FP / (TP + FP), NA_real_),
      f1_score              = if_else((2 * TP + FP + FN) > 0, (2 * TP) / (2 * TP + FP + FN), NA_real_)
    )
  return(output)
}





#' Plot metrics that range from 0 to 1 in radar perform with optional coloring and faceting, all input must have 1 obs per df (+/- groupings)
#'
#' @param df The summarized metrics table 
#' @param color_col an optional coloring if there are multiple groups in data 
#' @param facet_col an optional faceting if there are multiple groups in data
#' @param metric_map a map to turn the metrics into nice lables, can also be a 1:1 table, is also used to indicate what metrics are being plotted 
#' @param palette the coloring palette being used if coloring, will be passed to scale_color_manual
#' @param title the title for the plot
#' @param outer_pad the padding used to make sure labels aren't clipped off 
#' @param ring_breaks where to put the breaks in the rignt 
#'
#' @returns a ggplot object of the plot which can be further modified 
#' @export
radar_by_metrics <- function(df,
                             color_col = NULL,          
                             facet_col = NULL,          
                             metric_map,                # named chr: "Nice_Col_Name_Bro" = "df_col"
                             palette = NULL,
                             title = "Performance radar",
                             outer_pad = 0.12,
                             ring_breaks = seq(0, 1, 0.25)) {
  # Requires: ggplot2, dplyr, tidyr, tibble, scales
  is_named <- function(x) !is.null(names(x)) && all(nzchar(names(x)))
  stopifnot(is.character(metric_map), is_named(metric_map))
  stopifnot(all(unname(metric_map) %in% names(df)))
  
  has_color <- is.character(color_col) && length(color_col) == 1L && nzchar(color_col)
  has_facet <- is.character(facet_col) && length(facet_col) == 1L && nzchar(facet_col)
  
  sel_cols <- c(if (has_color) color_col,
                if (has_facet) facet_col,
                unname(metric_map))
  
  metrics <- names(metric_map); k <- length(metrics)
  
  # Long format
  radar_long_ix <- df %>%
    dplyr::select(dplyr::all_of(sel_cols)) %>%
    {
      if (has_color && has_facet) dplyr::rename(., .color = dplyr::all_of(color_col),
                                                .facet = dplyr::all_of(facet_col))
      else if (has_color && !has_facet) dplyr::rename(., .color = dplyr::all_of(color_col)) %>%
        dplyr::mutate(.facet = "(all)")
      else if (!has_color && has_facet) dplyr::rename(., .facet = dplyr::all_of(facet_col)) %>%
        dplyr::mutate(.color = "(all)")
      else dplyr::mutate(., .color = "(all)", .facet = "(all)")
    } %>%
    tidyr::pivot_longer(cols = dplyr::all_of(unname(metric_map)),
                        names_to = "metric_key", values_to = "value") %>%
    dplyr::mutate(
      metric = factor(metrics[match(metric_key, metric_map)], levels = metrics),
      value  = as.numeric(value),
      x_idx  = as.integer(metric)
    ) %>%
    dplyr::select(-metric_key)
  
  # Strict checks 
  ## Strict range check: finite & in [0,1]
  bad <- radar_long_ix %>% dplyr::filter(!is.na(value) & (!is.finite(value) | value < 0 | value > 1))
  if (nrow(bad) > 0) {
    summary_msg <- bad %>%
      dplyr::group_by(metric) %>%
      dplyr::summarise(n = dplyr::n(),
                       min_val = min(value, na.rm = TRUE),
                       max_val = max(value, na.rm = TRUE),
                       .groups = "drop") %>%
      dplyr::mutate(msg = sprintf("%s (n=%d, range=[%g, %g])", metric, n, min_val, max_val)) %>%
      dplyr::pull(msg) %>%
      paste(collapse = "; ")
    stop(sprintf("Out-of-bounds metric values (expected finite in [0,1]). Offenders: %s", summary_msg),
         call. = FALSE)
  }
  
  ## Uniqueness check: exactly 0/1 non-NA per point (facet × color × metric), throw if not 
  g_cols <- c(".facet", ".color", "metric", "x_idx")
  dup_check <- radar_long_ix %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(g_cols))) %>%
    dplyr::summarise(n_non_na = sum(!is.na(value)), .groups = "drop") %>%
    dplyr::filter(n_non_na > 1)
  
  if (nrow(dup_check) > 0) {
    offenders <- dup_check %>%
      dplyr::mutate(msg = sprintf("[facet=%s, color=%s, metric=%s, count=%d]",
                                  .facet, .color, as.character(metric), n_non_na)) %>%
      dplyr::pull(msg)
    stop(
      paste0(
        "Multiple observations per point detected. The input must be pre-summarized to one row ",
        "per (facet, color, metric). Offenders (showing up to first 10): ",
        paste(head(offenders, 10), collapse = "; "),
        if (nrow(dup_check) > 10) sprintf(" ... and %d more.", nrow(dup_check) - 10) else ""
      ),
      call. = FALSE
    )
  }
  
  # Keep at most one non-NA row per point (ok after uniqueness check)
  radar_long_ix <- radar_long_ix %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(g_cols))) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  # Close each line: append first vertex at x = k+1
  close_lines <- function(d) {
    d <- dplyr::arrange(d, x_idx)
    if (nrow(d) == 0) return(d)
    first <- d[1, ]; first$x_idx <- k + 1
    dplyr::bind_rows(d, first)
  }
  
  group_vars <- c(".facet", ".color")
  radar_closed_ix <- radar_long_ix %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::group_modify(~ close_lines(.x)) %>%
    dplyr::ungroup()
  
  # Background bands
  x_seq <- c(seq_len(k), k + 1)
  bands <- tibble::tibble(
    band_id = factor(seq_len(length(ring_breaks) - 1)),
    ymin = head(ring_breaks, -1),
    ymax = tail(ring_breaks, -1)
  ) %>% tidyr::crossing(tibble::tibble(x = x_seq))
  
  # Build plot
  p <- ggplot() +
    geom_ribbon(
      data = bands,
      aes(x = x, ymin = ymin, ymax = ymax, group = band_id),
      inherit.aes = FALSE, alpha = 0.08, fill = "grey50", color = NA
    ) +
    {
      if (has_color) {
        list(
          geom_line(
            data = dplyr::arrange(radar_closed_ix, .facet, .color, x_idx),
            aes(x = x_idx, y = value, color = .color, group = .color),
            linewidth = 0.8, show.legend = TRUE
          ),
          geom_point(
            data = radar_long_ix,
            aes(x = x_idx, y = value, color = .color, group = .color),
            size = 2, show.legend = TRUE
          )
        )
      } else {
        list(
          geom_line(
            data = dplyr::arrange(radar_closed_ix, .facet, x_idx),
            aes(x = x_idx, y = value, group = 1),
            linewidth = 0.9, color = "black", show.legend = FALSE
          ),
          geom_point(
            data = radar_long_ix,
            aes(x = x_idx, y = value),
            size = 2, color = "black", show.legend = FALSE
          )
        )
      }
    } +
    coord_polar(clip = "off") +
    scale_x_continuous(limits = c(1, k + 1), breaks = seq_len(k), labels = metrics) +
    scale_y_continuous(limits = c(0, 1 + outer_pad),
                       breaks = ring_breaks,
                       labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0))) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3),
      axis.text.x = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 9),
      legend.position = if (has_color) "bottom" else "none",
      strip.text = element_text(face = "bold"),
      plot.margin = margin(10, 30, 10, 30)
    )
  
  # Facet if requested
  if (has_facet) p <- p + facet_wrap(~ .facet)
  
  # Palette / legend title
  if (has_color && !is.null(palette)) {
    p <- p + scale_color_manual(values = palette, name = color_col)
  } else if (has_color) {
    p <- p + labs(color = color_col)
  }
  
  return(p)
}




# custom hex plot function for GGally
my_hex <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_hex(aes(fill = after_stat(count)), ...) +  # map fill here
    scale_fill_viridis_c(trans = "log10", name = "count") +  # log scale
    theme_minimal() +
    geom_smooth() + 
    # geom_smooth(method = 'gam', formula = 'y ~ s(x, bs = "cs")') +
    geom_abline(intercept = 0,
                slope = 1,
                color = "red")
}

# custom correlation and CCC 
my_corr_ccc <- function(data, mapping,
                        method = "pearson",
                        digits = 3,
                        ccc_conf = 0.95,
                        size = 4, ...) {
  # extract columns from mapping
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  
  label <- "Not enough data"
  if (length(x) >= 2 && stats::sd(x) > 0 && stats::sd(y) > 0) {
    r <- suppressWarnings(stats::cor(x, y, method = method, use = "complete.obs"))
    ccc <- try(DescTools::CCC(x, y, ci = "z-transform", conf.level = ccc_conf), silent = TRUE)
    
    if (inherits(ccc, "try-error")) {
      label <- sprintf("cor = %s\n\nCCC = NA", signif(r, digits))
    } else {
      label <- paste0(
        method, " cor = ", signif(r, digits), "\n\n",
        "CCC = ",
        paste0(
          signif(ccc$rho.c$est, digits), " \n(95% CI =",
          signif(ccc$rho.c$lwr.ci, digits), "–",
          signif(ccc$rho.c$upr.ci, digits), ")"
        )
      )
    }
  }
  
  # draw centered text in the panel
  ggplot(data.frame(x = 0.5, y = 0.5, label = label), aes(x, y, label = label)) +
    geom_text(size = size, lineheight = 1.05, ...) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void()
}




#' Add columns of mean of observed and expected and the observed minus the expected 
#'
#' @param df the dataframe to add to, needs an observation column and an expected column 
#' @param expected_col character string of the expected column
#' @param observed_col character string of the observed column
#'
#' @returns the input with the difference from expected and mean of expected and observed
#' @export
add_diff_mean <- function(df, expected_col, observed_col) {
  expected_col <- enquo(expected_col)
  observed_col <- enquo(observed_col)
  
  obs_name <- as_name(observed_col)
  
  df %>%
    rowwise() %>%
    mutate(
      !!paste0(obs_name, "_diff_from_exp") := !!observed_col - !!expected_col,
      !!paste0(obs_name, "_mean_with_exp") := mean(c(!!observed_col, !!expected_col))
    ) %>%
    ungroup()
}



#' Create a Bland–Altman style plot from observed vs expected COI measures
#'
#' This function generates a hex-binned Bland–Altman plot for a given COI measure
#' compared against an expected COI. It automatically looks up the mean and diff
#' columns based on a provided label (e.g., `"ecoi"`, `"coi"`, `"naive_coi"`, 
#' `"offset_naive_coi"`).
#'
#' @param df A data frame containing the mean and diff columns.
#' @param label A character string specifying the observed measure name. 
#'   The function will look for columns named `paste0(label, mean_suffix)` 
#'   and `paste0(label, diff_suffix)`.
#' @param mean_suffix A character string appended to the label to form the
#'   mean column name. Default is `"_mean_with_exp"`.
#' @param diff_suffix A character string appended to the label to form the
#'   diff column name. Default is `"_diff_from_exp"`.
#' @param expected_label A character string used in axis labels to describe
#'   the expected measure. Default is `"Expected COI"`.
#' @param bins Integer, number of bins for the hex plot. Default is `30`.
#'
#' @returns A `ggplot` object representing the Bland–Altman style plot.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example with ecoi
#' p_ecoi <- make_bland_altman_plot_by_label(moire_all_coi_summary, "ecoi")
#' 
#' # Example with naive_coi
#' p_naive <- make_bland_altman_plot_by_label(moire_all_coi_summary, "naive_coi")
#' 
#' # Example with custom suffixes (if column names differ)
#' p_custom <- make_bland_altman_plot_by_label(
#'   df = moire_all_coi_summary,
#'   label = "ecoi",
#'   mean_suffix = "_mean",
#'   diff_suffix = "_diff"
#' )
#' }
make_bland_altman_plot_by_label <- function(
    df,
    label,
    mean_suffix = "_mean_with_exp",
    diff_suffix = "_diff_from_exp",
    expected_label = "Expected",
    bins = 30
) {
  stopifnot(is.data.frame(df), is.character(label), length(label) == 1)
  
  mean_col <- paste0(label, mean_suffix)
  diff_col <- paste0(label, diff_suffix)
  
  if (!(mean_col %in% names(df))) stop(sprintf("Column '%s' not found in df.", mean_col))
  if (!(diff_col %in% names(df))) stop(sprintf("Column '%s' not found in df.", diff_col))
  
  ggplot(df) +
    geom_hex(
      aes(x = .data[[mean_col]], y = .data[[diff_col]], fill = log10(after_stat(count))),
      bins = bins
    ) +
    scale_fill_viridis_c(name = "log10(count)") +
    theme_minimal() +
    geom_smooth(
      aes(x = .data[[mean_col]], y = .data[[diff_col]]),
      method = "gam", se = FALSE
    ) +
    # +2 SD per facet
    stat_summary(
      aes(x = 1, y = .data[[diff_col]], yintercept = after_stat(y)),
      fun = function(y)  2 * stats::sd(y, na.rm = TRUE),
      geom = "hline",
      inherit.aes = TRUE,
      linetype = "dashed",
      color = "grey30"
    ) +
    # -2 SD per facet
    stat_summary(
      aes(x = 1, y = .data[[diff_col]], yintercept = after_stat(y)),
      fun = function(y) -2 * stats::sd(y, na.rm = TRUE),
      geom = "hline",
      inherit.aes = TRUE,
      linetype = "dashed",
      color = "grey30"
    ) +
    labs(
      title = label,
      y = paste0(label, " - ", expected_label),
      x = paste0("mean (", label, ", ", expected_label, ")")
    ) +
    coord_equal()
}




#' @title Concordance / Error Metrics (single or multiple observed columns)
#' @keywords CCC RMSE MAE COI
#'
#' @description
#' Compute metrics between an expected column and one or more observed columns.
#' Each function accepts an optional set of grouping columns. If omitted, an
#' overall (single-row) summary is returned.
#'
#' @param df A data frame.
#' @param expected_col Unquoted column name of expected values.
#' @param observed_cols One symbol or \code{c(sym1, sym2, ...)} of observed columns.
#' @param group_cols Optional grouping columns supplied like
#'   \code{c(population_name, paramset_id)}. If omitted or \code{NULL},
#'   computes overall metrics (no grouping).
#'
#' @return A tibble with group columns (if any) and metric columns per observed variable.
#'
#' @importFrom dplyr group_by summarise distinct left_join bind_cols
#' @importFrom rlang enquo enexpr is_call call_args as_quosure caller_env as_name
#' @importFrom purrr map map_dbl
#' @name vectorized_metrics


# ---- CCC (DescTools::CCC) ----------------------------------------------------

#' @rdname vectorized_metrics
#' @export
compute_ccc <- function(df, expected_col, observed_cols, group_cols = NULL) {
  stopifnot(is.data.frame(df))
  expected_col <- rlang::enquo(expected_col)
  
  # helper: TRUE if vectors match perfectly (after dropping NA pairs)
  is_perfect_match <- function(x, y) {
    ok <- stats::complete.cases(x, y)
    x <- x[ok]; y <- y[ok]
    length(x) > 0 && isTRUE(all.equal(x, y, tolerance = 0))
  }
  
  # helper: compute CCC, but return 1/1/1 when perfectly matching (or when CCC yields NA)
  ccc_safe <- function(x, y) {
    if (is_perfect_match(x, y)) {
      return(list(est = 1, lwr = 1, upr = 1))
    }
    res <- DescTools::CCC(x, y)
    est <- res$rho.c$est
    lwr <- res$rho.c$lwr.ci
    upr <- res$rho.c$upr.ci
    
    # guard: if DescTools gives NA even though it's a perfect match (or for numerical edge cases)
    if (is.na(est) && is_perfect_match(x, y)) {
      return(list(est = 1, lwr = 1, upr = 1))
    }
    
    list(est = est, lwr = lwr, upr = upr)
  }
  
  # capture observed columns (symbol or c(...))
  obs_expr <- rlang::enexpr(observed_cols)
  obs_quos <- if (rlang::is_call(obs_expr, "c")) {
    purrr::map(rlang::call_args(obs_expr), rlang::as_quosure, env = rlang::caller_env())
  } else {
    list(rlang::as_quosure(obs_expr, env = rlang::caller_env()))
  }
  
  # capture grouping columns (NULL/omitted => no groups)
  grp_expr <- rlang::enexpr(group_cols)
  grp_quos <- if (missing(group_cols) || is.null(grp_expr)) {
    list()
  } else if (rlang::is_call(grp_expr, "c")) {
    purrr::map(rlang::call_args(grp_expr), rlang::as_quosure, env = rlang::caller_env())
  } else {
    list(rlang::as_quosure(grp_expr, env = rlang::caller_env()))
  }
  
  if (length(grp_quos) > 0) {
    out <- df %>% dplyr::distinct(!!!grp_quos)
    join_by_cols <- names(out)
    
    for (q in obs_quos) {
      obs_name <- rlang::as_name(q)
      
      one <- df %>%
        dplyr::group_by(!!!grp_quos) %>%
        dplyr::summarise(
          .ccc = list(ccc_safe(rlang::eval_tidy(expected_col, data = dplyr::cur_data()),
                               rlang::eval_tidy(q,            data = dplyr::cur_data()))),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          !!paste0(obs_name, "_ccc")       := purrr::map_dbl(.ccc, ~ .x$est),
          !!paste0(obs_name, "_ccc_lower") := purrr::map_dbl(.ccc, ~ .x$lwr),
          !!paste0(obs_name, "_ccc_upper") := purrr::map_dbl(.ccc, ~ .x$upr)
        ) %>%
        dplyr::select(-.ccc)
      
      out <- dplyr::left_join(out, one, by = join_by_cols)
    }
  } else {
    pieces <- purrr::map(obs_quos, function(q) {
      obs_name <- rlang::as_name(q)
      
      x <- rlang::eval_tidy(expected_col, data = df)
      y <- rlang::eval_tidy(q,            data = df)
      v <- ccc_safe(x, y)
      
      tibble::tibble(
        !!paste0(obs_name, "_ccc")       := v$est,
        !!paste0(obs_name, "_ccc_lower") := v$lwr,
        !!paste0(obs_name, "_ccc_upper") := v$upr
      )
    })
    
    out <- dplyr::bind_cols(pieces)
  }
  
  out
}

# ---- RMSE --------------------------------------------------------------------

#' @rdname vectorized_metrics
#' @export
compute_rmse <- function(df, expected_col, observed_cols, group_cols = NULL) {
  stopifnot(is.data.frame(df))
  expected_col <- rlang::enquo(expected_col)
  
  obs_expr <- rlang::enexpr(observed_cols)
  obs_quos <- if (rlang::is_call(obs_expr, "c")) {
    purrr::map(rlang::call_args(obs_expr), rlang::as_quosure, env = rlang::caller_env())
  } else {
    list(rlang::as_quosure(obs_expr, env = rlang::caller_env()))
  }
  
  grp_expr <- rlang::enexpr(group_cols)
  grp_quos <- if (missing(group_cols) || is.null(grp_expr)) {
    list()
  } else if (rlang::is_call(grp_expr, "c")) {
    purrr::map(rlang::call_args(grp_expr), rlang::as_quosure, env = rlang::caller_env())
  } else {
    list(rlang::as_quosure(grp_expr, env = rlang::caller_env()))
  }
  
  metric_fun <- function(q) {
    obs_name <- rlang::as_name(q)
    if (length(grp_quos) > 0) {
      df %>%
        dplyr::group_by(!!!grp_quos) %>%
        dplyr::summarise(
          !!paste0(obs_name, "_rmse") := {
            x <- !!q; y <- !!expected_col
            ok <- !is.na(x) & !is.na(y)
            n  <- sum(ok)
            if (n == 0) NA_real_ else sqrt(sum((x[ok] - y[ok])^2) / n)
          },
          .groups = "drop"
        )
    } else {
      dplyr::summarise(df,
                       !!paste0(obs_name, "_rmse") := {
                         x <- !!q; y <- !!expected_col
                         ok <- !is.na(x) & !is.na(y)
                         n  <- sum(ok)
                         if (n == 0) NA_real_ else sqrt(sum((x[ok] - y[ok])^2) / n)
                       }
      )
    }
  }
  
  if (length(grp_quos) > 0) {
    out <- df %>% dplyr::distinct(!!!grp_quos)
    join_by_cols <- names(out)
    for (q in obs_quos) {
      out <- dplyr::left_join(out, metric_fun(q), by = join_by_cols)
    }
  } else {
    out <- dplyr::bind_cols(purrr::map(obs_quos, metric_fun))
  }
  
  out
}

# ---- MAE ---------------------------------------------------------------------

#' @rdname vectorized_metrics
#' @export
compute_mae <- function(df, expected_col, observed_cols, group_cols = NULL) {
  stopifnot(is.data.frame(df))
  expected_col <- rlang::enquo(expected_col)
  
  obs_expr <- rlang::enexpr(observed_cols)
  obs_quos <- if (rlang::is_call(obs_expr, "c")) {
    purrr::map(rlang::call_args(obs_expr), rlang::as_quosure, env = rlang::caller_env())
  } else {
    list(rlang::as_quosure(obs_expr, env = rlang::caller_env()))
  }
  
  grp_expr <- rlang::enexpr(group_cols)
  grp_quos <- if (missing(group_cols) || is.null(grp_expr)) {
    list()
  } else if (rlang::is_call(grp_expr, "c")) {
    purrr::map(rlang::call_args(grp_expr), rlang::as_quosure, env = rlang::caller_env())
  } else {
    list(rlang::as_quosure(grp_expr, env = rlang::caller_env()))
  }
  
  metric_fun <- function(q) {
    obs_name <- rlang::as_name(q)
    if (length(grp_quos) > 0) {
      df %>%
        dplyr::group_by(!!!grp_quos) %>%
        dplyr::summarise(
          !!paste0(obs_name, "_mean_abs_error") := {
            x <- !!q; y <- !!expected_col
            ok <- !is.na(x) & !is.na(y)
            if (!any(ok)) NA_real_ else mean(abs(x[ok] - y[ok]))
          },
          .groups = "drop"
        )
    } else {
      dplyr::summarise(df,
                       !!paste0(obs_name, "_mean_abs_error") := {
                         x <- !!q; y <- !!expected_col
                         ok <- !is.na(x) & !is.na(y)
                         if (!any(ok)) NA_real_ else mean(abs(x[ok] - y[ok]))
                       }
      )
    }
  }
  
  if (length(grp_quos) > 0) {
    out <- df %>% dplyr::distinct(!!!grp_quos)
    join_by_cols <- names(out)
    for (q in obs_quos) {
      out <- dplyr::left_join(out, metric_fun(q), by = join_by_cols)
    }
  } else {
    out <- dplyr::bind_cols(purrr::map(obs_quos, metric_fun))
  }
  
  out
}



quick_summary <- function(x, probs = c(0, 0.25, 0.5, 0.75, 1)) {
  q <- quantile(x, probs = probs, na.rm = TRUE)
  tibble(
    n = length(x[!is.na(x)]), 
    sd = sd(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    min  = min(x, na.rm = TRUE),
    max  = max(x, na.rm = TRUE),
    !!!setNames(as.list(q), paste0("q", probs*100))
  )
}




