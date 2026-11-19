# Shared data prep for the between-tool MAD4HATTER pages.
# Parameter sets are chosen by Lin's concordance. The benchmark root can be
# overridden with PGEFORGE_BENCHMARK_ROOT when a newer snapshot is rendered.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(jsonlite)
})

options(readr.show_progress = FALSE, dplyr.summarise.inform = FALSE)

.bench <- new.env(parent = emptyenv())

tool_colours <- c(
  "Naive COI" = "#325b74",
  "MOIRE" = "#c82f61",
  "TheREALMcCOIL" = "#ebab35",
  "coiaf" = "#2a9d8f",
  "SNP-slice" = "#6d597a",
  "Naive single-locus allele frequency" = "#325b74",
  "Incomplete data model" = "#e07a3d",
  "dcifer" = "#3d5a80",
  "Naive microhaplotype frequency" = "#1d7874",
  "Naive multi-locus frequency" = "#325b74",
  "FreqEstimationModel" = "#c82f61",
  "MultiLociBiallelicModel" = "#ebab35"
)

tool_levels <- names(tool_colours)

benchmark_root <- function() {
  env <- Sys.getenv("PGEFORGE_BENCHMARK_ROOT", unset = "")
  candidates <- c(
    env,
    "data/benchmarking/MAD4HATTER",
    file.path("..", "data", "benchmarking", "MAD4HATTER")
  )
  candidates <- candidates[nzchar(candidates) & dir.exists(candidates)]
  if (!length(candidates)) {
    stop("Benchmark data not found. Set PGEFORGE_BENCHMARK_ROOT to the MAD4HATTER snapshot.")
  }
  normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
}

panel_dir <- function(panel) {
  if (!panel %in% c("full", "subselection")) stop("Unknown panel: ", panel)
  file.path(benchmark_root(), panel)
}

coi_tool_order <- c("Naive COI", "MOIRE", "TheREALMcCOIL", "coiaf", "SNP-slice")

# Colours from the manuscript tool table for the COI comparison.
coi_colours <- c(
  "Naive COI" = "#009F81",
  "MOIRE" = "#FF5AAF",
  "TheREALMcCOIL" = "#FF6E3A",
  "coiaf" = "#008DF9",
  "SNP-slice" = "#8400CD"
)
variant_tool_order <- c(
  "Naive single-locus allele frequency", "Incomplete data model", "dcifer"
)
# Colours from the manuscript tool table for the single-locus frequency comparison.
variant_colours <- c(
  "Naive single-locus allele frequency" = "#009F81",
  "Incomplete data model" = "#008DF9",
  "dcifer" = "#9F0162"
)
mhap_tool_order <- c("MOIRE", "dcifer", "Naive microhaplotype frequency")
# Colours from the manuscript tool table for the microhaplotype comparison.
mhap_colours <- c(
  "MOIRE" = "#FF5AAF",
  "dcifer" = "#9F0162",
  "Naive microhaplotype frequency" = "#009F81"
)
multi_tool_order <- c(
  "Naive multi-locus frequency", "FreqEstimationModel",
  "MultiLociBiallelicModel", "SNP-slice"
)
# Colours from the manuscript tool table for the multi-locus comparison.
multi_colours <- c(
  "Naive multi-locus frequency" = "#009F81",
  "FreqEstimationModel" = "#FF5AAF",
  "MultiLociBiallelicModel" = "#00FCCF",
  "SNP-slice" = "#8400CD"
)

order_tools <- function(x, preferred = tool_levels) {
  present <- unique(as.character(x))
  levels <- c(preferred[preferred %in% present], setdiff(present, preferred))
  factor(x, levels = levels)
}

theme_bench <- function(legend = "none") {
  theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "#f9eaf1", colour = NA),
      panel.grid.minor = element_blank(),
      legend.position = legend
    )
}

lin_ccc <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3) return(NA_real_)
  sx <- stats::sd(x)
  sy <- stats::sd(y)
  if (!is.finite(sx) || !is.finite(sy) || sx == 0 || sy == 0) return(NA_real_)
  r <- stats::cor(x, y)
  if (!is.finite(r)) return(NA_real_)
  2 * r * sx * sy / (sx^2 + sy^2 + (mean(x) - mean(y))^2)
}

sample_key <- function(x) sub("_S[0-9]+$", "", as.character(x))

below_one_percent <- function(x) !is.finite(x) | x < (0.01 - 1e-8)

read_chars <- function(path, columns) {
  read_tsv(
    path,
    col_types = cols(.default = col_character()),
    col_select = any_of(columns),
    progress = FALSE,
    show_col_types = FALSE
  )
}

kept_populations <- function() {
  if (!is.null(.bench$pops)) return(.bench$pops)
  root <- file.path(benchmark_root(), "simulation_setup")
  dirs <- list.dirs(root, recursive = FALSE, full.names = FALSE)
  dirs <- dirs[grepl("^simulated_(Ghana|Cambodia)_pop_[0-9]+$", dirs)]
  meta <- lapply(dirs, function(d) {
    j <- fromJSON(file.path(root, d, "simulation_parameters.json"))
    data.frame(
      population_name = as.character(j$population_name),
      coi_cat = as.character(j$coi_cat),
      cotransmission_cat = as.character(j$cotransmission_cat),
      pop_alpha = as.numeric(j$pop_alpha),
      stringsAsFactors = FALSE
    )
  })
  pops <- bind_rows(meta)
  drop <- pops$cotransmission_cat == "No_Cotransmission" |
    (pops$coi_cat == "MeanCOI=1.26" & pops$cotransmission_cat == "High_Cotransmission") |
    (pops$coi_cat == "MeanCOI=5.47" & pops$cotransmission_cat == "Low_Cotransmission")
  pops <- pops[!drop, , drop = FALSE]
  pops$relatedness <- case_when(
    abs(pops$pop_alpha - 999) < 1e-6 ~ "0.1% related",
    abs(pops$pop_alpha - 99) < 1e-6 ~ "1% related",
    abs(pops$pop_alpha - 2.33) < 1e-2 ~ "30% related",
    TRUE ~ NA_character_
  )
  if (anyNA(pops$relatedness)) {
    stop(
      "Unrecognised pop_alpha values: ",
      paste(unique(pops$pop_alpha[is.na(pops$relatedness)]), collapse = ", ")
    )
  }
  pops$mean_coi <- sub("^MeanCOI=", "Mean COI ", pops$coi_cat)
  pops$cotransmission <- case_when(
    pops$cotransmission_cat == "Low_Cotransmission" ~ "Low cotransmission",
    pops$cotransmission_cat == "High_Cotransmission" ~ "High cotransmission",
    TRUE ~ pops$cotransmission_cat
  )
  .bench$pops <- pops
  pops
}

load_sample_truth <- function() {
  if (!is.null(.bench$sample_truth)) return(.bench$sample_truth)
  path <- file.path(benchmark_root(), "simulation_setup", "expected_ecoi.tsv")
  truth <- read_chars(path, c("population_name", "library_sample_name", "COI", "eCOI"))
  truth <- truth %>%
    mutate(
      sample_id = sample_key(library_sample_name),
      COI = as.numeric(COI),
      eCOI = as.numeric(eCOI)
    ) %>%
    inner_join(kept_populations(), by = "population_name")
  .bench$sample_truth <- truth
  truth
}

coi_piece <- function(panel, tool_dir, file, tool, column, label = NULL) {
  raw <- read_chars(
    file.path(panel_dir(panel), tool_dir, file),
    c("population_name", "paramset_id", "specimen_name", column)
  )
  estimate <- as.numeric(raw[[column]])
  if (identical(column, "coi_freq")) estimate[is.na(estimate)] <- 1
  paramset_id <- raw$paramset_id
  if (!is.null(label)) paramset_id <- paste0(paramset_id, " (", label, ")")
  data.frame(
    tool = tool,
    paramset_id = paramset_id,
    population_name = raw$population_name,
    sample_id = sample_key(raw$specimen_name),
    estimate = estimate,
    stringsAsFactors = FALSE
  )
}

load_coi_estimates <- function(panel) {
  key <- paste0("coi_", panel)
  if (!is.null(.bench[[key]])) return(.bench[[key]])
  parts <- list(
    coi_piece(panel, "naive_coi", "coi.tsv.gz", "Naive COI", "coi"),
    coi_piece(panel, "MOIRE", "coi_summary.tsv.gz", "MOIRE", "coi"),
    coi_piece(panel, "THEREALMcCOIL", "coi.tsv.gz", "TheREALMcCOIL", "coi"),
    coi_piece(panel, "coiaf", "coi.tsv.gz", "coiaf", "coi_freq", "frequency score"),
    coi_piece(panel, "coiaf", "coi.tsv.gz", "coiaf", "coi_variant", "variant score"),
    coi_piece(panel, "SNPSlice", "coi.tsv.gz", "SNP-slice", "coi", "coi"),
    coi_piece(panel, "SNPSlice", "coi.tsv.gz", "SNP-slice", "coi_cons_weighted", "consensus weighted")
  )
  kept <- kept_populations()$population_name
  out <- bind_rows(parts) %>%
    filter(population_name %in% kept, is.finite(estimate))
  .bench[[key]] <- out
  out
}

pick_paramset <- function(df, estimate, truth) {
  scored <- df %>%
    group_by(tool, paramset_id) %>%
    summarise(
      n = sum(is.finite(.data[[estimate]]) & is.finite(.data[[truth]])),
      n_pop = n_distinct(population_name),
      ccc = lin_ccc(.data[[estimate]], .data[[truth]]),
      .groups = "drop"
    ) %>%
    filter(n >= 3, is.finite(ccc)) %>%
    group_by(tool) %>%
    filter(n_pop == max(n_pop)) %>%
    arrange(desc(ccc), paramset_id) %>%
    slice(1) %>%
    ungroup()
  if (!nrow(scored)) stop("No parameter set could be scored.")
  scored
}

winner_rows <- function(df, winners, preferred = tool_levels) {
  df %>%
    inner_join(winners %>% select(tool, paramset_id), by = c("tool", "paramset_id")) %>%
    mutate(tool = order_tools(tool, preferred))
}

sample_coi_results <- function(panel) {
  key <- paste0("sample_results_", panel)
  if (!is.null(.bench[[key]])) return(.bench[[key]])
  joined <- load_coi_estimates(panel) %>%
    inner_join(
      load_sample_truth() %>% select(population_name, sample_id, COI),
      by = c("population_name", "sample_id")
    )
  if (!nrow(joined)) stop("Sample COI estimates did not match the simulated samples.")
  scored_on <- joined %>% filter(COI >= 2, COI <= 9)
  winners <- pick_paramset(scored_on, "estimate", "COI")
  chosen <- winner_rows(joined, winners, coi_tool_order)
  metrics <- chosen %>%
    filter(COI >= 2, COI <= 9) %>%
    group_by(tool, paramset_id) %>%
    summarise(
      `Lin's CCC` = lin_ccc(estimate, COI),
      MAE = mean(abs(estimate - COI)),
      `Mean Bias` = mean(estimate - COI),
      Samples = n(),
      .groups = "drop"
    )
  mono <- chosen %>%
    filter(abs(COI - 1) < 1e-6) %>%
    group_by(tool) %>%
    summarise(`Monoclonal recall` = mean(estimate < 1.5), .groups = "drop")
  metrics <- metrics %>%
    left_join(mono, by = "tool") %>%
    arrange(tool) %>%
    transmute(
      Tool = as.character(tool),
      `Best parameter set` = paramset_id,
      `Lin's CCC`, MAE, `Mean Bias`, Samples, `Monoclonal recall`
    )
  ccc_by_pop <- chosen %>%
    filter(COI >= 2, COI <= 9) %>%
    group_by(tool, population_name) %>%
    summarise(ccc = lin_ccc(estimate, COI), .groups = "drop") %>%
    filter(is.finite(ccc))
  mae_by_pop <- chosen %>%
    filter(COI >= 1, COI <= 9) %>%
    mutate(truth = factor(as.integer(round(COI)), levels = 1:9)) %>%
    group_by(tool, population_name, truth) %>%
    summarise(
      mae = mean(abs(estimate - COI)),
      error = mean(estimate - COI),
      .groups = "drop"
    ) %>%
    filter(is.finite(mae), is.finite(error), mae < 6)
  plot_df <- chosen %>%
    filter(COI >= 1, COI <= 9) %>%
    mutate(truth = factor(as.integer(round(COI)), levels = 1:9))
  out <- list(
    metrics = metrics,
    ccc_by_pop = ccc_by_pop,
    mae_by_pop = mae_by_pop,
    plot_df = plot_df,
    panel = panel
  )
  .bench[[key]] <- out
  out
}

plot_sample_coi <- function(results) {
  ref <- data.frame(truth = factor(1:9, levels = 1:9), y = 1:9)
  ggplot(results$plot_df, aes(truth, estimate, fill = tool)) +
    geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.25, linewidth = 0.3) +
    geom_line(
      data = ref,
      aes(truth, y, group = 1),
      inherit.aes = FALSE,
      colour = "red",
      linewidth = 0.4
    ) +
    facet_wrap(~tool, ncol = 3) +
    scale_fill_manual(values = coi_colours, guide = "none") +
    scale_y_continuous(breaks = seq(0, 15, by = 5)) +
    coord_cartesian(ylim = c(0, 15)) +
    labs(x = "Simulated COI", y = "Estimated COI") +
    theme_bench()
}

plot_sample_ccc <- function(results) {
  d <- results$ccc_by_pop %>%
    mutate(tool = order_tools(tool, coi_tool_order))
  ggplot(d, aes(tool, ccc, fill = tool)) +
    geom_boxplot(outlier.size = 0.8, linewidth = 0.3) +
    geom_hline(yintercept = 0.95, linetype = "dotted", colour = "red") +
    scale_fill_manual(values = coi_colours, guide = "none") +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(x = NULL, y = "Lin's CCC") +
    theme_bench() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

plot_sample_mae <- function(results) {
  d <- results$mae_by_pop %>%
    mutate(tool = order_tools(tool, coi_tool_order))
  ggplot(d, aes(truth, mae, fill = tool)) +
    geom_boxplot(outlier.size = 0.6, linewidth = 0.3) +
    scale_fill_manual(values = coi_colours, name = NULL) +
    labs(x = "Simulated COI", y = "Mean absolute error") +
    theme_bench("bottom")
}

plot_sample_error <- function(results) {
  d <- results$mae_by_pop %>%
    mutate(tool = order_tools(tool, coi_tool_order))
  ggplot(d, aes(truth, error, fill = tool)) +
    geom_boxplot(outlier.size = 0.6, linewidth = 0.3) +
    geom_hline(yintercept = 0, colour = "red", linewidth = 0.4) +
    scale_fill_manual(values = coi_colours, name = NULL) +
    labs(x = "Simulated COI", y = "Mean error") +
    theme_bench("bottom")
}

plot_monoclonal <- function(results) {
  d <- results$metrics %>%
    mutate(
      tool = order_tools(Tool, coi_tool_order),
      label = sprintf("%.1f%%", `Monoclonal recall` * 100)
    )
  ggplot(d, aes(`Monoclonal recall`, tool, fill = tool)) +
    geom_col(width = 0.72) +
    geom_label(
      aes(label = label),
      hjust = 1.08,
      fill = "white",
      linewidth = 0,
      size = 3.2
    ) +
    scale_y_discrete(limits = rev(levels(d$tool))) +
    scale_fill_manual(values = coi_colours, guide = "none") +
    scale_x_continuous(
      limits = c(0, 1),
      labels = function(x) paste0(round(100 * x), "%"),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(x = "Monoclonal recall", y = NULL) +
    theme_bench()
}

mean_coi_results <- function(panel) {
  key <- paste0("mean_results_", panel)
  if (!is.null(.bench[[key]])) return(.bench[[key]])
  joined <- load_coi_estimates(panel) %>%
    inner_join(load_sample_truth(), by = c("population_name", "sample_id")) %>%
    filter(is.finite(estimate), is.finite(COI))
  means <- joined %>%
    group_by(tool, paramset_id, population_name, mean_coi, relatedness, cotransmission) %>%
    summarise(
      estimate = mean(estimate),
      truth = mean(COI),
      .groups = "drop"
    )
  winners <- pick_paramset(means, "estimate", "truth")
  chosen <- winner_rows(means, winners, coi_tool_order)
  metrics <- chosen %>%
    group_by(tool, paramset_id) %>%
    summarise(
      `Lin's CCC` = lin_ccc(estimate, truth),
      MAE = mean(abs(estimate - truth)),
      `Mean Bias` = mean(estimate - truth),
      Populations = n(),
      .groups = "drop"
    ) %>%
    arrange(tool) %>%
    transmute(
      Tool = as.character(tool),
      `Best parameter set` = paramset_id,
      `Lin's CCC`, MAE, `Mean Bias`, Populations
    )
  by_mean <- stratum_table(chosen, "mean_coi", "Mean COI") %>%
    left_join(
      chosen %>% distinct(mean_coi, relatedness),
      by = c("Mean COI" = "mean_coi")
    ) %>%
    transmute(
      Tool, `Mean COI`, Relatedness = relatedness,
      `Lin's CCC`, MAE, `Mean Bias`
    )
  by_cot <- stratum_table(chosen, "cotransmission", "Cotransmission")
  out <- list(
    metrics = metrics,
    by_mean = by_mean,
    by_cotransmission = by_cot,
    plot_df = chosen,
    panel = panel
  )
  .bench[[key]] <- out
  out
}

stratum_table <- function(chosen, column, label) {
  chosen %>%
    group_by(tool, .data[[column]]) %>%
    summarise(
      `Lin's CCC` = lin_ccc(estimate, truth),
      MAE = mean(abs(estimate - truth)),
      `Mean Bias` = mean(estimate - truth),
      .groups = "drop"
    ) %>%
    arrange(tool, .data[[column]]) %>%
    transmute(
      Tool = as.character(tool),
      !!label := .data[[column]],
      `Lin's CCC`, MAE, `Mean Bias`
    )
}

mean_coi_levels <- c("Mean COI 1.26", "Mean COI 3.41", "Mean COI 5.47")

mean_coi_plot_df <- function(results) {
  results$plot_df %>%
    mutate(
      tool = order_tools(tool, coi_tool_order),
      mean_coi = factor(mean_coi, levels = mean_coi_levels)
    )
}

plot_mean_coi <- function(results) {
  ggplot(mean_coi_plot_df(results), aes(truth, estimate, colour = tool)) +
    geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.4) +
    geom_point(size = 1.8, alpha = 0.85) +
    facet_wrap(~tool, ncol = 3) +
    scale_colour_manual(values = coi_colours, guide = "none") +
    labs(x = "Simulated mean COI", y = "Estimated mean COI") +
    theme_bench() +
    coord_fixed()
}

mean_coi_level_df <- function(results) {
  mean_coi_plot_df(results) %>%
    mutate(
      level = factor(
        case_when(
          mean_coi == "Mean COI 1.26" ~ "Low",
          mean_coi == "Mean COI 3.41" ~ "Medium",
          mean_coi == "Mean COI 5.47" ~ "High"
        ),
        levels = c("Low", "Medium", "High")
      ),
      abs_err = abs(estimate - truth)
    )
}

plot_mean_by_level <- function(results) {
  d <- mean_coi_level_df(results)
  ref <- d %>%
    group_by(level) %>%
    summarise(y = mean(truth), .groups = "drop") %>%
    mutate(
      x = as.numeric(level) - 0.4,
      xend = as.numeric(level) + 0.4
    )
  ggplot(d, aes(level, estimate, fill = tool)) +
    geom_boxplot(outlier.size = 0.5, linewidth = 0.3) +
    geom_segment(
      data = ref,
      aes(x = x, xend = xend, y = y, yend = y),
      inherit.aes = FALSE,
      colour = "red",
      linewidth = 0.6
    ) +
    scale_fill_manual(values = coi_colours, name = NULL) +
    labs(x = "Simulated mean COI category", y = "Estimated mean COI") +
    theme_bench("bottom")
}

plot_mean_mae_by_level <- function(results) {
  ggplot(mean_coi_level_df(results), aes(level, abs_err, fill = tool)) +
    geom_boxplot(outlier.size = 0.5, linewidth = 0.3) +
    scale_fill_manual(values = coi_colours, name = NULL) +
    labs(x = "Simulated mean COI category", y = "Mean absolute error") +
    theme_bench("bottom")
}

ecoi_results <- function(panel) {
  key <- paste0("ecoi_results_", panel)
  if (!is.null(.bench[[key]])) return(.bench[[key]])
  sample_choice <- sample_coi_results(panel)$metrics
  paramset <- sample_choice$`Best parameter set`[sample_choice$Tool == "MOIRE"]
  if (length(paramset) != 1) stop("MOIRE sample-level parameter set was not found.")
  raw <- read_chars(
    file.path(panel_dir(panel), "MOIRE", "effective_coi_summary.tsv.gz"),
    c("population_name", "paramset_id", "specimen_name", "ecoi")
  )
  estimates <- raw %>%
    filter(paramset_id == paramset, population_name %in% kept_populations()$population_name) %>%
    transmute(
      tool = "MOIRE",
      paramset_id = paramset,
      population_name = population_name,
      sample_id = sample_key(specimen_name),
      estimate = as.numeric(ecoi)
    )
  chosen <- estimates %>%
    inner_join(
      load_sample_truth() %>% select(population_name, sample_id, eCOI),
      by = c("population_name", "sample_id")
    ) %>%
    filter(is.finite(estimate), is.finite(eCOI))
  if (!nrow(chosen)) stop("Effective COI estimates did not match the simulated samples.")
  metrics <- tibble(
    Tool = "MOIRE",
    `Best parameter set` = paramset,
    `Lin's CCC` = lin_ccc(chosen$estimate, chosen$eCOI),
    MAE = mean(abs(chosen$estimate - chosen$eCOI)),
    `Mean Bias` = mean(chosen$estimate - chosen$eCOI),
    Samples = nrow(chosen)
  )
  out <- list(metrics = metrics, plot_df = chosen, panel = panel)
  .bench[[key]] <- out
  out
}

plot_ecoi <- function(results) {
  ggplot(results$plot_df, aes(eCOI, estimate)) +
    geom_bin2d(bins = 40) +
    geom_abline(slope = 1, intercept = 0, colour = "#4a4a4a", linewidth = 0.4) +
    scale_fill_gradient(low = "#f9eaf1", high = "#c82f61", trans = "log10") +
    coord_fixed() +
    labs(x = "Simulated effective COI", y = "MOIRE effective COI", fill = "Samples") +
    theme_bench("right")
}

panel_loci <- function() {
  if (!is.null(.bench$loci)) return(.bench$loci)
  root <- file.path(benchmark_root(), "simulation_setup")
  aa <- read_chars(file.path(root, "aa_variants_positions.tsv"), c("ID", "aaPosition"))
  from_aa <- paste(aa$ID, aa$aaPosition, sep = ":")
  bed <- read_tsv(
    file.path(root, "mad4hatter_panel_info", "amino_acid_codon_positions.bed"),
    col_names = FALSE,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE,
    progress = FALSE
  )
  from_bed <- sub("-AA", ":", bed[[4]], fixed = TRUE)
  loci <- intersect(from_aa, from_bed)
  if (!length(loci)) stop("No single-locus panel positions were found.")
  .bench$loci <- loci
  loci
}

read_population_files <- function(filename, columns) {
  root <- file.path(benchmark_root(), "simulation_setup")
  parts <- lapply(kept_populations()$population_name, function(pop) {
    path <- file.path(root, pop, filename)
    if (!file.exists(path)) return(NULL)
    read_chars(path, columns)
  })
  bind_rows(parts)
}

single_locus_truth <- function() {
  if (!is.null(.bench$sl_truth)) return(.bench$sl_truth)
  raw <- read_population_files(
    "population_aa_variants_freqs_prevs.tsv.gz",
    c("population_name", "transcript_id", "aa_position", "aa", "freq")
  )
  loci <- panel_loci()
  truth <- raw %>%
    mutate(
      locus = paste(transcript_id, aa_position, sep = ":"),
      allele_id = paste(transcript_id, aa_position, aa, sep = ":"),
      expected = as.numeric(freq)
    ) %>%
    filter(locus %in% loci, is.finite(expected)) %>%
    group_by(population_name, allele_id) %>%
    summarise(expected = mean(expected), .groups = "drop")
  .bench$sl_truth <- truth
  truth
}

microhap_truth <- function() {
  if (!is.null(.bench$mh_truth)) return(.bench$mh_truth)
  raw <- read_population_files(
    "population_trimmed_mhap_freqs_prevs.tsv.gz",
    c("population_name", "target", "seq", "freq")
  )
  truth <- raw %>%
    mutate(
      allele_id = paste(trimws(target), trimws(seq), sep = "|"),
      expected = as.numeric(freq)
    ) %>%
    filter(is.finite(expected)) %>%
    group_by(population_name, allele_id) %>%
    summarise(expected = mean(expected), .groups = "drop")
  .bench$mh_truth <- truth
  truth
}

aa_count <- function(variant) {
  loci <- strsplit(variant, ";", fixed = TRUE)
  vapply(loci, function(parts) {
    sum(vapply(parts, function(part) {
      bits <- strsplit(part, ":", fixed = TRUE)[[1]]
      if (length(bits) < 2) return(0L)
      nchar(gsub("[^_]", "", bits[[2]], perl = FALSE)) + 1L
    }, integer(1)))
  }, integer(1))
}

multi_locus_truth <- function(kind) {
  key <- paste0("ml_truth_", kind)
  if (!is.null(.bench[[key]])) return(.bench[[key]])
  file <- if (kind == "mixed") {
    "population_aa_multiloci_variants_freqs_prevs_mixed.tsv.gz"
  } else {
    "population_aa_multiloci_variants_freqs_prevs.tsv.gz"
  }
  raw <- read_population_files(
    file,
    c("population_name", "variant", "aas_in_variant", "freq", "group_id")
  )
  groups <- multi_locus_groups(kind)
  truth <- raw %>%
    mutate(
      aas = as.integer(aas_in_variant),
      expected = as.numeric(freq)
    ) %>%
    filter(group_id %in% groups, is.finite(expected), is.finite(aas))
  full_aas <- truth %>%
    group_by(group_id) %>%
    summarise(full_aas = max(aas), .groups = "drop")
  truth <- truth %>%
    inner_join(full_aas, by = "group_id") %>%
    filter(aas == full_aas) %>%
    mutate(allele_id = paste(group_id, variant, sep = "|")) %>%
    group_by(population_name, allele_id, group_id) %>%
    summarise(expected = mean(expected), .groups = "drop")
  .bench[[paste0("ml_aas_", kind)]] <- full_aas
  .bench[[key]] <- truth
  truth
}

multi_locus_groups <- function(kind) {
  if (kind == "drug") {
    return(c("pfdhfr", "pfdhps", "pfdhfr_pfdhps_quintuple", "pfdhfr_pfdhps_octuple"))
  }
  if (kind == "mixed") {
    raw <- read_population_files(
      "population_aa_multiloci_variants_freqs_prevs_mixed.tsv.gz",
      "group_id"
    )
    return(unique(raw$group_id[grepl("mixed", raw$group_id, fixed = TRUE)]))
  }
  path <- file.path(
    benchmark_root(), "simulation_setup", "mad4hatter_panel_info",
    "biallelic_only_loci_groups_aa_level.tsv"
  )
  groups <- read_chars(path, "group_id")$group_id
  unique(groups)
}

allele_estimates <- function(panel, tool_dir, file, tool, id_columns) {
  columns <- unique(c("population_name", "paramset_id", "freq", id_columns))
  raw <- read_chars(file.path(panel_dir(panel), tool_dir, file), columns)
  raw$freq <- as.numeric(raw$freq)
  raw <- raw[is.finite(raw$freq) & raw$population_name %in% kept_populations()$population_name, , drop = FALSE]
  if ("target_name" %in% names(raw)) {
    allele_id <- paste(trimws(raw$target_name), trimws(raw$seq), sep = "|")
  } else if ("group_id" %in% names(raw)) {
    allele_id <- paste(raw$group_id, raw$variant, sep = "|")
  } else {
    allele_id <- raw$variant
  }
  data.frame(
    tool = tool,
    paramset_id = raw$paramset_id,
    population_name = raw$population_name,
    allele_id = allele_id,
    group_id = if ("group_id" %in% names(raw)) raw$group_id else NA_character_,
    variant = if ("variant" %in% names(raw)) raw$variant else NA_character_,
    estimate = raw$freq,
    stringsAsFactors = FALSE
  )
}

prepare_frequency <- function(estimates, truth, preferred = tool_levels, apply_rule = TRUE) {
  coverage <- estimates %>%
    distinct(tool, paramset_id, population_name)
  expanded <- coverage %>%
    inner_join(truth, by = "population_name", relationship = "many-to-many") %>%
    left_join(
      estimates %>% select(tool, paramset_id, population_name, allele_id, estimate),
      by = c("tool", "paramset_id", "population_name", "allele_id")
    )
  if (apply_rule) {
    rule <- apply_frequency_rule(expanded$expected, expanded$estimate)
    expanded$estimate <- rule$estimate
    scored <- expanded[rule$keep, , drop = FALSE]
  } else {
    scored <- expanded[is.finite(expanded$expected) & is.finite(expanded$estimate), , drop = FALSE]
  }
  n_pop <- coverage %>% count(tool, paramset_id, name = "n_pop")
  winners <- scored %>%
    group_by(tool, paramset_id) %>%
    summarise(
      n = n(),
      ccc = lin_ccc(estimate, expected),
      .groups = "drop"
    ) %>%
    inner_join(n_pop, by = c("tool", "paramset_id")) %>%
    filter(n >= 3, is.finite(ccc)) %>%
    group_by(tool) %>%
    filter(n_pop == max(n_pop)) %>%
    arrange(desc(ccc), paramset_id) %>%
    slice(1) %>%
    ungroup()
  if (!nrow(winners)) stop("No frequency parameter set could be scored.")
  chosen <- winner_rows(scored, winners, preferred) %>%
    left_join(
      kept_populations() %>% select(population_name, mean_coi),
      by = "population_name"
    )
  metrics <- chosen %>%
    group_by(tool, paramset_id) %>%
    group_modify(~ frequency_metrics(.x)) %>%
    ungroup() %>%
    arrange(tool) %>%
    transmute(
      Tool = as.character(tool),
      `Best parameter set` = paramset_id,
      `Lin's CCC`, MAE, Detection, Alleles
    )
  by_freq <- chosen %>%
    mutate(bin = frequency_bin(expected)) %>%
    group_by(tool, bin) %>%
    group_modify(~ frequency_metrics(.x)) %>%
    ungroup() %>%
    arrange(tool, bin) %>%
    transmute(
      Tool = as.character(tool),
      `Simulated frequency` = as.character(bin),
      `Lin's CCC`, MAE, Detection
    )
  by_mean <- chosen %>%
    group_by(tool, mean_coi) %>%
    group_modify(~ frequency_metrics(.x)) %>%
    ungroup() %>%
    arrange(tool, mean_coi) %>%
    transmute(
      Tool = as.character(tool),
      `Mean COI` = mean_coi,
      `Lin's CCC`, MAE, Detection
    )
  list(metrics = metrics, by_freq = by_freq, by_mean = by_mean, plot_df = chosen)
}

apply_frequency_rule <- function(expected, estimate) {
  est_below <- below_one_percent(estimate)
  exp_below <- is.finite(expected) & expected < (0.01 - 1e-8)
  exp_above <- is.finite(expected) & expected > (0.01 + 1e-8)
  at_one <- is.finite(expected) & !exp_below & !exp_above
  estimate[exp_above & est_below] <- 0
  estimate[at_one & !is.finite(estimate)] <- 0
  keep <- !(exp_below & est_below) & is.finite(expected) & is.finite(estimate)
  list(estimate = estimate, keep = keep)
}

frequency_bin <- function(expected) {
  cut(
    expected,
    breaks = c(-Inf, 0.01, 0.05, 0.1, 0.25, Inf),
    labels = c("Below 1%", "1-5%", "5-10%", "10-25%", "25% and above"),
    right = FALSE
  )
}

frequency_metrics <- function(df) {
  detectable <- is.finite(df$expected) & df$expected >= (0.01 - 1e-8)
  detection <- if (any(detectable)) mean(df$estimate[detectable] > 0) else NA_real_
  tibble(
    `Lin's CCC` = lin_ccc(df$estimate, df$expected),
    MAE = mean(abs(df$estimate - df$expected)),
    Detection = detection,
    Alleles = nrow(df)
  )
}

plot_frequency_mae <- function(results, colours, preferred) {
  d <- results$plot_df %>%
    filter(is.finite(expected), expected >= 0, expected <= 1 + 1e-8) %>%
    mutate(
      expected = pmin(expected, 1),
      bin = cut(
        expected,
        breaks = seq(0, 1, by = 0.2),
        include.lowest = TRUE,
        right = FALSE,
        labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1")
      ),
      tool = order_tools(tool, preferred)
    ) %>%
    filter(!is.na(bin)) %>%
    group_by(tool, population_name, bin) %>%
    summarise(mae = mean(abs(estimate - expected)), .groups = "drop") %>%
    filter(is.finite(mae))
  ggplot(d, aes(bin, mae, fill = tool)) +
    geom_boxplot(outlier.size = 0.5, linewidth = 0.3) +
    scale_fill_manual(values = colours, name = NULL) +
    labs(x = "Simulated frequency", y = "Mean absolute error") +
    theme_bench("bottom")
}

presence_calls <- function(estimates, truth, metrics) {
  winners <- metrics %>%
    transmute(tool = as.character(Tool), paramset_id = as.character(`Best parameter set`))
  chosen <- estimates %>%
    inner_join(winners, by = c("tool", "paramset_id"))
  coverage <- chosen %>% distinct(tool, paramset_id, population_name)
  from_truth <- coverage %>%
    inner_join(truth, by = "population_name", relationship = "many-to-many") %>%
    left_join(
      chosen %>% select(tool, paramset_id, population_name, allele_id, estimate),
      by = c("tool", "paramset_id", "population_name", "allele_id")
    ) %>%
    select(tool, group_id, estimate, expected)
  from_tool <- chosen %>%
    anti_join(truth %>% select(population_name, allele_id), by = c("population_name", "allele_id")) %>%
    transmute(tool, group_id, estimate, expected = 0)
  bind_rows(from_truth, from_tool) %>%
    mutate(
      estimate = ifelse(is.finite(estimate), estimate, 0),
      expected = ifelse(is.finite(expected), expected, 0)
    )
}

precision_recall <- function(calls) {
  calls %>%
    group_by(tool, group_id) %>%
    summarise(
      tp = sum(estimate != 0 & expected != 0),
      fp = sum(estimate != 0 & expected == 0),
      fn = sum(estimate == 0 & expected != 0),
      .groups = "drop"
    ) %>%
    mutate(
      Precision = ifelse(tp + fp > 0, tp / (tp + fp), NA_real_),
      Recall = ifelse(tp + fn > 0, tp / (tp + fn), NA_real_)
    )
}

filter_presence_calls <- function(calls) {
  rule <- apply_frequency_rule(calls$expected, calls$estimate)
  calls$estimate <- rule$estimate
  calls[rule$keep, , drop = FALSE]
}

multi_group_levels <- c(
  "2_low_div", "2_med_div", "4_low_div", "4_med_div", "6_low_div", "6_med_div",
  "2_low_div_mixed", "2_med_div_mixed", "2_high_div_mixed",
  "4_low_div_mixed", "4_med_div_mixed", "4_high_div_mixed",
  "6_low_div_mixed", "6_med_div_mixed", "6_high_div_mixed",
  "pfdhfr", "pfdhps", "pfdhfr_pfdhps_quintuple", "pfdhfr_pfdhps_octuple"
)

plot_precision_recall <- function(pr) {
  present <- unique(as.character(pr$group_id))
  levels <- c(multi_group_levels[multi_group_levels %in% present], setdiff(present, multi_group_levels))
  d <- bind_rows(
    pr %>% transmute(tool, group_id, measure = "Precision", value = Precision),
    pr %>% transmute(tool, group_id, measure = "Recall", value = Recall)
  ) %>%
    mutate(
      tool = order_tools(tool, multi_tool_order),
      group_id = factor(group_id, levels = levels),
      measure = factor(measure, levels = c("Precision", "Recall"))
    )
  ggplot(d, aes(group_id, value, fill = tool)) +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.72,
      colour = "black",
      linewidth = 0.2
    ) +
    facet_wrap(~measure) +
    scale_fill_manual(values = multi_colours, name = NULL) +
    scale_x_discrete(expand = expansion(add = 0.6)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    labs(x = NULL, y = NULL) +
    theme_bench("bottom") +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      plot.margin = margin(6, 10, 6, 22)
    )
}

plot_frequency <- function(results) {
  ggplot(results$plot_df, aes(expected, estimate)) +
    geom_bin2d(bins = 40) +
    geom_abline(slope = 1, intercept = 0, colour = "#4a4a4a", linewidth = 0.3) +
    scale_fill_gradient(low = "#f4f7f9", high = "#325b74", trans = "log10") +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    facet_wrap(~tool, ncol = 2) +
    labs(x = "Simulated frequency", y = "Estimated frequency", fill = "Alleles") +
    theme_bench("right")
}

single_variant_results <- function() {
  if (!is.null(.bench$sl_results)) return(.bench$sl_results)
  estimates <- bind_rows(
    allele_estimates("full", "naive_slaf", "slaf.tsv.gz", "Naive single-locus allele frequency", "variant"),
    allele_estimates("full", "IDM", "IDM_summary.tsv.gz", "Incomplete data model", "variant"),
    allele_estimates("full", "dcifer_slaf_from_mhap_freqs", "slaf.tsv.gz", "dcifer", "variant")
  )
  out <- prepare_frequency(
    estimates, single_locus_truth(), variant_tool_order, apply_rule = FALSE
  )
  .bench$sl_results <- out
  out
}

microhap_results <- function() {
  if (!is.null(.bench$mh_results)) return(.bench$mh_results)
  estimates <- bind_rows(
    allele_estimates(
      "full", "MOIRE", "allele_freq_summary.tsv.gz", "MOIRE",
      c("target_name", "seq")
    ),
    allele_estimates(
      "full", "dcifer_mhap_freqs", "mhaps_freqs.tsv.gz", "dcifer",
      c("target_name", "seq")
    ),
    allele_estimates(
      "full", "naive_mhap_freqs", "mhaps_freqs.tsv.gz", "Naive microhaplotype frequency",
      c("target_name", "seq")
    )
  )
  out <- prepare_frequency(
    estimates, microhap_truth(), mhap_tool_order, apply_rule = FALSE
  )
  .bench$mh_results <- out
  out
}

multi_locus_results <- function(kind) {
  if (!kind %in% c("biallelic", "drug", "mixed")) stop("Unknown multi-locus comparison: ", kind)
  key <- paste0("ml_results_", kind)
  if (!is.null(.bench[[key]])) return(.bench[[key]])
  groups <- multi_locus_groups(kind)
  pieces <- list(
    allele_estimates("full", "naive_mlafp", "mlafp.tsv.gz", "Naive multi-locus frequency", c("variant", "group_id")),
    allele_estimates("full", "FEM", "FEM_summary.tsv.gz", "FreqEstimationModel", c("variant", "group_id")),
    allele_estimates("full", "MLBM", "MLBM_summary.tsv.gz", "MultiLociBiallelicModel", c("variant", "group_id")),
    allele_estimates("full", "SNPSlice", "mlaf.tsv.gz", "SNP-slice", c("variant", "group_id"))
  )
  if (kind == "mixed") {
    pieces <- lapply(pieces, function(d) d[grepl("mixed", d$paramset_id, fixed = TRUE), , drop = FALSE])
  } else {
    pieces <- lapply(pieces, function(d) d[!grepl("mixed", d$paramset_id, fixed = TRUE), , drop = FALSE])
  }
  estimates <- bind_rows(pieces) %>%
    filter(group_id %in% groups)
  if (kind %in% c("drug", "mixed")) {
    estimates <- estimates %>%
      filter(tool %in% c("Naive multi-locus frequency", "SNP-slice"))
  }
  aas <- multi_locus_truth(kind)
  full_aas <- .bench[[paste0("ml_aas_", kind)]]
  estimates <- estimates %>%
    mutate(aa_n = aa_count(variant)) %>%
    inner_join(full_aas, by = "group_id") %>%
    filter(aa_n == full_aas) %>%
    select(-aa_n, -full_aas)
  if (!nrow(estimates)) stop("No multi-locus estimates were found for ", kind, ".")
  out <- prepare_frequency(estimates, aas, multi_tool_order)
  calls <- presence_calls(estimates, aas, out$metrics)
  out$precision_recall <- precision_recall(calls)
  out$precision_recall_filtered <- precision_recall(filter_presence_calls(calls))
  .bench[[key]] <- out
  out
}
