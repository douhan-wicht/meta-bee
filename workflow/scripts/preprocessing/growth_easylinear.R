# ============================================================
#  growthrates_easylinear.R
#  Estimate growth parameters using growthrates::fit_easylinear()
#  from an aggregated growth-data TSV file.
#
#  CLI-driven for Snakemake:
#    --input-tsv   path/to/growth-data.tsv
#    --output-csv  path/to/growth_rates.csv
#    --plot-dir    path/to/plots_root
#    [--time-unit hours]
#    [--h-window  8]
#    [--quota     0.95]
# ============================================================

# --- Install required packages ---
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(ggplot2)
})

# ---- Ensure growthrates is installed ----
if (!requireNamespace("growthrates", quietly = TRUE)) {
  cat("Package 'growthrates' not found. Installing from CRAN...\n")
  install.packages("growthrates", repos = "https://cloud.r-project.org")
}
library(growthrates)

# ---- parse command line options ----
option_list <- list(
  make_option(c("-i", "--input-tsv"), type = "character", help = "Input TSV with columns species, replicate, media, cond, time, value.", metavar = "file"),
  make_option(c("-o", "--output-csv"), type = "character", help = "Output CSV with growth parameters.", metavar = "file"),
  make_option(c("-p", "--plot-dir"), type = "character", help = "Output directory for plots.", metavar = "dir"),
  make_option(c("-t", "--time-unit"), type = "character", default = "hours",
              help = "Time unit for mumax (default: %default). If 'hours', time is converted from minutes to hours."),
  make_option(c("--h-window"), type = "integer", default = 8,
              help = "Moving window size h for easylinear (default: %default)."),
  make_option(c("--quota"), type = "double", default = 0.95,
              help = "Quota parameter for easylinear (default: %default).")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$`input-tsv`) || is.null(opt$`output-csv`) || is.null(opt$`plot-dir`)) {
  stop("You must supply --input-tsv, --output-csv and --plot-dir", call. = FALSE)
}

input_file  <- opt$`input-tsv`
output_file <- opt$`output-csv`
plot_root   <- opt$`plot-dir`
time_unit   <- opt$`time-unit`
h_window    <- opt$`h-window`
quota_val   <- opt$quota

cat("Input TSV   :", input_file, "\n")
cat("Output CSV  :", output_file, "\n")
cat("Plot dir    :", plot_root, "\n")
cat("Time unit   :", time_unit, "\n")
cat("h_window    :", h_window, "\n")
cat("quota_val   :", quota_val, "\n\n")

dir.create(plot_root, showWarnings = FALSE, recursive = TRUE)

# ---- read data ----
df <- readr::read_tsv(input_file, show_col_types = FALSE)

# --- check that the tsv follows the required GrowthRates guidelines ---
required_cols <- c("species", "replicate", "media", "cond", "time", "value")
if (!all(required_cols %in% names(df))) {
  stop("Input TSV must contain columns: ", paste(required_cols, collapse = ", "))
}

# ---- convert time to hours if needed ----
if (tolower(time_unit) %in% c("hour", "hours", "h")) {
  df <- df %>% mutate(time = time / 60)
}

# ---- clean data ----
df <- df %>%
  filter(is.finite(time), is.finite(value), value > 0) %>%
  arrange(species, media, cond, replicate, time)

# ---- fitting helper ----
fit_one <- function(dat) {
  dat <- dat %>% distinct(time, .keep_all = TRUE) %>% arrange(time)
  if (nrow(dat) < 5) return(NULL)
  tryCatch({
    fit <- growthrates::fit_easylinear(dat$time, dat$value, h = h_window, quota = quota_val)
    coefs <- coef(fit)
    tibble(
      mumax = coefs[["mumax"]],
      lag   = coefs[["lag"]],
      y0    = coefs[["y0"]],
      y0_lm = coefs[["y0_lm"]],
      r2    = growthrates::rsquared(fit)[["r2"]],
      rss   = deviance(fit)
    )
  }, error = function(e) {
    warning(sprintf("Fit failed: %s", conditionMessage(e)))
    tibble(mumax = NA_real_, lag = NA_real_, y0 = NA_real_,
           y0_lm = NA_real_, r2 = NA_real_, rss = NA_real_)
  })
}

# ---- run fits ----
cat("Running fits per (species, media, cond, replicate)...\n")
results <- df %>%
  group_by(species, media, cond, replicate) %>%
  group_modify(~ fit_one(.x)) %>%
  ungroup() %>%
  mutate(time_unit = time_unit)

cat("Preview of results:\n")
print(head(results, 10))

# ---- write output CSV ----
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
readr::write_csv(results, output_file)
cat("\nWrote results to:", output_file, "\n")

# ---- optional: quick aggregate summary ----
summary_table <- results %>%
  group_by(species, media, cond) %>%
  summarise(mean_mumax = mean(mumax, na.rm = TRUE),
            sd_mumax   = sd(mumax, na.rm = TRUE),
            n_reps     = sum(!is.na(mumax)),
            .groups = "drop")

cat("\nMean growth rates (mumax) summary:\n")
print(summary_table)

# ---- plotting: per-species overview + per-combination plots ----
cat("\nCreating growth-curve plots (one per species + one per combination)...\n")

df <- df %>% mutate(replicate = as.factor(replicate))

species_list <- unique(df$species)

for (sp in species_list) {
  cat("  Species:", sp, "\n")
  df_sp <- df %>% filter(species == sp)

  # ---------- 1) Overview plot per species ----------
  df_sp <- df_sp %>% mutate(panel = interaction(media, cond, drop = TRUE))
  n_panels <- dplyr::n_distinct(df_sp$panel)
  facet_rows <- if (n_panels > 5) 2 else 1

  p_overview <- ggplot(df_sp,
                       aes(x = time,
                           y = value,
                           group = replicate,
                           colour = replicate)) +
    geom_line(alpha = 0.8, linewidth = 0.4) +
    facet_wrap(~ media + cond,
               scales = "free_y",
               nrow = facet_rows) +
    labs(
      title  = paste("Growth curves for species:", sp),
      x      = paste("Time [", time_unit, "]"),
      y      = "Measurement (value)",
      colour = "Replicate"
    ) +
    theme_bw() +
    theme(
      strip.background = element_rect(colour = NA, fill = "grey90"),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      panel.spacing    = grid::unit(0.6, "lines")
    )

  # Create species folder (if not already created)
  species_dir <- file.path(plot_root, sp)
  dir.create(species_dir, showWarnings = FALSE, recursive = TRUE)

  # Save overview inside the species folder
  overview_file <- file.path(species_dir, paste0(sp, "_overview.png"))
  ggsave(overview_file, p_overview, width = 14, height = 7, dpi = 300)
  cat("    Saved overview:", overview_file, "\n")

  # ---------- 2) Detailed plots per (media, cond) combination ----------
  species_dir <- file.path(plot_root, sp)
  dir.create(species_dir, showWarnings = FALSE, recursive = TRUE)

  combs <- df_sp %>%
    distinct(media, cond) %>%
    arrange(media, cond)

  purrr::pwalk(combs, function(media, cond) {
    df_comb <- df_sp %>% filter(media == !!media, cond == !!cond)

    p_comb <- ggplot(df_comb,
                     aes(x = time,
                         y = value,
                         group = replicate,
                         colour = replicate)) +
      geom_line(alpha = 0.8, linewidth = 0.7) +
      labs(
        title  = paste("Growth curves:", sp, "|", media, "|", cond),
        x      = paste("Time [", time_unit, "]"),
        y      = "Measurement (value)",
        colour = "Replicate"
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )

    base_name <- paste("growth", sp, media, cond, sep = "_")
    base_name <- gsub("[^[:alnum:]_]+", "_", base_name)
    out_file <- file.path(species_dir, paste0(base_name, ".png"))

    ggsave(out_file, p_comb, width = 7, height = 5, dpi = 300)
    cat("      Saved combo plot:", out_file, "\n")
  })
}

cat("\nAll done.\n")
