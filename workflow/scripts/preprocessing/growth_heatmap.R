#!/usr/bin/env Rscript

# ============================================================
#  growth_heatmap_species_media.R
#
#  Read growth_binary.csv (species, strains, replicate, cond,
#  media, growth) and create a species x media heatmap where
#  each tile is split diagonally:
#
#    - bottom-left  triangle = anaerobic
#    - top-right    triangle = microaerobic
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(purrr)
})

# ---- CLI options ----
option_list <- list(
  make_option(c("-i", "--input"), type = "character", metavar = "file",
              help = "Input growth_binary.csv"),
  make_option(c("--output-png"), type = "character", metavar = "file",
              help = "Output PNG file"),
  make_option(c("--output-pdf"), type = "character", metavar = "file",
              help = "Optional output PDF file"),
  make_option(c("--species-order"), type = "character", default = NULL,
              help = "Optional comma-separated order of species on the y-axis"),
  make_option(c("--media-order"), type = "character", default = NULL,
              help = "Optional comma-separated order of media on the x-axis")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || is.null(opt$`output-png`)) {
  stop("You must supply --input and --output-png (and optionally --output-pdf).",
       call. = FALSE)
}

input_file        <- opt$input
output_png        <- opt$`output-png`
output_pdf        <- opt$`output-pdf`
species_order_arg <- opt$`species-order`
media_order_arg   <- opt$`media-order`

cat("Input      :", input_file, "\n")
cat("Output PNG :", output_png, "\n")
if (!is.null(output_pdf)) cat("Output PDF :", output_pdf, "\n")

# ---- read data ----
dat <- readr::read_csv(input_file, show_col_types = FALSE)

required_cols <- c("species", "media", "cond", "growth")
missing_cols  <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Input file is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

dat <- dat %>% mutate(growth = as.numeric(growth))

# ---- check conditions ----
cond_levels <- sort(unique(dat$cond))
cat("Conditions found:", paste(cond_levels, collapse = ", "), "\n")

expected_conds <- c("anaerobic", "microaerobic")
if (!all(expected_conds %in% cond_levels)) {
  stop(
    "This script expects conditions 'anaerobic' and 'microaerobic'.\n",
    "Found: ", paste(cond_levels, collapse = ", ")
  )
}

# ---- aggregate to species x media x cond ----
agg <- dat %>%
  group_by(species, media, cond) %>%
  summarise(
    any_non_na = any(!is.na(growth)),
    any_growth = any(growth == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    growth_sm = case_when(
      !any_non_na        ~ NA_real_,
      any_growth         ~ 1.0,
      TRUE               ~ 0.0
    )
  ) %>%
  select(species, media, cond, growth_sm)

# ---- build full grid species x media x cond ----
species_all <- sort(unique(dat$species))
media_all   <- sort(unique(dat$media))

grid <- expand.grid(
  species = species_all,
  media   = media_all,
  cond    = expected_conds,
  stringsAsFactors = FALSE
)

heat_df <- grid %>%
  left_join(agg, by = c("species", "media", "cond")) %>%
  mutate(
    growth_factor = case_when(
      is.na(growth_sm) ~ NA_character_,
      growth_sm == 1   ~ "growth",
      growth_sm == 0   ~ "no_growth",
      TRUE             ~ NA_character_
    )
  )

# ---- ordering ----
if (!is.null(species_order_arg)) {
  sp_order <- strsplit(species_order_arg, ",")[[1]] |> trimws()
  heat_df$species <- factor(heat_df$species, levels = sp_order)
} else {
  heat_df$species <- factor(heat_df$species, levels = sort(unique(heat_df$species)))
}

if (!is.null(media_order_arg)) {
  med_order <- strsplit(media_order_arg, ",")[[1]] |> trimws()
  heat_df$media <- factor(heat_df$media, levels = med_order)
} else {
  heat_df$media <- factor(heat_df$media, levels = sort(unique(heat_df$media)))
}

heat_df <- heat_df %>%
  mutate(
    x_center = as.numeric(media),
    y_center = as.numeric(species)
  )

# ---- triangles ----
df_ana <- heat_df %>%
  filter(cond == "anaerobic") %>%
  mutate(
    x = map2(x_center, y_center, ~ c(.x - 0.5, .x + 0.5, .x - 0.5)),
    y = map2(x_center, y_center, ~ c(.y - 0.5, .y - 0.5, .y + 0.5))
  ) %>%
  unnest(c(x, y))

df_micro <- heat_df %>%
  filter(cond == "microaerobic") %>%
  mutate(
    x = map2(x_center, y_center, ~ c(.x + 0.5, .x + 0.5, .x - 0.5)),
    y = map2(x_center, y_center, ~ c(.y + 0.5, .y - 0.5, .y + 0.5))
  ) %>%
  unnest(c(x, y))

tri_df <- bind_rows(df_ana, df_micro)

# ---- plot ----
p <- ggplot(tri_df, aes(x = x, y = y)) +
  geom_polygon(
    aes(
      group = interaction(species, media, cond),
      fill  = growth_factor
    ),
    color = "white"
  ) +
  scale_fill_manual(
    values = c(
      "growth"    = "green",
      "no_growth" = "red"
    ),
    na.value = "grey80",
    name = "Growth"
  ) +
  scale_x_continuous(
    breaks = seq_along(levels(heat_df$media)),
    labels = levels(heat_df$media),
    expand = expansion(mult = 0.02)
  ) +
  scale_y_continuous(
    breaks = seq_along(levels(heat_df$species)),
    labels = levels(heat_df$species),
    expand = expansion(mult = 0.02)
  ) +
  coord_fixed() +
  labs(
    title = "Growth presence/absence by species and media\n(Anaerobic = bottom-left, Microaerobic = top-right)",
    x = "Media",
    y = "Species"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid  = element_blank()
  )

# ---- save ----
dir.create(dirname(output_png), showWarnings = FALSE, recursive = TRUE)
ggsave(output_png, p, width = 8, height = 6, dpi = 300)

if (!is.null(output_pdf)) {
  dir.create(dirname(output_pdf), showWarnings = FALSE, recursive = TRUE)
  ggsave(output_pdf, p, width = 8, height = 6)
}

cat("Saved heatmap to:", output_png, "\n")
if (!is.null(output_pdf)) cat("Saved heatmap PDF to:", output_pdf, "\n")
cat("Done.\n")
