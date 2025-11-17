#!/usr/bin/env Rscript

# ============================================================
#  growth_present_from_easylinear.R
#
#  Read growth_rates.csv and classify growth/no-growth
#  per (species, [strain], cond, media).
#
#  Output columns:
#    species, strains, replicate, cond, media, growth
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(tibble)
  library(rlang)
})

# ---- CLI options ----
option_list <- list(
  make_option(c("-i", "--input"),  type = "character", metavar = "file",
              help = "Input growth_rates.csv"),
  make_option(c("-o", "--output"), type = "character", metavar = "file",
              help = "Output growth_binary.csv"),
  make_option(c("--r2-min"),       type = "double",   default = 0.95,
              help  = "Minimum R^2 for a replicate to be called growth (default: %default)")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || is.null(opt$output)) {
  stop("You must supply --input and --output", call. = FALSE)
}

input_file  <- opt$input
output_file <- opt$output
r2_min      <- opt$`r2-min`

cat("Input : ", input_file,  "\n")
cat("Output: ", output_file, "\n")
cat("R2 min: ", r2_min,      "\n\n")

# ---- read data ----
df <- readr::read_csv(input_file, show_col_types = FALSE)

required_cols <- c("species", "media", "cond", "replicate", "mumax", "r2")
missing_cols  <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Input file is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

# Optional strain / strains column
has_strain  <- "strain"  %in% names(df)
has_strains <- "strains" %in% names(df)

if (has_strain && has_strains) {
  stop("Input has both 'strain' and 'strains' columns. Please keep only one.")
}

if (has_strain) {
  df <- df %>% rename(strain = .data$strain)
} else if (has_strains) {
  df <- df %>% rename(strain = .data$strains)
} else {
  df <- df %>% mutate(strain = NA_character_)
}

# ---- classify growth per replicate ----
df_rep <- df %>%
  mutate(
    growth_rep = case_when(
      is.na(mumax) ~ 0L,  # no positive segment found -> no growth
      is.na(r2)   ~ if_else(mumax > 0, 1L, 0L),
      r2 < r2_min ~ 0L,
      mumax <= 0  ~ 0L,
      TRUE        ~ 1L
    )
  )

# ---- aggregate to species / strain / cond / media level ----
group_vars <- c("species", "strain", "cond", "media")

binary <- df_rep %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(
    replicate = n_distinct(replicate),
    growth    = if_else(any(growth_rep == 1L), 1, 0),
    .groups   = "drop"
  )

binary_out <- binary %>%
  select(species,
         strains = strain,
         replicate,
         cond,
         media,
         growth) %>%
  arrange(species, strains, media, cond)

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
readr::write_csv(binary_out, output_file)

cat("Wrote binary growth calls to:", output_file, "\n")
