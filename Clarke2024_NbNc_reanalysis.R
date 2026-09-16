# ==============================================================================
# Re-analysis of Clarke et al. (2024): co-occurrence of low Nb/Nc and small Nb
# ==============================================================================
#
# Purpose
# -------
# This script re-analyzes the non-replicated dataset associated with
# Clarke et al. (2024) to examine whether low effective-to-census size ratios
# (Nb/Nc) tend to coincide with small absolute effective numbers of breeders
# (Nb), and to summarize the taxonomic composition of the available records.
#
# Data source
# -----------
# The input file, "MEC-23-0962_NonReplicated.csv", is available from the Dryad
# repository associated with Clarke et al. (2024):
#
#   https://datadryad.org/dataset/doi:10.5061/dryad.p2ngf1vzm
#   Dryad DOI: https://doi.org/10.5061/dryad.p2ngf1vzm
#
# The non-replicated dataset removes spatial and temporal replication, so that
# repeatedly sampled populations do not receive disproportionate weight in the
# descriptive summaries below.
#
# Reference
# ---------
# Clarke et al. (2024). Global assessment of effective population sizes:
# Consistent taxonomic differences in meeting the 50/500 rule.
# Molecular Ecology. https://doi.org/10.1111/mec.17353
#
# Notes on interpretation
# -----------------------
# 1. Only records with both Nb and an exact census-size estimate (Nc) are used.
#    Records in which Nc is missing or reported only as a lower bound (e.g.,
#    ">5000" or ">6500") are excluded because an exact Nb/Nc ratio cannot be
#    calculated.
# 2. Nb/Nc is recalculated directly from the reported Nb and Nc values rather
#    than relying on the ratio field supplied in the Dryad dataset.
# 3. The empirical Nb estimates are not standardized to female-only Nb, whereas
#    the theoretical framework in the associated manuscript defines Nb as the
#    effective number of female breeders. Likewise, Nc in the comparative
#    dataset is not necessarily equivalent to the number of mature females in
#    the theoretical model. This analysis is therefore descriptive rather than
#    a direct quantitative validation of the theoretical model.
# 4. The thresholds Nb = 50, 100, and 200 are descriptive reference values and
#    should not be interpreted as universal biological thresholds.
#
# Requirements
# ------------
# R package: dplyr
#
# Expected input file
# -------------------
# MEC-23-0962_NonReplicated.csv
#
# ==============================================================================

library(dplyr)

# ==============================================================================
# 1. Load and clean the non-replicated dataset
# ==============================================================================

input_file <- "MEC-23-0962_NonReplicated.csv"

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Input file not found: ", input_file, "\n",
      "Download it from: ",
      "https://datadryad.org/dataset/doi:10.5061/dryad.p2ngf1vzm"
    )
  )
}

dat <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Remove columns with missing or empty names (the source CSV contains one such
# column in some versions of the downloaded file).
dat <- dat[, !is.na(names(dat)) & names(dat) != ""]

# Confirm that the variables required for this analysis are present.
required_columns <- c(
  "Genus/Species", "Class", "Population", "PopulationID",
  "Nb", "Nc", "ratio", "RatioType"
)

missing_columns <- setdiff(required_columns, names(dat))
if (length(missing_columns) > 0) {
  stop(
    "Required column(s) not found: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Treat empty census-size entries as missing values.
dat$Nc[dat$Nc == ""] <- NA

# Inspect census-size values reported only as lower bounds. These records are
# subsequently excluded because an exact Nb/Nc ratio cannot be calculated.
censored_Nc <- dat %>%
  filter(grepl("^>", Nc)) %>%
  select(
    `Genus/Species`,
    Population,
    Nb,
    Nc,
    ratio,
    RatioType
  )

print(censored_Nc)

# Convert Nb, Nc, and the reported ratio to numeric values. Non-numeric Nc
# entries such as ">5000" become NA and are excluded below.
dat <- dat %>%
  mutate(
    Nb = suppressWarnings(as.numeric(Nb)),
    Nc = suppressWarnings(as.numeric(Nc)),
    ratio_reported = suppressWarnings(as.numeric(ratio))
  )

# ==============================================================================
# 2. Retain records with both Nb and exact Nc estimates
# ==============================================================================

dat_nb <- dat %>%
  filter(
    !is.na(Nb),
    !is.na(Nc),
    Nb > 0,
    Nc > 0
  ) %>%
  mutate(
    Nb_Nc = Nb / Nc
  )

# Number of records included in the analysis.
cat("Number of records with both Nb and exact Nc estimates:", nrow(dat_nb), "\n")

# Inspect the variables used in the analysis.
dat_nb %>%
  select(
    `Genus/Species`,
    Class,
    Population,
    PopulationID,
    Nb,
    Nc,
    Nb_Nc,
    ratio_reported,
    RatioType
  ) %>%
  head() %>%
  print()

# ==============================================================================
# 3. Check agreement with the ratio reported in the Dryad dataset
# ==============================================================================

# Compare the recalculated Nb/Nc ratio with the ratio field provided in Dryad.
# This is a diagnostic check only; all analyses below use Nb_Nc calculated here.
print(summary(dat_nb$Nb_Nc - dat_nb$ratio_reported))

# Inspect the RatioType labels among the retained records. The analysis does not
# filter on RatioType because Nb/Nc is calculated directly from Nb and Nc.
print(table(dat_nb$RatioType, useNA = "ifany"))

# ==============================================================================
# 4. Identify populations with Nb/Nc < 0.1
# ==============================================================================

low_ratio <- dat_nb %>%
  filter(Nb_Nc < 0.1)

cat("Number of populations with Nb/Nc < 0.1:", nrow(low_ratio), "\n")
cat(
  "Percentage of analysed populations with Nb/Nc < 0.1:",
  round(100 * nrow(low_ratio) / nrow(dat_nb), 1), "%\n"
)

# Inspect populations meeting the Nb/Nc < 0.1 criterion.
low_ratio %>%
  select(
    `Genus/Species`,
    Class,
    Population,
    PopulationID,
    Nb,
    Nc,
    Nb_Nc
  ) %>%
  arrange(Nb_Nc)

# ==============================================================================
# 5. Frequency of small absolute Nb among populations with Nb/Nc < 0.1
# ==============================================================================

low_ratio_summary <- low_ratio %>%
  summarise(
    n = n(),
    n_Nb_lt_50 = sum(Nb < 50),
    p_Nb_lt_50 = 100 * mean(Nb < 50),
    n_Nb_lt_100 = sum(Nb < 100),
    p_Nb_lt_100 = 100 * mean(Nb < 100),
    n_Nb_lt_200 = sum(Nb < 200),
    p_Nb_lt_200 = 100 * mean(Nb < 200)
  )

print(low_ratio_summary)

# ==============================================================================
# 6. Sensitivity to the Nb/Nc threshold
# ==============================================================================

ratio_thresholds <- c(0.10, 0.05, 0.02)

summary_table <- bind_rows(
  lapply(ratio_thresholds, function(threshold) {

    tmp <- dat_nb %>%
      filter(Nb_Nc < threshold)

    data.frame(
      NbNc_threshold = paste0("< ", threshold),
      n = nrow(tmp),
      Nb_lt_50_n = sum(tmp$Nb < 50),
      Nb_lt_50_percent = round(100 * mean(tmp$Nb < 50), 1),
      Nb_lt_100_n = sum(tmp$Nb < 100),
      Nb_lt_100_percent = round(100 * mean(tmp$Nb < 100), 1),
      Nb_lt_200_n = sum(tmp$Nb < 200),
      Nb_lt_200_percent = round(100 * mean(tmp$Nb < 200), 1)
    )
  })
)

print(summary_table)

# ==============================================================================
# 7. Taxonomic composition of populations with Nb/Nc < 0.1
# ==============================================================================

# Composition of the low-ratio subset. These percentages describe
# P(taxon | Nb/Nc < 0.1) and should not be interpreted as taxon-specific
# frequencies of low Nb/Nc.
low_ratio_taxonomic_composition <- low_ratio %>%
  count(Class, sort = TRUE) %>%
  mutate(
    percent = 100 * n / sum(n)
  )

print(low_ratio_taxonomic_composition)

# Taxon-specific frequency of Nb/Nc < 0.1 among records available for each
# taxonomic group, i.e. P(Nb/Nc < 0.1 | taxon).
taxon_low_ratio_frequency <- dat_nb %>%
  group_by(Class) %>%
  summarise(
    n_total = n(),
    n_low_ratio = sum(Nb_Nc < 0.1),
    percent_low_ratio = 100 * mean(Nb_Nc < 0.1),
    .groups = "drop"
  ) %>%
  arrange(desc(percent_low_ratio))

print(taxon_low_ratio_frequency)

# ==============================================================================
# 8. Small-Nb frequencies within taxonomic groups
# ==============================================================================

# Number of records within each taxonomic group that simultaneously meet the
# low-ratio criterion and each descriptive Nb threshold.
taxon_summary <- dat_nb %>%
  group_by(Class) %>%
  summarise(
    n_total = n(),
    n_low_ratio = sum(Nb_Nc < 0.1),
    p_low_ratio = 100 * mean(Nb_Nc < 0.1),
    n_low_ratio_Nb_lt_50 = sum(Nb_Nc < 0.1 & Nb < 50),
    n_low_ratio_Nb_lt_100 = sum(Nb_Nc < 0.1 & Nb < 100),
    n_low_ratio_Nb_lt_200 = sum(Nb_Nc < 0.1 & Nb < 200),
    .groups = "drop"
  ) %>%
  arrange(desc(n_low_ratio))

print(taxon_summary)

# ==============================================================================
# 9. Nb distribution within low-ratio populations by taxon
# ==============================================================================

Nb_distribution_by_taxon <- low_ratio %>%
  group_by(Class) %>%
  summarise(
    n = n(),
    median_Nb = median(Nb),
    min_Nb = min(Nb),
    max_Nb = max(Nb),
    p_Nb_lt_50 = 100 * mean(Nb < 50),
    p_Nb_lt_100 = 100 * mean(Nb < 100),
    p_Nb_lt_200 = 100 * mean(Nb < 200),
    .groups = "drop"
  ) %>%
  arrange(median_Nb)

print(Nb_distribution_by_taxon)

# ==============================================================================
# End of script
# ==============================================================================
