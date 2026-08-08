# =============================================================================
# GWIscan - per-family ProtParam distributions (figure + summary in one script)
#
# Reads the compiled GWIscan results table (gwiscan_results.tsv), which already
# carries the family assignment (family) alongside the ProtParam columns, and
# for every family with at least MIN_N members produces:
#   * a faceted boxplot of the selected properties (600 dpi PNG/TIFF + PDF),
#   * a per-family summary table (N, range, mean, SD),
#   * the list of values excluded as outliers.
# The boxplot and the summary are built from the same values (identical Tukey
# outlier rule), so the figure and the table always agree.
#
# Any family reaching MIN_N members is included. Set FAMILY_ORDER to fix the
# top-to-bottom order in the figure.
# =============================================================================

# Required packages, installed on first run if missing.
required <- c("ggplot2", "dplyr", "tidyr", "RColorBrewer")
missing  <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

library(ggplot2)
library(dplyr)
library(tidyr)

setwd(".")

# -- Paths --------------------------------------------------------------------
# Read the compiled table from the working directory (final_results/). The FIGURES
# (the boxplot images) go into a figures/ subfolder of final_results/. The stats and
# outlier TABLES are intermediate working files, so they go to the intermediate
# directory passed as the first command-line argument (falls back to the cwd).
results_file <- "gwiscan_results.tsv"   # GWIscan compiled table (family + ProtParam)

args <- commandArgs(trailingOnly = TRUE)
inter_dir <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "."

fig_dir <- "figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

out_png  <- file.path(fig_dir, "protparam_boxplots.600dpi.png")
out_tiff <- file.path(fig_dir, "protparam_boxplots.600dpi.tiff")
out_pdf  <- file.path(fig_dir, "protparam_boxplots.pdf")
out_stats    <- file.path(inter_dir, "protparam_stats.csv")
out_outliers <- file.path(inter_dir, "protparam_outliers.csv")

# -- Settings -----------------------------------------------------------------
# Grouping column: "family" for families, "superfamily" for the superfamily
# rollup (present only when GWIscan was run in superfamily mode).
FAMILY_COL <- "family"

# A family enters the figure and the table once it reaches this many members.
MIN_N <- 5

# Optional fixed order (top to bottom). Leave empty to order by family size,
# largest at the top; any family reaching MIN_N is included either way.
FAMILY_ORDER <- character(0)

# Column -> panel/label. Edit to change which properties are shown.
# molecularWeight is already reported in kDa by GWIscan.
properties <- c(
  lengthAa        = "Length (aa)",
  molecularWeight = "MW (kDa)",
  isoelectricPoint = "pI",
  gravy           = "GRAVY"
)

# -- Step 1: Load the compiled table ------------------------------------------
read_table_auto <- function(fp) {
  first <- readLines(fp, n = 1)
  sep <- if (grepl("\t", first)) "\t" else ","
  read.delim(fp, header = TRUE, sep = sep, quote = "\"", comment.char = "",
             strip.white = TRUE, stringsAsFactors = FALSE)
}

res <- read_table_auto(results_file)

for (col in c("proteinId", FAMILY_COL, names(properties))) {
  if (!col %in% colnames(res)) {
    stop(sprintf("column '%s' not found in %s", col, results_file))
  }
}

# ProtParam values are per protein; keep one row per protein per family.
dat <- res %>%
  rename(Family = all_of(FAMILY_COL)) %>%
  filter(!is.na(Family), trimws(Family) != "", Family != "-") %>%
  distinct(proteinId, Family, .keep_all = TRUE)

# GWIscan writes "-" where a value could not be computed; treat those as missing.
for (col in names(properties)) {
  if (!is.numeric(dat[[col]])) {
    dat[[col]][trimws(dat[[col]]) %in% c("-", "", "NA", "N/A")] <- NA
    dat[[col]] <- suppressWarnings(as.numeric(dat[[col]]))
  }
}

# -- Step 2: Select families reaching MIN_N members ---------------------------
family_sizes <- dat %>% count(Family, name = "n_members") %>% arrange(desc(n_members))
kept <- family_sizes %>% filter(n_members >= MIN_N) %>% pull(Family)

if (length(kept) == 0) {
  stop(sprintf("no family reaches MIN_N = %d members; nothing to plot", MIN_N))
}

# Order: the given FAMILY_ORDER (restricted to kept families) first, then any
# remaining kept families by descending size.
if (length(FAMILY_ORDER) > 0) {
  ordered_families <- c(intersect(FAMILY_ORDER, kept),
                        setdiff(kept, FAMILY_ORDER))
} else {
  ordered_families <- kept   # already sorted by size
}

dropped <- family_sizes %>% filter(n_members < MIN_N)
if (nrow(dropped) > 0) {
  cat(sprintf("Families below MIN_N = %d (excluded):\n", MIN_N))
  print(as.data.frame(dropped))
}
cat("\nMembers per included family:\n")
print(as.data.frame(family_sizes %>% filter(Family %in% ordered_families)))

dat <- dat %>% filter(Family %in% ordered_families)

# -- Step 3: Long format, one row per gene per property -----------------------
long <- dat %>%
  select(proteinId, Family, all_of(names(properties))) %>%
  pivot_longer(cols = all_of(names(properties)),
               names_to = "Property", values_to = "Value") %>%
  filter(!is.na(Value))

# Flag Tukey outliers (below Q1 - 1.5*IQR or above Q3 + 1.5*IQR), the same
# fences the boxplot whiskers use. Not applied to groups with n < MIN_N.
flagged <- long %>%
  group_by(Family, Property) %>%
  mutate(
    n_group = n(),
    Q1 = quantile(Value, 0.25, names = FALSE),
    Q3 = quantile(Value, 0.75, names = FALSE),
    IQR_val = Q3 - Q1,
    Outlier = n_group >= MIN_N & (Value < Q1 - 1.5 * IQR_val | Value > Q3 + 1.5 * IQR_val)
  ) %>%
  ungroup()

# -- Step 4: Summary statistics on the retained values ------------------------
stats <- flagged %>%
  group_by(Family, Property) %>%
  summarise(
    N_total    = n(),
    N_outliers = sum(Outlier),
    N_used     = sum(!Outlier),
    Min        = min(Value[!Outlier]),
    Max        = max(Value[!Outlier]),
    Mean       = mean(Value[!Outlier]),
    SD         = sd(Value[!Outlier]),
    .groups    = "drop"
  ) %>%
  mutate(
    across(c(Min, Max, Mean, SD), ~ round(.x, 2)),
    Range = paste0(Min, "-", Max),
    Family   = factor(Family, levels = ordered_families),
    Property = factor(Property, levels = names(properties), labels = properties)
  ) %>%
  arrange(Property, Family) %>%
  select(Property, Family, N_total, N_outliers, N_used, Min, Max, Range, Mean, SD)

write.csv(stats, out_stats, row.names = FALSE)

outliers <- flagged %>%
  filter(Outlier) %>%
  mutate(Property = factor(Property, levels = names(properties), labels = properties)) %>%
  select(Property, Family, proteinId, Value) %>%
  arrange(Property, Family, Value)

write.csv(outliers, out_outliers, row.names = FALSE)

# -- Step 5: Faceted boxplot (600 dpi) ----------------------------------------
plot_data <- flagged %>%
  filter(!Outlier) %>%
  mutate(
    Family   = factor(Family, levels = rev(ordered_families)),
    Property = factor(Property, levels = names(properties), labels = properties)
  )

p <- ggplot(plot_data, aes(x = Value, y = Family)) +
  geom_boxplot(aes(fill = Family), outlier.shape = NA, linewidth = 0.35, width = 0.6) +
  geom_jitter(height = 0.18, width = 0, shape = 21, fill = "white",
              color = "grey25", size = 1.1, stroke = 0.25, alpha = 0.8) +
  facet_grid(. ~ Property, scales = "free_x") +
  scale_fill_brewer(palette = "Set3", guide = "none") +
  scale_x_continuous(n.breaks = 4) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5),
    axis.text.x = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.spacing = unit(0.35, "lines")
  )

# Height scales with the number of families so panels stay legible.
fig_h <- max(3.5, 0.55 * length(ordered_families) + 1.5)
ggsave(out_png, p, width = 9, height = fig_h, dpi = 600, device = "png")
ggsave(out_tiff, p, width = 9, height = fig_h, dpi = 600, compression = "lzw", device = "tiff")
ggsave(out_pdf, p, width = 9, height = fig_h)

cat(sprintf("\nWrote:\n  %s\n  %s\n  %s\n  %s (%d rows)\n  %s (%d excluded values)\n",
            out_png, out_tiff, out_pdf, out_stats, nrow(stats),
            out_outliers, nrow(outliers)))
print(as.data.frame(stats))
