setwd("/Users/bfj994/Documents/barcodeMiner/results/")
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(tibble)
library(purrr)

dir_path <- "/Users/bfj994/Documents/barcodeMiner/results/"
lca_files <- list.files(path = dir_path, pattern = "\\.lca\\.gz$", full.names = TRUE)
# Read and combine all files with an extra column indicating the base filename
lca_data <- rbindlist(lapply(lca_files, function(file) {
  df <- fread(file)  # use read.table(file, header = TRUE) if not using data.table
  file_label <- sub("\\.lca\\.gz$", "", basename(file))
  df$sample <- file_label
  return(df)
}))
lca_data <- as.data.table(lca_data)

# Split LCA into taxid, name, and rank
lca_data[, c("taxid", "name", "rank") := tstrsplit(lca, ':"|":"', perl=TRUE)]

# Create the label column for rank + name
lca <- lca_data |>
  mutate(
    # Convert "50k" to 50000 and "1m" to 1000000
    amount = case_when(
      str_detect(sample, "\\d+(?=k)") ~ as.numeric(str_extract(sample, "\\d+(?=k)")) * 1000,
      str_detect(sample, "\\d+(?=m)") ~ as.numeric(str_extract(sample, "\\d+(?=m)")) * 1e6,
      TRUE ~ NA_real_
    ),
    
    # Compute coverage
    coverage = (amount * 45) / 160000,
    
    # Extract DB (method)
    db = str_extract(sample, "(?<=_)[^_]+$"),
    
    # Combine rank and name into a single label
    rank_name = paste0(rank, ": ", name)
  )

lca <- lca |>
  filter(!is.na(name) & !is.na(rank)) |>
  mutate(
    rank_name = paste0(rank, ": ", name)
  )

# Step 2: Simulated total reads per coverage
coverage_reads_lookup <- lca |>
  distinct(coverage, amount) |>
  deframe()  # turns into named vector: coverage -> amount

# Step 3: Consistent color mapping for DBs
db_levels <- unique(lca$db)
db_colors <- setNames(RColorBrewer::brewer.pal(n = length(db_levels), name = "Set2"), db_levels)

# Step 4: Rank counts with filtering
rank_counts <- lca |>
  group_by(rank_name, db, coverage) |>
  summarise(n = n(), .groups = "drop") |>
  filter(n > 1)

# Step 5: Generate plots per coverage
plots_by_coverage <- rank_counts |>
  group_split(coverage) |>
  map(function(df) {
    coverage_val <- unique(df$coverage)
    sim_reads <- coverage_reads_lookup[as.character(coverage_val)]
    
    ggplot(df, aes(x = reorder(rank_name, -n), y = n, fill = db)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = db_colors) +  # Consistent colors
      theme_bw(base_size = 16) +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        panel.spacing = unit(1, "lines")
      ) +
      labs(
        title = paste0(round(coverage_val), "x coverage (Simulated reads: ", format(sim_reads, big.mark = ","), ")"),
        x = "Taxonomic Rank: Name",
        y = "Number of Reads",
        fill = "Method (DB)"
      )
  })

# Save each plot to a file
walk2(plots_by_coverage, seq_along(plots_by_coverage), function(plot, i) {
  ggsave(
    filename = paste0("rank_assignments_", i, ".png"),
    plot = plot,
    width = 12, height = 8, dpi = 300
  )
})


